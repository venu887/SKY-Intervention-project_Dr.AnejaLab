#@@@@@@@@@@@@ 
# Analyses of gene expression patterns showed that stress-responsive transcripts 
## Modified Intent-to-Treat (mITT) Analysis
# Here, we developed a Modified Intent-to-Treat (mITT) analysis algorithm to interpret the clinical data of the SKY and Control groups.
# ==============================================================================
rm(list = ls())
library(limma)
library(tidyverse)
library(org.Hs.eg.db)
library(enrichR)
library(ggrepel)
library(biomaRt)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(ggplot2)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[1]. Differentiall Expression analysis 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#====================
# Step 1: Data Preprocessing of RNAseq TPM Values
#====================
tpm_raw <- read.csv("tpm_allsamples.csv", row.names = 1)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
tpm_df <- tpm_raw %>% 
  rownames_to_column("ensembl_gene_id") %>%
  mutate(ensembl_gene_id_clean = sub("\\..*", "", ensembl_gene_id))

genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
               filters = "ensembl_gene_id",
               values = tpm_df$ensembl_gene_id_clean,
               mart = mart)

# Merge and apply fallback logic (If Symbol is missing, use Ensembl ID)
tpm_mapped <- tpm_df %>%
  left_join(genes, by = c("ensembl_gene_id_clean" = "ensembl_gene_id")) %>%
  mutate(Gene_Symbol = ifelse(is.na(external_gene_name) | external_gene_name == "", 
                              ensembl_gene_id, 
                              external_gene_name)) %>%
  select(ensembl_gene_id, Gene_Symbol, everything(), 
         -ensembl_gene_id_clean, -external_gene_name)

# COLLAPSE DUPLICATES (AVERAGING)
# In this step, we compress the ~60k Ensembl IDs into unique Gene Symbols.
# Multiple IDs for the same Symbol (like Y_RNA) will be averaged.
cat("Original gene count (Ensembl IDs):", nrow(tpm_mapped), "\n")
cat("Extra duplicate rows to be averaged:", sum(duplicated(tpm_mapped$Gene_Symbol)), "\n")

tpm_final <- tpm_mapped %>%
  select(-ensembl_gene_id) %>% 
  group_by(Gene_Symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

cat("Final unique gene count (Symbols):", nrow(tpm_final), "\n")

tpm_final<-as.data.frame(tpm_final)
rownames(tpm_final)<-tpm_final$Gene_Symbol
tpm_final$Gene_Symbol<-NULL

file1<-"tpm_gene_clear_data.csv"
write.csv(tpm_final, file1, row.names = T)





#====================
# Step 2: Intent-to-Treat Analysis Using Linear Mixed Models (limma) using RNAseq data
#====================
load("combined_counts_and_metadata.RData") # Sample Meta data
file1<-"tpm_gene_clear_data.csv"
tpm_raw <-read.csv(file1, row.names = 1)
colnames(tpm_raw)
dim(tpm_raw)

sample_info <- combined_sample_info %>%
  mutate(
    condition = ifelse(str_starts(sample, "Pre"), "Pre", "Post"),
    group = case_when(
      str_detect(sample, "(?i)C_") ~ "control",
      str_detect(sample, "(?i)S")  ~ "sky"
    ),
    SubjectID = str_extract(sample, "(?<=Pre|Post)[0-9]+")
  ) %>%
  filter(!is.na(SubjectID))

# Identify subjects who have both Pre and Post samples
paired_subjects <- sample_info %>%
  group_by(SubjectID) %>%
  filter(n() == 2) %>%
  pull(SubjectID) %>% unique()

# Filter and align metadata
meta <- sample_info %>% 
  filter(SubjectID %in% paired_subjects) %>%
  arrange(SubjectID, condition) 
table(meta$condition, meta$group)
head(meta)
# Create Log2 TPM matrix
expr_log2 <- log2(tpm_raw[, meta$sample] + 1)

# DELTA CALCULATION & STRESS-RESPONSIVE FILTER
# Separate Pre and Post to subtract
meta_pre  <- meta %>% filter(condition == "Pre")
meta_post <- meta %>% filter(condition == "Post")
# Calculate Response: Delta = Post - Pre
expr_delta <- expr_log2[, meta_post$sample] - expr_log2[, meta_pre$sample]
colnames(expr_delta) <- meta_pre$SubjectID
# --- CLEAN NAs FIRST ---
expr_delta_clean <- na.omit(expr_delta)
# --- APPLY NIERATSCHKER FILTER ---
# Difference > 1.2 FC (log2(1.2) ~ 0.263) in at least 20% of subjects
fc_thresh <- log2(1.2)
pass_filter <- rowSums(abs(expr_delta_clean) > fc_thresh) >= (0.20 * ncol(expr_delta_clean))
expr_filtered <- expr_delta_clean[pass_filter, ]
message(paste("Genes analyzed after filtering:", nrow(expr_filtered)))

# LIMMA INTERACTION MODELING
# On Delta data, the Group effect is the Interaction Effect (DiD)
design <- model.matrix(~ group, data = meta_pre)
fit <- lmFit(expr_filtered, design)
fit <- eBayes(fit)

res_df <- topTable(fit, coef = "groupsky", number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Gene")

# Map Gene Symbols
res_df$Symbol <- mapIds(org.Hs.eg.db, keys=res_df$Gene, column="SYMBOL", keytype="ENSEMBL")

# --- Results: Intervention Effect (Interaction) ---
# This is the "groupsky" coefficient from your DiD model
res_intervention <- topTable(fit, coef = "groupsky", number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Gene") %>%
  rename_with(~paste0("Intervention_", .), -Gene) # Prefix columns for clarity

# --- Results: Within-Sky Group Change ---
# We fit a simple model to just the Sky subjects' deltas
sky_subjects <- as.character(meta_pre$SubjectID[meta_pre$group == "sky"])
fit_sky <- eBayes(lmFit(expr_filtered[, sky_subjects]))
res_sky <- topTable(fit_sky, coef = 1, number = Inf) %>%
  rownames_to_column("Gene") %>%
  rename_with(~paste0("Sky_", .), -Gene)

# --- Results: Within-Control Group Change ---
# We fit a simple model to just the Control subjects' deltas
ctrl_subjects <- as.character(meta_pre$SubjectID[meta_pre$group == "control"])
fit_ctrl <- eBayes(lmFit(expr_filtered[, ctrl_subjects]))
res_ctrl <- topTable(fit_ctrl, coef = 1, number = Inf) %>%
  rownames_to_column("Gene") %>%
  rename_with(~paste0("Ctrl_", .), -Gene)

# --- COMBINE EVERYTHING ---
master_results <- res_intervention %>%
  left_join(res_sky, by = "Gene") %>%
  left_join(res_ctrl, by = "Gene") %>%
  mutate(Symbol = mapIds(org.Hs.eg.db, keys=Gene, column="SYMBOL", keytype="ENSEMBL")) %>%
  select(Gene, starts_with("Intervention_"), starts_with("Sky_"), starts_with("Ctrl_"))
head(master_results)
master_results$Symbol<-master_results$Gene
write.csv(master_results, "RCT_Intervension_effect.csv", row.names = T)




#====================
# Step 3: Volcano plot for Intent-to-Treat Analysis Using Linear Mixed Models (limma) 
#====================
master_results <-read.csv("RCT_Intervension_effect.csv", row.names = 1)
head(master_results)
# Top 10 Upregulated (positive logFC, sorted by t-statistic)
top_up <- master_results %>%
  filter(Intervention_logFC > 0, Intervention_P.Value < 0.05) %>%
  arrange(desc(Intervention_t)) %>%
  head(10)
# Top 10 Downregulated (negative logFC, sorted by t-statistic)
top_down <- master_results %>%
  filter(Intervention_logFC < 0, Intervention_P.Value < 0.05) %>%
  arrange(Intervention_t) %>%
  head(10)

top_up[, c("Gene", "Intervention_logFC", "Intervention_t", "Intervention_P.Value")]
top_down[, c("Gene", "Intervention_logFC", "Intervention_t", "Intervention_P.Value")]

fc_thresh <- log2(1.2)
# Prepare data for plotting
res_plot <- master_results %>%
  mutate(
    # Use the new Intervention column names
    logFC_val = Intervention_logFC,
    P_val = Intervention_P.Value,
    
    # Create the labeling logic
    change = case_when(
      P_val < 0.05 & logFC_val > fc_thresh  ~ "Up-regulated Interaction",
      P_val < 0.05 & logFC_val < -fc_thresh ~ "Down-regulated Interaction",
      TRUE ~ "Not Significant"
    ),
    # Create a display label: use Symbol, but if NA, use Gene (Ensembl ID)
    plot_label = ifelse(is.na(Symbol) | Symbol == "" | Symbol == "0", Gene, Symbol)
  )

# Dynamic Counts for Legend
counts <- table(res_plot$change)
up_n   <- ifelse("Up-regulated Interaction" %in% names(counts), counts["Up-regulated Interaction"], 0)
dn_n   <- ifelse("Down-regulated Interaction" %in% names(counts), counts["Down-regulated Interaction"], 0)

res_plot <- res_plot %>%
  mutate(change_labeled = case_when(
    change == "Up-regulated Interaction"   ~ paste0("Up (", up_n, ")"),
    change == "Down-regulated Interaction" ~ paste0("Down (", dn_n, ")"),
    TRUE                                   ~ "Not Significant"
  ))

# Define Colors
color_values <- c("#d73027", "#4575b4", "grey80")
names(color_values) <- c(paste0("Up (", up_n, ")"), paste0("Down (", dn_n, ")"), "Not Significant")

# Create Plot
did_formula <- "Delta*Delta*Expression == (Sky[Post] - Sky[Pre]) - (Control[Post] - Control[Pre])"

volcano_p <- ggplot(res_plot, aes(x = logFC_val, y = -log10(P_val))) +
  geom_point(aes(color = change_labeled), alpha = 0.6, size = 1) +
  # Add labels for top 10 significant genes
  geom_text_repel(data = filter(res_plot, change != "Not Significant") %>% head(10),
                  aes(label = plot_label), size = 2.5, box.padding = 0.5) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "gray30") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray30") +
  scale_color_manual(values = color_values) +
  theme_classic() +
  labs(title = "Effect of SKY trearmene using gene expression",
       x = expression(log[2]~Fold~Change~(Interaction)),
       y = expression(-log[10]~P-value),
       color = "Significance") +
  annotate("text", x = -Inf, y = Inf, label = did_formula,
           parse = TRUE, hjust = -0.05, vjust = 0.5, size = 4, fontface = "italic") +
  theme(legend.position = "top", 
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title.x = element_text(size = 10,color = "black"))

print(volcano_p)

pdf("SKY_volcano_master.pdf", width = 4.5, height = 4.5)
print(volcano_p)
dev.off()


#====================
# Step 5: Box plots of Top genes with highly significant nominal p-value
#====================
genes_to_plot <- c("CD163", "LTC4S", "KCTD17")
present_genes <- genes_to_plot[genes_to_plot %in% rownames(expr_delta)]
print(paste("Genes present:", paste(present_genes, collapse = ", ")))

delta_combined <- data.frame()
for(gene in present_genes) {
  gene_stats <- master_results[master_results$Gene == gene, ]
  
  gene_df <- data.frame(
    SubjectID = colnames(expr_delta),
    Delta = as.numeric(expr_delta[gene, ]),
    Group = meta_pre$group,
    Gene = gene,
    logFC = round(gene_stats$Intervention_logFC, 3),
    p_value = format(gene_stats$Intervention_P.Value, scientific = TRUE, digits = 2),
    t_stat = round(gene_stats$Intervention_t, 2)
  )
  delta_combined <- rbind(delta_combined, gene_df)
}

# --- PLOT 
p2_faceted <- ggplot(delta_combined, aes(x = Group, y = Delta, fill = Group)) +
  geom_boxplot(alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  stat_compare_means(method = "wilcox.test", label.x.npc = "center", size = 3) +
  scale_fill_manual(values = c("control" = "#1f77b4", "sky" = "#ff7f0e")) +
  facet_wrap(~Gene, scales = "free_y") +
  labs(title = "M2 Polarization Signature: CD163, LTC4S, and KCTD17",
       subtitle = "All three M2-associated genes show SKY-specific upregulation",
       x = "", y = "Δ log2(TPM+1) (Post - Pre)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(size = 12, face = "bold"))

print(p2_faceted)









# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2] fgsea: GSEA Rank based pathway analysis 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(fgsea)
library(msigdbr)
library(tidyverse)
library(ggplot2)

master_results <- read.csv("RCT_Intervension_effect.csv", row.names = 1)
head(master_results)
  # Create ranked gene list (using t-statistic)
gene_list <- master_results %>%
    filter(!is.na(Intervention_t), Gene != "") %>%
    distinct(Gene, .keep_all = TRUE) %>%
    arrange(desc(Intervention_t))
  
ranks <- gene_list$Intervention_t
names(ranks) <- gene_list$Gene
length(ranks)
min(ranks)
max(ranks)

# ==============================================================================
# 1. HALLMARK GENESETS (H) - SIMPLE VERSION
# ==============================================================================
hm <- msigdbr(species = "Homo sapiens", category = "H") %>%
    split(x = .$gene_symbol, f = .$gs_name)
  
# Run fgseaMultilevel
set.seed(42)
fgsea_hm <- fgseaMultilevel(
    pathways = hm,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 1e-50,
    nPermSimple = 10000)
  
  # Filter significant (FDR < 0.05)
fgsea_hm_sig <- fgsea_hm %>% 
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
print(fgsea_hm_sig %>% select(pathway, NES, padj, size))
  
  
# ==============================================================================
# 2. REACTOME (C2:CP:REACTOME) - SIMPLE VERSION
# ==============================================================================
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
    split(x = .$gene_symbol, f = .$gs_name)
# Run fgseaMultilevel
set.seed(42)
fgsea_reactome <- fgseaMultilevel(
    pathways = reactome,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 1e-50,
    nPermSimple = 10000)
  
  # Filter significant (FDR < 0.05)
fgsea_reactome_sig <- fgsea_reactome %>% 
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
print(fgsea_reactome_sig %>% select(pathway, NES, padj, size) %>% head(10))
  
  
# ==============================================================================
# 3. GO:BP (C5:GO:BP) - SIMPLE VERSION
# ==============================================================================
gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
    split(x = .$gene_symbol, f = .$gs_name)
  
# Run fgseaMultilevel
set.seed(42)
fgsea_gobp <- fgseaMultilevel(
    pathways = gobp,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 1e-50,
    nPermSimple = 10000)
  
# Filter significant (FDR < 0.05)
fgsea_gobp_sig <- fgsea_gobp %>% 
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
print(fgsea_gobp_sig %>% select(pathway, NES, padj, size) %>% head(10))
  
  
# ==============================================================================
# 4. GO:MF (C5:GO:MF) - SIMPLE VERSION
# ==============================================================================
gomf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
    split(x = .$gene_symbol, f = .$gs_name)
  
# Run fgseaMultilevel
set.seed(42)
fgsea_gomf <- fgseaMultilevel(
    pathways = gomf,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 1e-50,
    nPermSimple = 10000)
  
# Filter significant (FDR < 0.05)
fgsea_gomf_sig <- fgsea_gomf %>% 
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
print(fgsea_gomf_sig %>% select(pathway, NES, padj, size) %>% head(10))
  
  
# ==============================================================================
# 5. GO:CC (C5:GO:CC) - SIMPLE VERSION
# ==============================================================================
gocc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>%
    split(x = .$gene_symbol, f = .$gs_name)
  
# Run fgseaMultilevel
set.seed(42)
fgsea_gocc <- fgseaMultilevel(
    pathways = gocc,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 1e-50,
    nPermSimple = 10000)
  
  # Filter significant (FDR < 0.05)
fgsea_gocc_sig <- fgsea_gocc %>% 
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
print(fgsea_gocc_sig %>% select(pathway, NES, padj, size) %>% head(10))
  
  
# ==============================================================================
# 6. KEGG (C2:CP:KEGG) - SIMPLE VERSION
# ==============================================================================
kegg<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY") %>%
    split(x = .$gene_symbol, f = .$gs_name)
  # Run fgseaMultilevel
set.seed(42)
fgsea_kegg <- fgseaMultilevel(
    pathways = kegg,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 1e-50,
    nPermSimple = 10000)
  
# Filter significant (FDR < 0.05)
fgsea_kegg_sig <- fgsea_kegg %>% 
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
print(fgsea_kegg_sig %>% select(pathway, NES, padj, size) %>% head(10))
  
  
# ==============================================================================
# 7. ALL COMBINED (for summary table) and save the results
# ==============================================================================
all_results_complete <- bind_rows(
    fgsea_hm_sig %>% mutate(Database = "Hallmark"),
    fgsea_reactome_sig %>% mutate(Database = "Reactome"),
    fgsea_gobp_sig %>% mutate(Database = "GO:BP"),
    fgsea_gomf_sig %>% mutate(Database = "GO:MF"),
    fgsea_gocc_sig %>% mutate(Database = "GO:CC"),
    fgsea_kegg_sig %>% mutate(Database = "KEGG")
) %>%
# Keep ALL original fgsea columns plus new Database column
select(Database, everything()) %>%
arrange(Database, padj)
  

library(openxlsx)
# Create directory
dir.create("GSEA_Results", showWarnings = FALSE, recursive = TRUE)
  
# Function to convert list columns to text
fix_df <- function(df) {
    as.data.frame(lapply(df, function(x) {
      if(is.list(x)) sapply(x, function(y) paste(unlist(y), collapse = "; ")) else x
    }), stringsAsFactors = FALSE)}
  
# Combine all significant results with Database column
all_results <- bind_rows(
    fgsea_hm_sig %>% mutate(Database = "Hallmark"),
    fgsea_reactome_sig %>% mutate(Database = "Reactome"),
    fgsea_gobp_sig %>% mutate(Database = "GO:BP"),
    fgsea_gomf_sig %>% mutate(Database = "GO:MF"),
    fgsea_gocc_sig %>% mutate(Database = "GO:CC"),
    fgsea_kegg_sig %>% mutate(Database = "KEGG")) %>% 
    select(Database, everything()) %>%
    arrange(Database, padj)
  
# Create workbook and add sheet
wb <- createWorkbook()
addWorksheet(wb, "Significant_Pathways")
writeData(wb, "Significant_Pathways", fix_df(all_results))

saveWorkbook(wb, "fgseaMulti_SKYvsControl_Significant.xlsx", 
               overwrite = TRUE)
  
  
  



#@@@@@@@@@@@@@@@@@@@@@@@@@
# Extract the significant pathways and the pathways associated to srress and inflamation
#@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls())
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(stringr)
  
pathways <- read_xlsx("fgseaMulti_SKYvsControl_Significant.xlsx")

master_results <- read.csv("RCT_Intervension_effect.csv", row.names = 1)
head(master_results)
  
# Exclude, if pathways anything associated to other than cancer, virus, mouse, and other unwanted pathways 
# Because we are observing changes in normal individual who practiced SKY 
exclude_terms <- c(
    # Cancer-related
    "CANCER", "TUMOR", "TUMOUR", "ONCOGENE", "CARCINOMA", "METASTASIS", "MALIGNANT",
    "P53", "TP53", "MYC", "RAS", "SRC", "EGFR", "ERBB", "NEOPLASM", "LEUKEMIA",
    "GLIOMA", "MELANOMA", "SARCOMA", "LYMPHOMA", "BLASTOMA",
    # Cell cycle/proliferation (cancer-associated)
    "CELL_CYCLE", "MITOTIC", "E2F", "G2M", "M_PHASE", "S_PHASE",
    "DNA_REPLICATION", "CHROMOSOME_SEGREGATION", "SISTER_CHROMATID",
    "KINETOCHORE", "SPINDLE", "CYTOKINESIS", "CELL_DIVISION", "CENTROSOME",
    "CONDENSED_CHROMOSOME", "CHROMOSOMAL_REGION",
    # Virus/infection-related
    "VIRUS", "VIRAL", "INFLUENZA", "SARS", "COVID", "CORONA", "HIV", "HPV",
    "HEPATITIS", "EBV", "CYTOMEGALO", "HERPES", "ADENOVIRUS", "RETROVIRUS",
    "PARASITE", "INFECTION", "PATHOGEN", "BACTERIA", "BACTERIAL",
    # Mouse/model organism specific
    "MOUSE", "MURINE", "RAT", "MUS_MUSCULUS", "RODENT",
    # Development/morphogenesis (often not relevant for stress response)
    "EMBRYO", "MORPHOGENESIS", "ORGANOGENESIS", "HEART_DEVELOPMENT",
    "VENTRICULAR", "CARDIAC", "MUSCLE_TISSUE", "FOREBRAIN", "NEURON",
    "AXON", "DENDRITE", "SYNAPSE", "NEUROTRANSMITTER", "NERVOUS_SYSTEM",
    "CILIUM", "FLAGELLUM",
    # Translation/ribosome (housekeeping)
    "RIBOSOME", "RIBOSOMAL", "TRANSLATION", "RRNA", "TRNA", "PROTEASOME",
    # General metabolic (unless specifically stress-related)
    "GLYCOLYSIS", "GLUCONEOGENESIS", "FATTY_ACID", "LIPID_METABOLISM",
    "AMINO_ACID", "CARBOHYDRATE", "CHOLESTEROL", "STEROID",
    # Other unnecessary terms
    "KERATINIZATION", "CORNIFICATION", "EPITHELIAL", "SECRETION",
    "VESICLE", "GOLGI", "ENDOPLASMIC_RETICULUM", "ENDOSOME",
    "EXOCYTOSIS", "ENDOCYTOSIS", "TRANSPORT")
  
  # Terms associated stress and inflammation related terms to INCLUDE
  include_terms <- c(
    # Stress responses
    "STRESS", "OXIDATIVE", "HYPOXIA", "UV_RESPONSE", "DNA_REPAIR",
    "STARVATION", "NUTRIENT", "HORMONE", "INSULIN", "GROWTH_FACTOR",
    "APOPTOSIS", "NECROPTOSIS", "AUTOPHAGY", "MITOPHAGY", "UNFOLDED_PROTEIN",
    "HEAT_SHOCK", "OSMOTIC", "GENOTOXIC", "XENOBIOTIC",
    # Inflammation
    "INFLAMMATORY", "INFLAMMATION", "ACUTE", "CHRONIC", "COMPLEMENT",
    "COAGULATION", "WOUND_HEALING",
    # Cytokine/Chemokine
    "CYTOKINE", "CHEMOKINE", "INTERLEUKIN", "INTERFERON", "TNFA", "NFKB",
    "TGF_BETA", "JAK_STAT", "IL6", "IL1", "IL10", "IL2",
    # Immune response
    "IMMUNE", "T_CELL", "B_CELL", "MACROPHAGE", "NEUTROPHIL", "MONOCYTE",
    "LEUKOCYTE", "MYELOID", "LYMPHOCYTE", "DENDRITIC", "MAST_CELL",
    "PATTERN_RECOGNITION", "TOLL_LIKE", "NOD_LIKE", "RIG_I",
    "ANTIGEN", "MHC", "BCR", "TCR", "FC_RECEPTOR",
    "ANTIVIRAL", "ISG15", "IFN",
    # Cell signaling (stress-related)
    "MTORC1", "TOR", "MAPK", "P38", "JNK", "ERK", "AKT", "PI3K",
    "NF_KAPPA_B", "SIGNALING",
    # Differentiation (immune-related)
    "DIFFERENTIATION", "ACTIVATION", "POLARIZATION",
    # Metabolism (stress-related)
    "RESPIRATION", "ELECTRON_TRANSPORT", "ATP", "MITOCHONDRIAL",
    "OXIDATIVE_PHOSPHORYLATION")
  
all_pathways <- pathways$pathway
# First, exclude unwanted pathways
filtered_pathways <- all_pathways[!str_detect(all_pathways, paste(exclude_terms, collapse = "|"))]
  
stress_inflammation_pathways <- filtered_pathways[str_detect(filtered_pathways, paste(include_terms, collapse = "|"))]

# finally we got these 43 significant pathways, here are the lise
target_pathways <- c(
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_UV_RESPONSE_DN",
    "GOBP_RESPONSE_TO_STRESS",
    "GOBP_CELLULAR_RESPONSE_TO_STRESS",
    "GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS",
    "GOBP_RESPONSE_TO_STARVATION",
    "GOBP_CELLULAR_RESPONSE_TO_STARVATION",
    "GOBP_RESPONSE_TO_HORMONE",
    "GOBP_RESPONSE_TO_PEPTIDE_HORMONE",
    "GOBP_RESPONSE_TO_INSULIN",
    "GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS",
    "GOBP_CELLULAR_RESPONSE_TO_HORMONE_STIMULUS",
    "GOBP_CELLULAR_RESPONSE_TO_NUTRIENT_LEVELS",
    "GOBP_DNA_REPAIR",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_HYPOXIA",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    # Inflammation
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_COAGULATION",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "GOBP_INFLAMMATORY_RESPONSE",
    "GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE",
    "GOBP_ACUTE_INFLAMMATORY_RESPONSE",
    # Cytokine signaling
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "GOBP_REGULATION_OF_RECEPTOR_SIGNALING_PATHWAY_VIA_JAK_STAT",
    "REACTOME_FCERI_MEDIATED_MAPK_ACTIVATION",
    "REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS",
    # Immune cell function
    "GOBP_MYELOID_CELL_DIFFERENTIATION",
    "GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE",
    "GOBP_LEUKOCYTE_ACTIVATION",
    "GOBP_LEUKOCYTE_MIGRATION",
    "GOBP_LEUKOCYTE_CHEMOTAXIS",
    "GOBP_MACROPHAGE_ACTIVATION",
    "GOBP_NEUTROPHIL_ACTIVATION",
    "GOBP_NEUTROPHIL_CHEMOTAXIS",
    "GOBP_POSITIVE_REGULATION_OF_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_REGULATION_OF_TOLL_LIKE_RECEPTOR_3_SIGNALING_PATHWAY",
    "GOBP_POSITIVE_REGULATION_OF_TOLL_LIKE_RECEPTOR_4_SIGNALING_PATHWAY",
    "GOBP_ENDOLYSOSOMAL_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
    # Antigen presentation
    "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
    "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS",
    "REACTOME_CD22_MEDIATED_BCR_REGULATION",
    # Antiviral responses (if relevant to your study)
    "REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES",
    "REACTOME_ISG15_ANTIVIRAL_MECHANISM",
    # TGF-beta signaling (stress-related)
    "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION",
    "GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA",
    "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
    "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX",
    "REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS",
    # Metabolic stress - KEEP ORIGINAL NAMES WITH DATABASE PREFIXES
    "KEGG_OXIDATIVE_PHOSPHORYLATION",
    "GOBP_OXIDATIVE_PHOSPHORYLATION",
    "GOBP_ATP_BIOSYNTHETIC_PROCESS",
    "GOBP_AEROBIC_RESPIRATION",
    "GOBP_PROTON_TRANSMEMBRANE_TRANSPORT",
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
    # Cell signaling (stress-related)
    "HALLMARK_MTORC1_SIGNALING",
    "GOBP_TOR_SIGNALING",
    "GOBP_NEGATIVE_REGULATION_OF_TOR_SIGNALING",
    "GOBP_NEGATIVE_REGULATION_OF_TORC1_SIGNALING"
  )

target_pathways <- target_pathways[target_pathways %in% all_pathways]
plot_data_comp <- pathways %>%
    filter(pathway %in% target_pathways) %>% #  padj < 0.02
    rowwise() %>%
    mutate(
      Database = case_when(
        str_detect(pathway, "HALLMARK") ~ "Hallmark",
        str_detect(pathway, "GOBP") ~ "GO:BP",
        str_detect(pathway, "KEGG") ~ "KEGG",
        TRUE ~ "Reactome"
      ),
      Short_Name = pathway,
      Leading_Count = length(unlist(strsplit(as.character(leadingEdge), "; "))),
      Up = ifelse(NES > 0, Leading_Count, 0),
      Down = ifelse(NES < 0, Leading_Count, 0),
      Other = size - (Up + Down),
      Significance = case_when(padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ "")
    ) %>%
    ungroup() %>%
    filter(!str_detect(Short_Name, "KEGG_OXIDATIVE_PHOSPHORYLATION|GOBP_OXIDATIVE_PHOSPHORYLATION")) %>%
    # Reshape and Fix Factor Levels for correct stacking and color mapping
    pivot_longer(cols = c(Up, Down, Other), names_to = "Status", values_to = "Count") %>%
    mutate(
      Percentage = (Count / size) * 100,
      Status = factor(Status, levels = c("Other", "Down", "Up")), # 'Other' at the base
      # Text color: White for Up/Down bars, Black for the gray 'Other' bar
      Text_Col = ifelse(Status == "Other", "black", "white"),
      Plot_Label = reorder(paste0(Database, ": ", Short_Name), NES)
    )
  
# --- Generate Custom Plot ---
p_final <- ggplot(plot_data_comp, aes(x = Plot_Label, y = Percentage, fill = Status)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.1) +
    
    # Corrected Labels (Percentage + n)
    geom_text(aes(label = ifelse(Count > 0, paste0(round(Percentage, 0), "%\n(n=", Count, ")"), ""),
                  color = Text_Col), 
              position = position_stack(vjust = 0.5), 
              size = 2.0, fontface = "bold", lineheight = 0.8) +
    
    # Map text colors explicitly
    scale_color_identity() + 
    
    # External Info Labels (n, NES, Stars)
    geom_text(data = plot_data_comp %>% distinct(Plot_Label, size, NES, Significance),
              aes(x = Plot_Label, y = 102, 
                  label = paste0("n=", size, " | NES=", round(NES, 2), Significance)), 
              inherit.aes = FALSE, hjust = 0, size = 2.4, fontface = "italic", color = "grey30") +
    
    coord_flip() +
    theme_bw() +
    scale_fill_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "Other" = "grey90"),
                      labels = c("Other Genes",  "Down in SKY" ,"Up in SKY")) +
    scale_y_continuous(limits = c(0, 135), breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%")) +
    labs(title = "SKY treatment effect GSEA",
         x = NULL, y = "Proportion of Pathway Genes (%)", fill = "Regulation") +
    theme(
      axis.text.y = element_text(size = 7, face = "bold"),
      legend.position = "top",
      # --- Grid Removal ---
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
    )
  
print(p_final)

pdf("GSEA_Final.pdf", width = 10.5, height = 10)
print(p_final)
dev.off()





#@@@@@@@@@@ Save the results as csv
pathways_filtered <- pathways %>%
  filter(pathway %in% target_pathways) %>%
  # Add a cleaned name column for readability
  mutate(
    Clean_Name = case_when(
      str_detect(pathway, "OXIDATIVE_PHOSPHORYLATION") ~ {
        db <- case_when(
          str_detect(pathway, "^HALLMARK") ~ "HALLMARK",
          str_detect(pathway, "^KEGG") ~ "KEGG",
          str_detect(pathway, "^GOBP") ~ "GOBP",
          str_detect(pathway, "^REACTOME") ~ "REACTOME",
          TRUE ~ ""
        )
        paste0(db, ": Oxidative Phosphorylation")
      },
      TRUE ~ {
        cleaned <- str_replace_all(pathway, "HALLMARK_|GOBP_|KEGG_|REACTOME_", "")
        cleaned <- str_replace_all(cleaned, "_", " ")
        str_to_title(cleaned)
      }
    ),
    # Add direction based on NES
    Direction = ifelse(NES > 0, "Up-regulated", "Down-regulated"),
    # Add significance category
    Significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  # Select and reorder columns for better readability
  select(
    Pathway_ID = pathway,
    Clean_Name,
    Direction,
    NES,
    padj,
    Significance,
    size,
    leadingEdge,
    everything())

write.csv(pathways_filtered, 
          "/Users/mekalav/Documents/UAB/SKY/RNAseq_res/RCT_Limma/SKY_Stress_Inflammation_Pathways.csv", 
          row.names = T)
  

# ==============================================================================
# HALLMARK GSEA PLOT — using plot_data_comp1: Exclude only Hallmark pathways
# ==============================================================================
UP_COL   <- "#d73027"
DOWN_COL <- "#4575b4"

gene_status <- master_results %>%
    select(Gene, Intervention_logFC) %>%
    mutate(Status = case_when(
      Intervention_logFC > 0 ~ "Up",
      Intervention_logFC < 0 ~ "Down",
      TRUE                   ~ "Other"
    ))
  
plot_hallmark <- plot_data_comp1 %>%
    mutate(
      Significance = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE         ~ ""),
      Short_Name = reorder(pathway, NES)
    ) %>%
    rowwise() %>%
    mutate(probs = list({
      core_genes <- unlist(strsplit(leadingEdge, ";\\s*"))
      path_data  <- gene_status %>% filter(Gene %in% core_genes)
      up_c <- sum(path_data$Status == "Up")
      dw_c <- sum(path_data$Status == "Down")
      ot_c <- size - (up_c + dw_c)
      data.frame(Up = up_c, Down = dw_c, Other = ot_c)
    })) %>%
    unnest(probs) %>%
    pivot_longer(
      cols      = c(Up, Down, Other),
      names_to  = "Gene_Status",
      values_to = "Gene_Count") %>%
    mutate(
      Percentage  = (Gene_Count / size) * 100,
      Gene_Status = factor(Gene_Status, levels = c("Other", "Down", "Up")),
      Text_Col    = ifelse(Gene_Status == "Other", "black", "white")
    )
  
# ===============
#  PLOT
# ===============
p_hallmark <- ggplot(plot_hallmark,
                       aes(x = Short_Name,
                           y = Percentage,
                           fill = Gene_Status)) +
    geom_bar(stat = "identity", color = "white", linewidth = 0.1) +
    # Inside bar labels
    geom_text(aes(label = ifelse(Gene_Count > 0,
                         paste0(round(Percentage, 0), "%\n(n=", Gene_Count, ")"), ""),
          color = Text_Col),
      position   = position_stack(vjust = 0.5),
      size       = 2.0,
      fontface   = "bold",
      lineheight = 0.8) +
    scale_color_identity() +
    # External NES + significance + n labels
    geom_text(
      data = plot_hallmark %>%
        distinct(Short_Name, size, NES, Significance, padj),
      aes(x     = Short_Name,
          y     = 102,
          label = paste0("n=", size,
                         " | NES=", round(NES, 2),
                         " | FDR=", formatC(padj, format = "e", digits = 1),
                         Significance)),
      inherit.aes = FALSE,
      hjust    = 0,
      size     = 2.2,
      fontface = "italic",
      color    = "grey30") +
    # Vertical reference line at 100%
    geom_hline(yintercept = 100,
               linetype   = "dashed",
               color      = "grey60",
               linewidth  = 0.3) +
    coord_flip() +
    scale_fill_manual(
      values = c("Up"    = UP_COL,
                 "Down"  = DOWN_COL,
                 "Other" = "grey90"),
      labels = c("Other Genes",
                 "Down in SKY",
                 "Up in SKY")) +
    scale_y_continuous(
      limits = c(0, 145),
      breaks = seq(0, 100, 25),
      labels = function(x) paste0(x, "%")) +
    labs(title    = "GSEA SKY Treatment Effect",
      x        = NULL,
      y        = "Proportion of Hallmark Pathway Genes (%)",
      fill     = "Regulation") +
    
    theme_bw(base_size = 9) +
    theme(
      text             = element_text(family = "Helvetica"),
      plot.title       = element_text(size   = 9,   face = "bold",   hjust = 0),
      plot.subtitle    = element_text(size   = 7.5, face = "italic", hjust = 0, color  = "grey40"),
      axis.text.y      = element_text(size   = 7.5, face = "bold"),
      axis.text.x      = element_text(size   = 7),
      axis.title.x     = element_text(size   = 7.5),
      legend.position  = "top",
      legend.title     = element_text(size   = 7.5, face = "bold"),
      legend.text      = element_text(size   = 7),
      legend.key.size  = unit(3, "mm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "black",
                                      fill   = NA,
                                      linewidth = 0.5),
      plot.margin      = margin(5, 25, 5, 5, "mm"))
print(p_hallmark)

pdf(file.path(out, "GSEA_Hallmark_Final.pdf"),
      width = 8.5, height = 4.5, useDingbats = FALSE)
  print(p_hallmark)
dev.off()

  
  
  
  
  
  
  
  
  













# ==============================================================================
# Custom Heatmap: using the top Pathways from Hallmark pathway analysis
# ==============================================================================
rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

file1 <- "tpm_gene_clear_data.csv"
load("combined_counts_and_metadata.RData")
tpm_raw  <- read.csv(file1, row.names = 1)
pathways <- read.csv("SKY_Stress_Inflammation_Pathways.csv",
                       row.names = 1)
# ==============================================================================
# METADATA
# ==============================================================================
meta <- combined_sample_info %>%
    mutate(
      condition  = factor(ifelse(str_starts(sample, "Pre"), "Pre", "Post"),
                          levels = c("Pre", "Post")),
      group      = factor(case_when(
        str_detect(sample, "(?i)C_") ~ "Control",
        str_detect(sample, "(?i)S")  ~ "SKY"),
        levels = c("Control", "SKY")),
      SubjectID  = str_extract(sample, "(?<=Pre|Post)[0-9]+"),
      group_time = factor(paste(group, condition, sep = "_"),
                          levels = c("Control_Pre", "Control_Post",
                                     "SKY_Pre",     "SKY_Post"))
    ) %>% filter(!is.na(SubjectID)) %>%
    group_by(SubjectID) %>% filter(n() == 2) %>% ungroup() %>%
    arrange(group, SubjectID, condition)
  
expr_log2 <- log2(tpm_raw[, meta$sample] + 1)
  
# ==============================================================================
# DESIGN TOKENS
# ==============================================================================
NC_FONT      <- "Helvetica"
NC_BODY      <- 6
NC_LABEL     <- 7
NC_TITLE     <- 8
NC_PANEL     <- 9
NC_MAIN      <- 10
NC_LEG_TITLE <- 8
NC_LEG_BODY  <- 7
  
# Original colors restored
col_z <- colorRamp2(c(-2, 0, 2), c("#3A86FF", "white", "#FF006E"))
  
pal_anno <- c(
    Control_Pre  = "#E8D5C4",
    Control_Post = "#9C6644",
    SKY_Pre      = "#CAF0F8",
    SKY_Post     = "#0077B6")
  
# Plain ASCII — no Unicode encoding issues
panel_info <- list(
    a = list(id    = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
             title = "Oxidative Phosphorylation",
             col   = "#0077B6",
             nes   = "NES = +1.91, FDR = 0.001",
             file  = "Heatmap_A_OxPhos.pdf"),
    b = list(id    = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
             title = "Interferon-gamma Response",
             col   = "#0077B6",
             nes   = "NES = -1.64, FDR = 0.016",
             file  = "Heatmap_B_IFNg.pdf"),
    c = list(id    = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
             title = "TNF-alpha Signalling via NF-kB",
             col   = "#D62828",
             nes   = "NES = -1.76, FDR = 0.003",
             file  = "Heatmap_C_NFKB.pdf"),
    d = list(id    = "HALLMARK_MTORC1_SIGNALING",
             title = "mTORC1 Signalling",
             col   = "#D62828",
             nes   = "NES = -1.84, FDR = 0.001",
             file  = "Heatmap_D_mTORC1.pdf"))
  
# ==============================================================================
# TOP ANNOTATION
# ==============================================================================
ha_top <- HeatmapAnnotation(
    Group = meta$group_time,
    col   = list(Group = pal_anno),
    simple_anno_size   = unit(3, "mm"),
    border             = TRUE,
    annotation_name_gp = gpar(fontsize   = NC_LABEL,
                              fontfamily = NC_FONT,
                              fontface   = "plain"),
    show_legend = FALSE)
  
# ==============================================================================
# HEATMAP FACTORY
# ==============================================================================
make_ht <- function(info) {
    genes <- unlist(strsplit(
      pathways %>% filter(Pathway_ID == info$id) %>% pull(leadingEdge), ";\\s*"))
    genes <- intersect(genes, rownames(expr_log2))
    mat   <- t(scale(t(expr_log2[genes, ])))
    colnames(mat) <- meta$SubjectID
    
    Heatmap(mat,
            name                = info$title,
            col                 = col_z,
            top_annotation      = ha_top,
            column_split        = meta$group_time,
            cluster_columns     = FALSE,
            cluster_rows        = TRUE,
            show_row_dend       = FALSE,
            show_column_dend    = FALSE,
            show_column_names   = TRUE,
            column_names_gp     = gpar(fontsize   = NC_BODY,
                                       fontfamily = NC_FONT),
            column_names_rot    = 90,
            row_names_gp        = gpar(fontsize   = NC_BODY,
                                       fontfamily = NC_FONT,
                                       fontface   = "italic"),
            row_names_side      = "left",
            column_gap          = unit(c(0.8, 3, 0.8), "mm"),
            row_gap             = unit(0.5, "mm"),
            border              = TRUE,
            border_gp           = gpar(col = "grey30", lwd = 0.5),
            rect_gp             = gpar(col = "white",  lwd = 0.2),
            column_title        = info$nes,
            column_title_gp     = gpar(fontsize   = NC_TITLE,
                                       fontfamily = NC_FONT,
                                       fontface   = "italic",
                                       col        = "grey40"),
            show_heatmap_legend = FALSE)}
  
# ==============================================================================
# BUILD ALL FOUR HEATMAPS
# ==============================================================================
hts <- lapply(panel_info, make_ht)
  
# ==============================================================================
# LEGENDS — original colors
# ==============================================================================
lgd_zscore <- Legend(
    col_fun      = col_z,
    title        = "Row z-score",
    title_gp     = gpar(fontsize   = NC_LEG_TITLE,
                        fontfamily = NC_FONT,
                        fontface   = "bold"),
    labels_gp    = gpar(fontsize   = NC_LEG_BODY,
                        fontfamily = NC_FONT),
    legend_width = unit(3, "cm"),
    grid_height  = unit(3, "mm"),
    direction    = "horizontal",
    at           = c(-2, -1, 0, 1, 2),
    labels       = c("-2", "-1", "0", "1", "2")
  )
  
  lgd_group <- Legend(
    labels      = c("Control - Pre", "Control - Post",
                    "SKY - Pre",     "SKY - Post"),
    title       = "Group",
    title_gp    = gpar(fontsize   = NC_LEG_TITLE,
                       fontfamily = NC_FONT,
                       fontface   = "bold"),
    labels_gp   = gpar(fontsize   = NC_LEG_BODY,
                       fontfamily = NC_FONT),
    legend_gp   = gpar(fill = c("#E8D5C4", "#9C6644",
                                "#CAF0F8", "#0077B6")),
    grid_height = unit(3, "mm"),
    grid_width  = unit(3, "mm"),
    direction   = "horizontal"
  )
  
  packed_legend <- packLegend(lgd_zscore, lgd_group,
                              direction = "horizontal",
                              gap       = unit(15, "mm"))
  
out <- "RNAseq_res/RCT_Limma/HM"
  draw_single_page <- function(label, info, ht) {
    pushViewport(viewport(
      layout = grid.layout(
        nrow    = 3,
        ncol    = 1,
        heights = unit(c(0.40, 1.00, 9.0), "inches")  # reduced from 11.0 to 9.0
      )
    ))
    # Panel title + NES
    pushViewport(viewport(layout.pos.row = 1))
    grid.text(
      paste0(label, "   ", info$title, "   |   ", info$nes),
      x  = unit(0.02, "npc"), hjust = 0,
      gp = gpar(fontsize   = NC_PANEL,
                fontfamily = NC_FONT,
                fontface   = "bold",
                col        = info$col))
    popViewport()
    
    # Legend
    pushViewport(viewport(layout.pos.row = 2,
                          x = 0.5, y = 0.5, just = "centre"))
    draw(packed_legend)
    popViewport()
    
    # Heatmap
    pushViewport(viewport(layout.pos.row = 3))
    draw(ht,
         show_heatmap_legend    = FALSE,
         show_annotation_legend = FALSE,
         newpage                = FALSE,
         padding                = unit(c(3, 3, 3, 3), "mm"))
    popViewport()
    
    popViewport()
  }
  
  draw_label <- function(row, col, label, info) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.text(
      paste0(label, "   ", info$title),
      x  = unit(0.015, "npc"), hjust = 0,
      gp = gpar(fontsize   = NC_PANEL,
                fontfamily = NC_FONT,
                fontface   = "bold",
                col        = info$col))
    popViewport()
  }
  
  draw_ht <- function(row, col, ht) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    draw(ht,
         show_heatmap_legend    = FALSE,
         show_annotation_legend = FALSE,
         newpage                = FALSE,
         padding                = unit(c(3, 3, 3, 3), "mm"))
    popViewport()
  }
  
  
# ==============================================================================
# COMBINED 2x2 PDF — reduce heatmap row heights
# ==============================================================================
pdf(file.path(out, "Heatmap_Combined_2x2.pdf"),
      width = 18, height = 18, useDingbats = FALSE)  # reduced from 20 to 18
  pushViewport(viewport(
    layout = grid.layout(
      nrow    = 6,
      ncol    = 2,
      heights = unit(c(0.40,   # main title
                       1.00,   # legend
                       0.35,   # panel labels row 1
                       8.0,    # heatmaps row 1  ← reduced from 10.5
                       0.35,   # panel labels row 2
                       8.0),   # heatmaps row 2  ← reduced from 10.5
                     "inches"),
      widths  = unit(c(1, 1), "null"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
  grid.text( "SKY Breathing Intervention - Pathway-Level Gene Expression (ITT DiD Analysis)",
    gp = gpar(fontsize   = NC_MAIN,
              fontfamily = NC_FONT,
              fontface   = "bold",
              col        = "black"))
  popViewport()
# Legend
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:2,
                        x = 0.5, y = 0.5, just = "centre"))
draw(packed_legend)
popViewport()
  
# Draw all panels
draw_label(3, 1, "a", panel_info$a)
draw_label(3, 2, "b", panel_info$b)
draw_ht(4, 1, hts$a)
draw_ht(4, 2, hts$b)
draw_label(5, 1, "c", panel_info$c)
draw_label(5, 2, "d", panel_info$d)
draw_ht(6, 1, hts$c)
draw_ht(6, 2, hts$d)
  
popViewport()
dev.off()
  

  
  
  
  
  
# ==============================================================================
# Heatmap of Important up and down regulated genes in pathways
# ==============================================================================
rm(list = ls())
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

file1 <- "tpm_gene_clear_data.csv"
load("combined_counts_and_metadata.RData")
tpm_raw        <- read.csv(file1, row.names = 1)
master_results <- read.csv("RCT_Intervension_effect.csv", row.names = 1)
  
gene_groups <- list(
    "Complex I (NDUF)"     = c("NDUFA1", "NDUFA2", "NDUFA4", "NDUFB6", "MT-ND6"),
    "Complex III/IV"       = c("UQCRC1", "UQCRQ", "UQCR11", "COX4I1", "COX5B", "COX6A1",
                               "COX6C", "COX7A2", "COX7B", "COX7C", "COX8A"),
    "ATP Synthase"         = c("ATP5F1D", "ATP5F1E", "ATP5ME", "ATP5PD", "ATP5PF", "ATP5MC1"),
    "Mito Quality/Markers" = c("PINK1", "TAFAZZIN", "SURF1", "CD163"),
    "NF-kB/Inflam"         = c("TNFAIP3", "NFKBIA", "NFKBIZ", "BIRC2", "LITAF", "REL"),
    "Chemokines"           = c("CCL2", "CXCL2", "CXCL3", "PTGS2", "IL6ST"),
    "Stress/MAPK"          = c("FOSB", "EGR2", "EGR3", "NR4A2", "MAP3K8", "MAP2K3"),
    "mTOR/Glycolysis"      = c("HIF1A", "HK2", "LDHA", "SLC2A1", "PIK3CA", "PIK3R1",
                               "RPS6KB1", "MTOR"),
    "Cell Cycle"           = c("CCNA2", "AURKA", "MCM2", "MCM3", "MCM7"),
    "TGFb/Integrin"        = c("TGFBR1", "TGFBR2", "ITGB1", "SMAD3", "SMAD4"))
  
subgroup_colors <- c(
    "Complex I (NDUF)"     = "#1B9E77",
    "Complex III/IV"       = "#D95F02",
    "ATP Synthase"         = "#7570B3",
    "Mito Quality/Markers" = "#E7298A",
    "NF-kB/Inflam"         = "#66A61E",
    "Chemokines"           = "#E6AB02",
    "Stress/MAPK"          = "#A6761D",
    "mTOR/Glycolysis"      = "#666666",
    "Cell Cycle"           = "#4575B4",
    "TGFb/Integrin"        = "#D73027")
  
up_pathways <- c("Complex I (NDUF)", "Complex III/IV",
                   "ATP Synthase",     "Mito Quality/Markers")
  
# ==============================================================================
# METADATA
# ==============================================================================
meta_paired <- combined_sample_info %>%
    mutate(condition = factor(ifelse(str_starts(sample, "Pre"), "Pre", "Post"),
                         levels = c("Pre", "Post")),
      group     = factor(case_when(
        str_detect(sample, "(?i)C_") ~ "control",
        str_detect(sample, "(?i)S")  ~ "sky"),
        levels = c("control", "sky")),
      SubjectID = str_extract(sample, "(?<=Pre|Post)[0-9]+")
    ) %>%
    filter(!is.na(SubjectID)) %>%
    group_by(SubjectID) %>% filter(n() == 2) %>% ungroup() %>%
    mutate(group_time = factor(paste(group, condition, sep = "_"),
                               levels = c("control_Pre", "control_Post",
                                          "sky_Pre",     "sky_Post"))) %>%
    arrange(group_time, SubjectID)
  
# ==============================================================================
  subgroup_df <- enframe(gene_groups, name = "Subgroup", value = "Symbol") %>%
    unnest(Symbol)
  valid_genes <- intersect(subgroup_df$Symbol, rownames(tpm_raw))
  
  plot_data <- master_results %>%
    filter(Symbol %in% valid_genes) %>%
    left_join(subgroup_df, by = "Symbol") %>%
    mutate(
      Subgroup  = factor(Subgroup, levels = names(gene_groups)),
      Direction = factor(ifelse(Subgroup %in% up_pathways,
                                "Upregulated", "Downregulated"),
                         levels = c("Upregulated", "Downregulated"))
    ) %>%
    arrange(Direction, Subgroup)
  
  expr_log2   <- log2(tpm_raw[plot_data$Symbol, meta_paired$sample] + 1)
  data_scaled <- t(scale(t(expr_log2)))
  colnames(data_scaled) <- meta_paired$SubjectID
  
  mat_up   <- data_scaled[plot_data$Direction == "Upregulated",   ]
  mat_down <- data_scaled[plot_data$Direction == "Downregulated", ]
  
col_z    <- colorRamp2(c(-2, 0, 2), c("#3A86FF", "white", "#FF006E"))
  
  # Mean Z-score barplot colors — matched to heatmap scale
Z_POS <- "#FF006E"   # positive z-score
Z_NEG <- "#3A86FF"   # negative z-score
  
  # LFC barplot colors — your new colors
UP_COL   <-  "#D62828"    # upregulated LFC
DOWN_COL <-   "#0077B6" # downregulated LFC
  
# ==============================================================================
ha_up_top <- HeatmapAnnotation(
    "Mean Z-score\n(Up genes)" = anno_barplot(
      means_up, baseline = 0,
      gp     = gpar(fill = ifelse(means_up > 0, Z_POS, Z_NEG), border = NA),
      height = unit(1.2, "cm")),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 7, fontface = "bold")
  )
  
ha_down_top <- HeatmapAnnotation(
    "Mean Z-score\n(Down genes)" = anno_barplot(
      means_down, baseline = 0,
      gp     = gpar(fill = ifelse(means_down > 0, Z_POS, Z_NEG), border = NA),
      height = unit(1.2, "cm")),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 7, fontface = "bold")
  )
  
# ==============================================================================
# HEATMAPS — LFC bars use UP_COL / DOWN_COL
# ==============================================================================
ht_up <- Heatmap(mat_up,
                   name                = "Mean Z-score\n(Up)",
                   show_heatmap_legend = FALSE,
                   top_annotation      = ha_up_top,
                   left_annotation     = ha_left_up,
                   right_annotation    = rowAnnotation(
                     "Up LFC" = anno_barplot(
                       plot_data[plot_data$Direction == "Upregulated", ]$Intervention_logFC,
                       baseline = 0,
                       gp       = gpar(fill = UP_COL, border = NA))),   # blue
                   col               = col_z,
                   column_split      = meta_paired$group_time,
                   cluster_columns   = FALSE,
                   cluster_rows      = FALSE,
                   row_split         = plot_data[plot_data$Direction == "Upregulated", ]$Subgroup,
                   row_title         = NULL,
                   show_column_names = FALSE,
                   row_names_side    = "left",
                   row_names_gp      = gpar(fontsize = 7, fontface = "italic"))
  
ht_down <- Heatmap(mat_down,
                     name                = "Mean Z-score\n(Down)",
                     show_heatmap_legend = TRUE,
                     top_annotation      = ha_down_top,
                     left_annotation     = ha_left_down,
                     right_annotation    = rowAnnotation(
                       "Down LFC" = anno_barplot(
                         plot_data[plot_data$Direction == "Downregulated", ]$Intervention_logFC,
                         baseline = 0,
                         gp       = gpar(fill = DOWN_COL, border = NA))),  # red
                     col                   = col_z,
                     column_split          = meta_paired$group_time,
                     cluster_columns       = FALSE,
                     cluster_rows          = FALSE,
                     row_split             = plot_data[plot_data$Direction == "Downregulated", ]$Subgroup,
                     row_title             = NULL,
                     show_column_names     = TRUE,
                     column_names_rot      = 0,
                     column_names_gp       = gpar(fontsize = 8, fontface = "plain"),
                     column_names_centered = TRUE,
                     row_names_side        = "left",
                     row_names_gp          = gpar(fontsize = 7, fontface = "italic"))
# ==============================================================================
# EXPORT
# ==============================================================================
final_path <- "RNAseq_res/RCT_Limma/SKY_Final_ZScore_Integrated_new.pdf"
  
pdf(final_path, width = 12, height = 10, useDingbats = FALSE)
  draw(ht_up %v% ht_down,
       column_title    = "SKY Treatment Effect (ITT - Pathway Genes)",
       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
       merge_legends   = TRUE,
       ht_gap          = unit(10, "mm"))
dev.off()


  
  
  
  
  
  
  
  
  
