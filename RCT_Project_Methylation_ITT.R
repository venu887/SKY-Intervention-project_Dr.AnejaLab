#@@@@@@@@@@@@@@@@@@@@@@@@@@@
# DNA methylation analysis using Limma
# reference gene annotation files downloaded from 
# Dataset posted on 2023-05-15 https://figshare.com/articles/dataset/Official_annotation_files_for_the_Illumina_Methylation_BeadChips_EPIC_v1_0_and_v2_0_/22819970
# Published February 26, 2025 Re-annotated EPICV2 are found in here:https://zenodo.org/records/14933469

# Illumina :https://support.illumina.com/downloads/infinium-methylationepic-v2-0-product-files.html 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10919401/
# https://www.cd-genomics.com/resource-dmc-dmr-dmg-in-dna-methylation-analysis.html
# https://www.cd-genomics.com/resource-comprehensive-guide-dna-methylation-arrays-epigenetic-research.html

# DNA Methylation Analysis Pipeline — SKY vs. Control RCT
# Overview
# This pipeline performs genome-wide differential DNA methylation analysis using Illumina EPIC v2 array data to identify CpG sites and pathways altered by SKY (Sudarshan Kriya Yoga) intervention compared to controls, following an Intent-to-Treat (ITT) interaction design.
# 
# Pipeline Steps
# 1. Data Loading & Sample Harmonization
# 
# Loads EPIC v2 beta values and sample sheet from .rds file
# Removes incomplete/advanced samples and standardizes sample naming conventions
# Retains only subjects with complete Pre/Post paired samples (n = 30 pairs: 12 Control, 18 SKY)
# 
# 2. Probe Quality Control
# Applies sequential filtering steps:
#   
#   Sex/mitochondrial chromosomes: Removes probes on chrX, chrY, chrM
# Missing data: Excludes probes with >20% missing values
# Zero variance: Removes non-variable probes
# SNP confounding: Removes probes with known SNPs within 10bp at MAF > 0.01 (Strategy B, publication-grade stringency)
# Remaining NAs are imputed to zero
# 
# 3. M-value Conversion
# 
# Beta values are clamped to [0.0001, 0.9999] to avoid infinite M-values
# Converted to M-values using sesame::BetaValueToMValue() for downstream linear modeling
# 
# 4. Limma Linear Mixed Model (ITT Design)
# 
# Models within-subject correlation using duplicateCorrelation() with subject ID as blocking factor
# Fits a ~ 0 + Group_Time design matrix with four levels: Control_Pre, Control_Post, SKY_Pre, SKY_Post
# Three contrasts tested:
#   
#   SKY Post vs. Pre
# Control Post vs. Pre
# ITT Interaction: (SKY_Post − SKY_Pre) − (Control_Post − Control_Pre)
# 
# 
# Results annotated with re-annotated EPIC v2 manifest (GENCODEv47)
# 
# 5. Visualization
# 
# Volcano plot: ITT interaction effect with directional coloring (hypermethylated/hypomethylated), top 10 probe labels, and ΔΔβ formula annotation
# Heatmap: Top 100 ITT CpGs across all samples, ordered by group/timepoint with row-scaled M-values
# 
# 6. Pathway Integration
# 
# Significant CpGs (p < 0.05, |ΔM| > 0.2) are mapped to genes via GENCODEv47 annotations
# Genes are linked to stress/inflammation pathways derived from a paired RNA-seq GSEA analysis
# Output table connects each significant CpG to its pathway, gene, genomic context, and methylation direction


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ LIMMA DE cpg, modified ITT (mITT)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls())
library(sesame)
library(minfi)
library(limma)
library(missMethyl)
library(DMRcate)
library(Gviz)
library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(GenomicRanges)
library(vroom)

file = "SKY_EpicV2_bvals_and_sample_sheet_new.rds"
data<-readRDS(file)
methyl<-data$bvals
sample_sheet<-data$samples
#Methylation data, remove Advanced samples and make uniform names
xx <- c("Post_Control_46", "Pre_Control_46")
sample_sheet <- sample_sheet[!sample_sheet$Sample_Name %in% xx, ]
methyl <- methyl[, !colnames(methyl) %in% xx]
colnames(methyl) <- gsub("SKY", "Sky", colnames(methyl))
colnames(methyl) <- gsub("_0([0-9])", "_\\1", colnames(methyl))

sample_sheet$Sample_Name <- gsub("SKY", "Sky", sample_sheet$Sample_Name)
sample_sheet$Sample_Name <- gsub("_0([0-9])", "_\\1", sample_sheet$Sample_Name)

sample_sheet$Treatment <- ifelse(
  grepl("^Pre", sample_sheet$Sample_Name, ignore.case = TRUE),
  "Pre", ifelse(grepl("^Post", sample_sheet$Sample_Name, ignore.case = TRUE), "Post",NA))
table(sample_sheet$GROUP2, sample_sheet$Treatment)
sample_sheet$Group <- factor(paste0(sample_sheet$Treatment,"_",sample_sheet$GROUP2),
                             levels = c("Pre_Control", "Pre_SKY", "Post_Control", "Post_SKY"))  

# 1. Identify SUBJECT3 IDs that have a complete pair (count of 2)
subject_counts <- table(sample_sheet$SUBJECT3)
complete_subjects <- names(subject_counts[subject_counts == 2])
# 2. Subset the sample_sheet to keep only these 37 pairs
sample_sheet <- sample_sheet[sample_sheet$SUBJECT3 %in% complete_subjects, ]

common_samples <- intersect(colnames(methyl), sample_sheet$Sample_Name)
sample_sheet <- sample_sheet[sample_sheet$Sample_Name %in% common_samples, ]
subject_counts <- table(sample_sheet$SUBJECT3)
complete_subjects <- names(subject_counts[subject_counts == 2])
sample_sheet <- sample_sheet[sample_sheet$SUBJECT3 %in% complete_subjects, ]

methyl <- methyl[, sample_sheet$Sample_Name]
print(paste("Final paired samples for analysis:", ncol(methyl)))
table(sample_sheet$CONDITION, sample_sheet$Sample_Group6)
# Group_Post Group_Pre
# Control         12        12
# SKY             18        18

methyl_CpG <- readRDS("/Users/mekalav/Documents/UAB/SKY/Data/Epigenetic/EPICv2_reannotated.rds")
table(methyl_CpG$CHR)

# ============================================================
# PROBE QC: Remove sex chromosomes, low quality, 
#           zero variance, AND SNP-confounded probes
# ============================================================
# Step 1: Remove sex chromosome and unassigned probes
methyl_CpG$cg_id <- sub("_.*$", "", methyl_CpG$IlmnID)
remove_cpgs <- methyl_CpG$cg_id[
  methyl_CpG$CHR %in% c("chrX", "chrY", "chrM", "chr0")]
methyl <- methyl[!rownames(methyl) %in% remove_cpgs, ]
print(paste("After sex/mito chr removal:", nrow(methyl)))
# Expected: 898,467

# Step 2: Remove probes with >20% missing values
keep_probes <- rowMeans(is.na(methyl)) <= 0.2
methyl <- methyl[keep_probes, ]
print(paste("After NA filtering:", nrow(methyl)))

# Step 3: Remove zero-variance probes
library(matrixStats)
keep_var <- rowVars(as.matrix(methyl), na.rm = TRUE) > 0
methyl <- methyl[keep_var, ]
print(paste("After variance filtering:", nrow(methyl)))

# Step 4: Remove SNP-confounded probes
# -------------------------------------------------------
# Strategy A: Remove probes with ANY known SNP at the 
# CpG site (distance = 0) regardless of MAF
# -------------------------------------------------------
snp_info <- methyl_CpG %>%
  mutate(cg_id = sub("_.*$", "", IlmnID)) %>%
  filter(cg_id %in% rownames(methyl)) %>%
  dplyr::select(cg_id, SNP_ID, SNP_DISTANCE, SNP_MinorAlleleFrequency)

# Identify probes with SNP AT the CpG site (distance == 0)
# SNP_DISTANCE is semicolon-separated - check if any distance is 0
snp_at_cpg <- snp_info %>%
  filter(!is.na(SNP_ID) & SNP_ID != "") %>%
  rowwise() %>%
  mutate(
    distances = list(as.numeric(
      unlist(strsplit(as.character(SNP_DISTANCE), ";")))),
    mafs = list(as.numeric(
      unlist(strsplit(as.character(SNP_MinorAlleleFrequency), ";")))),
    # Flag if SNP is at position 0 (CpG site) with MAF > 0.01
    has_snp_at_cpg = any(distances == 0 & mafs > 0.01, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(has_snp_at_cpg) %>%
  pull(cg_id)

print(paste("Probes with SNP at CpG site (MAF>0.01):", 
            length(snp_at_cpg)))

# -------------------------------------------------------
# Strategy B (more stringent): Remove probes with any SNP
# within 10bp AND MAF > 0.01
# -------------------------------------------------------
snp_near_cpg <- snp_info %>%
  filter(!is.na(SNP_ID) & SNP_ID != "") %>%
  rowwise() %>%
  mutate(
    distances = list(as.numeric(
      unlist(strsplit(as.character(SNP_DISTANCE), ";")))),
    mafs = list(as.numeric(
      unlist(strsplit(as.character(SNP_MinorAlleleFrequency), ";")))),
    # Flag if SNP within 10bp with MAF > 0.01
    has_nearby_snp = any(distances <= 10 & mafs > 0.01, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(has_nearby_snp) %>%
  pull(cg_id)

print(paste("Probes with SNP within 10bp (MAF>0.01):", 
            length(snp_near_cpg)))

# Apply SNP removal - use Strategy B for publication
methyl <- methyl[!rownames(methyl) %in% snp_near_cpg, ]
print(paste("After SNP probe removal:", nrow(methyl)))

# Step 5: Impute remaining NAs and remove
methyl[is.na(methyl)] <- 0
methyl <- na.omit(methyl)
print(paste("Final probe count:", nrow(methyl)))

#@@@@@@@@@@@@@@@@@ DM CpGs and Regions in Sesame
library(sesame)
dim(methyl)
class(methyl)
# 1. Clamp the beta values to avoid exactly 0 or 1
# pmax(x, 0.0001) ensures no value is below 0.0001
# pmin(x, 0.9999) ensures no value is above 0.9999
methyl_clamped <- pmax(pmin(as.matrix(methyl), 0.9999), 0.0001)

m_value <- BetaValueToMValue(methyl_clamped)
sum(is.infinite(m_value))

# Identify which samples exist in BOTH objects
sample_sheet<-as.data.frame(sample_sheet)
row.names(sample_sheet)=sample_sheet$Sample_Name

sample_sheet$Sample_Group=as.factor(as.vector(sample_sheet$Sample_Group))
sample_sheet$X=as.factor(as.vector(sample_sheet$CONDITION))
saveRDS(m_value,"ITT_limma/m_value_limma.rds") # Save processed data

sample_sheet$Group_Fixed <- ifelse(grepl("Sky", sample_sheet$Sample_Name, ignore.case = TRUE), "SKY", "Control")
sample_sheet$Group_Fixed <- factor(sample_sheet$Group_Fixed, levels = c("Control", "SKY"))

sample_sheet$Timepoint_Fixed <- ifelse(grepl("Pre", sample_sheet$Sample_Name, ignore.case = TRUE), "Pre", "Post")
sample_sheet$Timepoint_Fixed <- factor(sample_sheet$Timepoint_Fixed, levels = c("Pre", "Post"))
table(sample_sheet$Group_Fixed, sample_sheet$Timepoint_Fixed)

sample_sheet$Group_Time <- factor(paste(sample_sheet$Group_Fixed, sample_sheet$Timepoint_Fixed, sep="_"))
design <- model.matrix(~ 0 + Group_Time, data = sample_sheet)
colnames(design) <- levels(sample_sheet$Group_Time)
corfit <- duplicateCorrelation(m_value, design, block = sample_sheet$SUBJECT3)
fit <- lmFit(m_value, design, block = sample_sheet$SUBJECT3, correlation = corfit$consensus)


# 4. Define specific Contrasts
# Sky_Post vs Sky_Pre
# Control_Post vs Control_Pre
# ITT (The difference of the differences)
cont.matrix <- makeContrasts(
  Sky_Post_vs_Pre = SKY_Post - SKY_Pre,
  Control_Post_vs_Pre = Control_Post - Control_Pre,
  ITT_Interaction = (SKY_Post - SKY_Pre) - (Control_Post - Control_Pre),
  levels = design)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

sky_res <- topTable(fit2, coef = "Sky_Post_vs_Pre", number = Inf, sort.by = "none")
con_res <- topTable(fit2, coef = "Control_Post_vs_Pre", number = Inf, sort.by = "none")
itt_res <- topTable(fit2, coef = "ITT_Interaction", number = Inf, sort.by = "none")

colnames(sky_res) <- paste0("Sky_", colnames(sky_res))
colnames(con_res) <- paste0("Control_", colnames(con_res))
colnames(itt_res) <- paste0("ITT_", colnames(itt_res))

combined_wide <- cbind(sky_res, con_res, itt_res)

combined_wide$Probe_ID <- rownames(combined_wide)

annotated_final <- combined_wide %>%
  left_join(methyl_CpG, by = c("Probe_ID" = "cg_id")) %>%
  select(Probe_ID, everything()) 

write.csv(annotated_final, "ITT_limma/ITT_limma_Volcano.csv", row.names = T)




# ============================================================
# VOLCANO PLOT: Limma ITT Differential Methylation
# ============================================================
library(ggplot2)
library(ggrepel)
library(dplyr)
rm(list = ls())
itt_results<-read.csv("ITT_limma/ITT_limma_Volcano.csv", row.names = 1)
head(annotated_ITT)
itt_results <- annotated_ITT

volcano_data <- itt_results %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      P.Value < 0.05 & logFC >  0.2 ~ "Hypermethylated",
      P.Value < 0.05 & logFC < -0.2 ~ "Hypomethylated",
      TRUE                           ~ "Not Significant"))

# ---- Dynamic counts ----
n_hyper <- sum(volcano_data$Direction == "Hypermethylated")
n_hypo  <- sum(volcano_data$Direction == "Hypomethylated")
n_ns    <- sum(volcano_data$Direction == "Not Significant")

cat("Hypermethylated (p<0.05, logFC>0.2):", n_hyper, "\n")
cat("Hypomethylated  (p<0.05, logFC<-0.2):", n_hypo,  "\n")
cat("Not Significant:", n_ns, "\n")

# ---- Add direction labels ----
volcano_data <- volcano_data %>%
  dplyr::mutate(
    Direction_label = dplyr::case_when(
      Direction == "Hypermethylated" ~ paste0("Hypermethylated (n=", n_hyper, ")"),
      Direction == "Hypomethylated"  ~ paste0("Hypomethylated (n=",  n_hypo,  ")"),
      TRUE                           ~ paste0("Not Significant (n=",  n_ns,   ")") ))

# ---- Define colors using setNames (fixes the parse error) ----
color_values <- setNames(c("dodgerblue3", "firebrick3", "grey80"),
  c(paste0("Hypermethylated (n=", n_hyper, ")"),
    paste0("Hypomethylated (n=",  n_hypo,  ")"),
    paste0("Not Significant (n=",  n_ns,   ")")))

# ---- Select top 10 significant probes to label ----
top_probes <- volcano_data %>%
  dplyr::filter(Direction != "Not Significant") %>%
  dplyr::arrange(P.Value) %>%
  dplyr::slice_head(n = 10)

# ---- Build volcano plot ----
p_volcano <- ggplot(volcano_data,
                    aes(x     = logFC,
                        y     = -log10(P.Value),
                        color = Direction_label)) +
  
  # Non-significant points first (background layer)
  geom_point(
    data  = dplyr::filter(volcano_data, Direction == "Not Significant"),
    alpha = 0.3, size = 0.6) +
  
  # Significant points on top (foreground layer)
  geom_point(
    data  = dplyr::filter(volcano_data, Direction != "Not Significant"),
    alpha = 0.7, size = 1.2) +
  
  # Significance threshold line
  geom_hline(
    yintercept = -log10(0.05),
    linetype   = "dashed",
    color      = "grey30",
    linewidth  = 0.5) +
  
  # Effect size threshold lines
  geom_vline(
    xintercept = c(-0.2, 0.2),
    linetype   = "dashed",
    color      = "grey30",
    linewidth  = 0.5) +
  
  # Top probe labels
  geom_text_repel(
    data          = top_probes,
    aes(label     = Probe_ID),
    size          = 3,
    fontface      = "italic",
    box.padding   = 0.5,
    max.overlaps  = 20,
    segment.color = "grey50",
    segment.size  = 0.3,
    color         = "black") +
  
  # DiD formula annotation (top right)
  annotate(
    "text",
    x        = max(volcano_data$logFC) * 0.55,
    y        = max(-log10(volcano_data$P.Value)) * 0.95,
    label    = "Delta*Delta*Beta==(SKY[Post]-SKY[Pre])-(Ctrl[Post]-Ctrl[Pre])",
    parse    = TRUE,
    size     = 3.2,
    color    = "grey20",
    fontface = "italic") +
  
  # Total probe count annotation (top left)
  annotate(
    "text",
    x     = min(volcano_data$logFC) * 0.85,
    y     = max(-log10(volcano_data$P.Value)) * 0.95,
    label = paste0("Total probes tested: ",
                   format(nrow(volcano_data), big.mark = ",")),
    size  = 3.2,
    color = "grey30") +
  
  scale_color_manual(values = color_values) +
  
  theme_classic(base_size = 13) +
  labs(
    title    = "Differential DNA Methylation CpGs: ITT Interaction Effect",
    subtitle = paste0("SKY vs Control  |  Nominal p < 0.05  |  ",
                      "|ΔM| > 0.2 threshold"),
    x        = expression(log~Fold~Change~(M-value~scale)),
    y        = expression(-log[10](P-value)),
    color    = "Methylation Direction"
  ) +
  theme(
    legend.position  = "top",
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 10, face = "bold"),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(color = "grey40", size = 9),
    panel.grid.major = element_line(color = "grey95"),
    axis.text        = element_text(size = 10)
  )

print(p_volcano)

ggsave("ITT_limma/ITT_limma_Volcano.pdf", p_volcano, width = 10, height = 8, dpi = 300)



#==============================================================================
# HEATMAP of Top CpGs
#==============================================================================
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(grid)
m_value <- readRDS("ITT_limma/m_value_limma.rds")
top_probes <- itt_results %>% 
  arrange(ITT_P.Value) %>% 
  head(100) %>% 
  pull(Probe_ID)

plot_data <- m_value[top_probes, ]

# 2. Create the "CpG ID (Gene Name)" labels
row_labels_df <- itt_results %>%
  filter(Probe_ID %in% rownames(plot_data)) %>%
  mutate(Short_Gene = sapply(strsplit(as.character(GENCODEv47_Gene_Name), ";"), `[`, 1),
    # Format: cg000... (Gene)
    Display_Label = ifelse(is.na(Short_Gene) | Short_Gene == "" | Short_Gene == " ", 
                           as.character(Probe_ID), 
                           paste0(Probe_ID, " (", Short_Gene, ")")))

# Named vector for mapping
label_vector <- row_labels_df$Display_Label
names(label_vector) <- row_labels_df$Probe_ID

# Build and Sort Annotations
anno_col <- data.frame(
  row.names = colnames(plot_data),
  Group = ifelse(grepl("Sky", colnames(plot_data), ignore.case = TRUE), "SKY", "Control"),
  Timepoint = ifelse(grepl("Pre", colnames(plot_data), ignore.case = TRUE), "Pre", "Post")
) %>%
  mutate(group_time = factor(paste(Group, Timepoint, sep = "_"), 
                             levels = c("Control_Pre", "Control_Post", "SKY_Pre", "SKY_Post")))

# Define specific colors
pal_anno <- c(
  Control_Pre  = "#E8D5C4",
  Control_Post = "#9C6644",
  SKY_Pre      = "#CAF0F8",
  SKY_Post     = "#0077B6")

# Order columns so groups sit together
sample_order <- rownames(anno_col[order(anno_col$group_time), ])
plot_data_sorted <- plot_data[, sample_order]
anno_col_sorted <- anno_col[sample_order, "group_time", drop = FALSE]

pdf("ITT_limma/ITT_limma.pdf", width = 11, height = 9)
pheatmap(plot_data_sorted,
         main = "Top 100 Probes: SKY vs Control ITT",
         cluster_rows = TRUE,         
         cluster_cols = FALSE,        
         show_colnames = FALSE,       
         annotation_col = anno_col_sorted,   
         annotation_colors = list(group_time = pal_anno),
         labels_row = label_vector[rownames(plot_data_sorted)], 
         scale = "row",               
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = NA,
         fontsize_row = 6,
         gaps_col = which(diff(as.numeric(anno_col_sorted$group_time)) != 0))
dev.off()




#==============================================================================
# Significant CpGs annotation with associated pathway from RNAseq
#==============================================================================
rm(list = ls())
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggrepel)
itt_results<- read.csv("ITT_limma/ITT_limma_Volcano.csv", row.names=1)
head(itt_results)
pathways_filtered <- read.csv("RCT_Limma/SKY_Stress_Inflammation_Pathways.csv",row.names = 1)
head(pathways_filtered)
#==============================================================================
# 2. PATHWAY-GENE-CpG INTEGRATION
#    KEY FIX: Clean gene names BEFORE joining to fix missing CpG attachment
#==============================================================================
# Build pathway-gene map
pathway_gene_map <- pathways_filtered %>%
  #  dplyr::filter(Pathway_ID %in% target_pathways) %>%
  dplyr::select(Pathway_ID, Clean_Name, leadingEdge) %>%
  dplyr::mutate(
    Clean_Name = gsub("HALLMARK: ", "", Clean_Name),
    leadingEdge = gsub(";\\s*", "; ", leadingEdge)
  ) %>%
  tidyr::separate_rows(leadingEdge, sep = ";\\s*") %>%
  dplyr::mutate(leadingEdge = trimws(leadingEdge)) %>%
  dplyr::rename(Gene_Symbol = leadingEdge) %>%
  dplyr::filter(Gene_Symbol != "" & !is.na(Gene_Symbol))

cat("Pathway gene map rows:", nrow(pathway_gene_map), "\n")
cat("Unique genes:", n_distinct(pathway_gene_map$Gene_Symbol), "\n")

itt_expanded <- itt_results %>%
  # Select only needed columns to reduce memory
  dplyr::select(Probe_ID, logFC, P.Value, adj.P.Val, 
                CHR, MAPINFO,
                GENCODEv47_Gene_Name,
                Relation_to_UCSC_CpG_Island,
                GENCODEv47_Feature_Type) %>%
  # Standardize semicolons
  dplyr::mutate(
    GENCODEv47_Gene_Name = gsub(";\\s*", "; ", GENCODEv47_Gene_Name)
  ) %>%
  # Expand one row per gene
  tidyr::separate_rows(GENCODEv47_Gene_Name, sep = ";\\s*") %>%
  dplyr::mutate(
    GENCODEv47_Gene_Name = trimws(GENCODEv47_Gene_Name)
  ) %>%
  dplyr::filter(
    GENCODEv47_Gene_Name != "" & !is.na(GENCODEv47_Gene_Name)
  )

cat("Expanded ITT rows:", nrow(itt_expanded), "\n")

# Join CpGs to pathway genes
itt_pathway_cpgs <- itt_expanded %>%
  dplyr::inner_join(
    pathway_gene_map,
    by = c("GENCODEv47_Gene_Name" = "Gene_Symbol")
  ) %>%
  dplyr::distinct(Probe_ID, Pathway_ID, .keep_all = TRUE)

cat("CpGs mapped to pathways:", nrow(itt_pathway_cpgs), "\n")
cat("Per pathway:\n")
print(table(itt_pathway_cpgs$Clean_Name))

# Filter for significant CpGs only
itt_sig <- itt_pathway_cpgs %>%
  dplyr::filter(P.Value < 0.05) %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      logFC >  0.2 ~ "Hypermethylated",
      logFC < -0.2 ~ "Hypomethylated",
      TRUE         ~ "Minimal Change"
    )
  )

cat("\nSignificant CpGs per pathway:\n")
print(table(itt_sig$Clean_Name))

write.csv(itt_sig, "ITT_limma/Sig_CpGs_43_pathways.csv", row.names = T)
