#@@@@@@@@@@@@@@@@@@@@@@@@
# SKY RCT study
#@@@@@@@@@@@@@@@@@@@@@@@@
#++++++++++++++++++++++++++
# First panel Survey scores of participant groups
#++++++++++++++++++++++++++
library(tidyverse)
library(readxl)
library(ggplot2)
library(tidyverse)
library(readxl)
library(patchwork)
library(ggpubr)

df_raw <- read_excel("Between_Group_Estimates.xlsx")

all_data_final <- df_raw %>%
  filter(!is.na(.[[2]])) %>%
  mutate(
    Outcome_Raw = str_squish(as.character(.[[1]])),
    MD_raw      = as.character(.[[2]]),
    CI95_raw    = as.character(.[[3]]),
    P_raw       = as.character(.[[6]]),
    AdjP_raw    = as.character(.[[7]]),
    # Numeric conversions
    MD_num   = as.numeric(str_extract(MD_raw, "-?[0-9.]+")),
    AdjP_num = as.numeric(str_extract(AdjP_raw, "[0-9.]+")),
    RawP_num = as.numeric(str_extract(P_raw, "[0-9.]+")),
    CI_low   = as.numeric(str_extract(CI95_raw, "(?<=\\[)-?[0-9.]+")),
    CI_high  = as.numeric(str_extract(CI95_raw, "-?[0-9.]+(?=\\])")),
       Category = case_when(
      Outcome_Raw %in% c("Burnout", "Perceived stress", "Overall score", "General distress", 
                         "Anxious arousal", "Anhedonicdepression", "Negative affect") ~ "Mental Health",
      
      Outcome_Raw %in% c("Positive affect", "Eudaimonic well-being", "Life satisfaction", 
                         "Gratitude", "Optimism", "Self-esteem", "Self-compassion", 
                         "Social connectedness", "Mindfulness") ~ "Psychological Thriving",
      
      TRUE ~ "Personality, Coping & Emotion"
    ),
        Outcome_Display = case_when(
      Outcome_Raw == "Overall score" & Category == "Mental Health" ~ "MASQ D30 – Overall anxiety and depression",
      Outcome_Raw == "General distress" ~ "MASQ D30 Subscale – General distress",
      Outcome_Raw == "Anxious arousal" ~ "MASQ D30 Subscale – Anxious arousal",
      Outcome_Raw == "Anhedonicdepression" ~ "MASQ D30 Subscale – Anhedonic depression",
      Outcome_Raw == "Positive affect" ~ "PANAS - Positive affectivity",
      Outcome_Raw == "Negative affect" ~ "PANAS - Negative affectivity",
      Outcome_Raw == "Overall score" & Category == "Psychological Thriving" ~ "RYFF - Overall psychological well-being",
      Outcome_Raw == "Self-compassion" ~ "SCS-SF – Overall Self-Compassion",
      Outcome_Raw == "Mindfulness" ~ "FFMQ-15 – Overall Mindfulness",
      Outcome_Raw == "Avoidant coping" ~ "Coping styles - Avoidant coping",
      Outcome_Raw == "Problem-focused" ~ "Coping styles - Problem-focused",
      Outcome_Raw == "Emotion-focused" ~ "Coping styles - Emotion-Focused",
      Outcome_Raw == "Extraversion" ~ "Personality Factors - Extraversion",
      Outcome_Raw == "Agreeableness" ~ "Personality Factors - Agreeableness",
      Outcome_Raw == "Conscientiousness" ~ "Personality Factors - Conscientiousness",
      Outcome_Raw == "Neuroticism" ~ "Personality Factors - Neuroticism",
      Outcome_Raw == "Openness" ~ "Personality Factors - Openness",
      Outcome_Raw == "Emotional intelligence" ~ "Emotion Mindset - Emotional Intelligence",
      Outcome_Raw == "Beliefs/opinions about emotion" ~ "Emotion Mindset - Belief/Opinions about Emotion",
      TRUE ~ Outcome_Raw
    ),
    
    SigColor = ifelse(!is.na(RawP_num) & RawP_num < 0.05, "p < 0.05", "p >= 0.05"),
    PlotOrder = row_number()
  )

# =============================================================================
# FIGURE PANEL HRV: MENTAL HEALTH
# =============================================================================
df_A <- all_data_final %>% 
  filter(Category == "Mental Health")
FigA <- ggplot(df_4A, aes(x = MD_num, 
                           y = reorder(Outcome_Display, -PlotOrder), 
                           color = SigColor)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  # Directional Annotations
  annotate("label", x = -10, y = Inf, label = "← FAVORS SKY", fill = "#0077B6", 
           color = "white", fontface = "bold", vjust = 0.5, size = 4) +
  annotate("label", x = 10, y = Inf, label = "FAVORS CONTROL →", fill = "#9C6644", 
           color = "white", fontface = "bold", vjust = 0.5, size = 4) +
  annotate("text", x = 16, y = Inf, label = "Mean Diff [95% CI]", 
           color = "black", fontface = "bold", vjust = 0.5, size = 4) +
  # Horizontal Error Bars (Modern ggplot2 syntax)
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), height = 0.2, linewidth = 0.7) + 
  geom_point(size = 2.5) +
  # Text labels for the exact values on the right
  geom_text(aes(x = 16, label = sprintf("%.2f [%.2f, %.2f]", MD_num, CI_low, CI_high)), 
            color = "black", size = 3.5) +
  # Manual color mapping for significance
  scale_color_manual(values = c("p < 0.05" = "orangered", "p >= 0.05" = "black")) +
  theme_classic() + 
  coord_cartesian(xlim = c(-18, 22), clip = "off") + # 'clip = off' allows text in margins
  
  labs(title = "Effect of SKY on Mental Health", 
       subtitle = "Mean Differences (SKY - Control) | Linear regression model adjusted for baseline",
       x = "Mean Difference (95% CI)", y = "", color = "Significance") +
  
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, margin = margin(t = 10), color = "black"),
        plot.title = element_text(face = "bold", size = 14, color = "black"),
        legend.position = "top", 
        plot.margin = margin(40, 100, 20, 20)) # Extra space for the side labels

print(FigA)

ggsave("FigA_Mental_Health.pdf", Fig4A, width = 8, height = 4)


# =============================================================================
# FIGURE B: PSYCHOLOGICAL THRIVING
# =============================================================================
df_B <- all_data_final %>% 
  filter(Category == "Psychological Thriving")
FigB <- ggplot(df_4B, aes(x = MD_num, 
                           y = reorder(Outcome_Display, -PlotOrder), # Fixed: Use - instead of desc()
                           color = SigColor)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  # Directional Headers
  annotate("label", x = -2, y = Inf, label = "← FAVOR CONTROL", fill = "#9C6644", 
           color = "white", fontface = "bold", vjust = 0, size = 4) +
  annotate("label", x = 7, y = Inf, label = "FAVOR SKY →", fill = "#0077B6", 
           color = "white", fontface = "bold", vjust = 0, size = 4) +
  # Table Header
  annotate("text", x = 22, y = Inf, label = "Mean Diff [95% CI]", 
           color = "black", fontface = "bold", vjust = 0.5, size = 4) +
  # Fixed: Use geom_errorbar (no 'h')
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), height = 0.2, linewidth = 0.7) + 
  geom_point(size = 2.5) +
  # Table Content
  geom_text(aes(x = 22, label = sprintf("%.2f [%.2f, %.2f]", MD_num, CI_low, CI_high)), 
            color = "black", size = 3.5) +
  scale_color_manual(values = c("p < 0.05" = "orangered", "p >= 0.05" = "black")) +
  theme_classic() + 
  coord_cartesian(xlim = c(-5, 28), clip = "off") + # 'clip = off' is vital for your side-table
  
  labs(title = "Effect of SKY on Psychological Thriving", 
       subtitle = "Mean Differences (SKY - Control) | Linear regression model adjusted for baseline",
       x = "Mean Difference (95% CI)", y = "", color = "Significance") +
  
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, margin = margin(t = 10)),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "top", 
        plot.margin = margin(40, 100, 20, 20))


print(FigB)

ggsave("FigB_Thriving.pdf", FigB, width = 8, height = 4)

# =============================================================================
# FIGURE C: PERSONALITY, COPING & EMOTION
# =============================================================================
df_C <- all_data_final %>% 
  filter(Category == "Personality, Coping & Emotion") %>%
  filter(!Outcome_Raw %in% c("Negative affect", "Positive affect"))

FigC <- ggplot(df_4C, aes(x = MD_num, 
                           y = reorder(Outcome_Display, -PlotOrder), # Fixed: Use - instead of desc()
                           color = SigColor)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  
  # Table Header
  annotate("text", x = 15, y = Inf, label = "Mean Diff [95% CI]", 
           color = "black", fontface = "bold", vjust = 0.5, size = 4) +
  
  # Fixed: Use geom_errorbar (no 'h')
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), height = 0.2, linewidth = 0.7) + 
  geom_point(size = 2.5) +
  
  # Table Content
  geom_text(aes(x = 15, label = sprintf("%.2f [%.2f, %.2f]", MD_num, CI_low, CI_high)), 
            color = "black", size = 3.5) +
  
  scale_color_manual(values = c("p < 0.05" = "orangered", "p >= 0.05" = "black")) +
  theme_classic() + 
  coord_cartesian(xlim = c(-5, 21), clip = "off") + # Vital for side-table display
  
  labs(title = "Effect of SKY on Personality, Coping & Emotion", 
       subtitle = "Mean Differences (SKY - Control) | Linear regression model adjusted for baseline",
       x = "Mean Difference (95% CI)", y = "", color = "Significance") +
  
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, margin = margin(t = 10)),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "top", 
        plot.margin = margin(40, 100, 20, 20))

print(FigC)

ggsave("FigC_Personality.pdf", Fig4C, width = 10, height = 4)




# =============================================================================
# FIGURE D: Box plots of Significant Survey scores PERSONALITY, COPING & EMOTION
# =============================================================================
file3 <- "SKY_BehavioralData_Adele_07.03.24_SurveyScores.xlsx"
data3 <- read_excel(file3, sheet = 3)

data3 <- data3 %>% mutate(ID = row_number())
colnames(data3)

vars_sig <- list(
  "PerceivedStressScore_Total"       = "Perceived Stress",
  "MASQScore_Composite"              = "MASQ - Overall Anxiety/Depression",
  "MASQScore_AnxiousArousal"         = "MASQ - Anxious Arousal",
  "PANAS_PositiveAffectScore"        = "PANAS - Positive Affect",
  "RYFFScore_SumTotal"               = "Psychological Well-being (RYFF)",
  "LOTRScore_Total"                  = "Optimism (LOT-R)",
  "SelfCompassionScoreSum_Total"     = "Self-Compassion (SCS-SF)",
  "FFMQScoreSum_Total"               = "Mindfulness (FFMQ-15)",
  "EudWellBeingScore_Total"          = "Eudaimonic Well-being",
  "BFI44ScoreMean_Extraversion"      = "Extraversion"
)

traj_plots <- list()

# 2. Loop to Generate Plots
for(v in names(vars_sig)) {
  pre_col  <- paste0("Pre_", v)
  post_col <- paste0("Post_", v)
  
  if(!all(c(pre_col, post_col) %in% colnames(data3))) next
  
  # Prepare data (Explicitly using dplyr::select to avoid MASS conflict)
  df_long <- data3 %>%
    dplyr::select(ID, Treatment, all_of(c(pre_col, post_col))) %>%
    drop_na() %>%
    pivot_longer(cols = starts_with(c("Pre", "Post")), 
                 names_to = "Time", 
                 values_to = "Score") %>%
    mutate(
      Time = factor(gsub(paste0("_", v), "", Time), levels = c("Pre", "Post")),
      Treatment = factor(Treatment, levels = c("Control", "SKY"))
    )
  
  df_medians <- df_long %>%
    group_by(Treatment, Time) %>%
    summarise(MedianScore = median(Score, na.rm = TRUE), .groups = "drop")
  
  # Plotting with Wilcoxon aesthetics
  traj_plots[[v]] <- ggplot(df_long, aes(x = Time, y = Score, color = Treatment)) +
    geom_line(aes(group = ID), alpha = 0.5, linewidth = 0.5) + 
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.3, color = "black", fill = "gray95") +
    stat_compare_means(method = "wilcox.test", paired = TRUE, label = "p.format", 
                       label.x = 1.3, size = 3, color = "black") +
    geom_line(data = df_medians, aes(x = Time, y = MedianScore, group = Treatment), 
              color = "firebrick3", linewidth = 0.8) + 
    geom_point(data = df_medians, aes(x = Time, y = MedianScore), color = "firebrick3", size = 1.5) +
    facet_wrap(~Treatment) +
    scale_color_manual(values = c("Control" = "#9C6644", "SKY" = "#0077B6")) +
    labs(title = vars_sig[[v]], x = "", y = "Score") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
}

# 3. Assemble and Save Grid
final_grid <- wrap_plots(traj_plots, ncol = 2) + 
  plot_annotation(title = "Significant Outcomes: Pre-Post Individual Changes",
                  subtitle = "P-values: Paired Wilcoxon Signed-Rank Test | Red Line: Group Median")

ggsave("Significnat_Trajectories_with_Boxes_WT.pdf", 
       final_grid, width = 4, height = 7)




# =============================================================================
# FIGURE E: Treatment Effect by inflamation status
# =============================================================================
Psychosocial<-read.csv("sky_psychosocial_March2026.csv")
Psychosocial<-Psychosocial[Psychosocial$treatment %in% "SKY-Control",]
Psychosocial<-Psychosocial[Psychosocial$time %in% "post-pre",]
Psychosocial<-Psychosocial[!Psychosocial$inflam_grp %in% "Elevated-Low",]
str(Psychosocial)
library(tidyverse)

plot_data <- Psychosocial %>%
  mutate(p_numeric = as.numeric(p),
    LabelText = paste0("p=", sprintf("%.3f", p_numeric)),
    Measure_Clean = case_when(
      measure == "PerceivedStressScore_Total"    ~ "Perceived Stress Score",
      measure == "MASQScore_Composite"           ~ "MASQ Overall Score",
      measure == "MASQScore_GenDistress"         ~ "General Distress",
      measure == "MASQScore_AnxiousArousal"      ~ "Anxious Arousal",
      measure == "MASQScore_AnhedonicDepression" ~ "Anhedonic Depression",
      measure == "SCSRScoreSum_Total"            ~ "Social Connectedness",
      TRUE ~ as.character(measure)))

final_plot <- ggplot(plot_data, aes(x = estimate, y = Measure_Clean, 
                                    color = p_numeric < 0.05, 
                                    linetype = inflam_grp,
                                    group = interaction(Measure_Clean, inflam_grp))) +
  # Zero reference line
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.8, color = "black") +
  # Error bars: inherits color and custom linetype
  geom_errorbarh(aes(xmin = lower, xmax = upper), 
                 height = 0.5, 
                 linewidth = 0.5,
                 position = position_dodge(width = 0.7)) +
  # Point estimates
  geom_point(size = 2.5, position = position_dodge(width = 0.7)) +
  # COLOR: Red for significant, Black for others
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "p >= 0.05")) +
  # "solid" is a standard straight line
  scale_linetype_manual(values = c("Elevated" = "13", "Low" = "solid")) +
  
  theme_classic() +
  labs(
    title = "Treatment effects by inflammation status",
    subtitle = "Significant (p < 0.05) | Wide Dotted: Elevated, Straight: Low",
    x = "Mean Difference (95% CI)",
    y = "",
    color = "Significance",
    linetype = "Inflammation Status") +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 10,color = "black"),
        axis.title.x = element_text(size = 10), 
        plot.margin = margin(10, 10, 10, 10, "pt") # Reduced right margin since text is gone
  ) + coord_cartesian(clip = 'off')

print(final_plot)

pdf("Psychosocial_Inflam_status.pdf", height=3, width=4.5)
print(final_plot)
dev.off()
