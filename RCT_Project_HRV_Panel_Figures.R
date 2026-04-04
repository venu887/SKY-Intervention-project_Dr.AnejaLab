
#@@@@@@@@@@@@@@ 
# FIgure Panel Session HRV, SLeep HRV, inflamation HRV
# ==============================================================================
# SKY PROJECT: HRV & SLEEP DATA 
# ==============================================================================
rm(list = ls())
library(tidyverse)
mean_hrv<-read.csv("~/Documents/UAB/SKY/Data/HRV_Survay/FW_HRV_files/treatment_diff_75_data.csv")
sleep_elevate<-read.csv("~/Documents/UAB/SKY/Data/HRV_Survay/FW_HRV_files/sleep_hrv.csv")
hrv           <- readRDS("~/Documents/UAB/SKY/Data/HRV_Survay/FW_HRV_files/hrv.RDS")
overallsleep  <- readRDS("~/Documents/UAB/SKY/Data/HRV_Survay/FW_HRV_files/overallsleep.RDS")

clean_study_data <- function(df) {
  df %>%
    mutate(
      # Standardize Time Labels
      time = str_replace(time, "WS: ", "Day "),   # Handles "WS: 1" -> "Day 1" # Work Session (WS)
      time = str_replace(time, "^WS$", "Day 1"),  # Handles "WS" -> "Day 1" # Work Session (WS)
      time = str_replace(time, "FU: ", "Week "),  # Handles "FU: 1" -> "Week 1" # Follow up (FU)
      
      time = factor(time, levels = c(
        "Baseline", 
        paste("Day", 1:3), 
        paste("Week", 1:8)
      )),
      treatment = factor(treatment, levels = c("Control", "Sky", "Sky-Control"))
    ) %>%

    filter(!is.na(Estimate)) %>%
    arrange(treatment, time)
}

hrv_clean    <- clean_study_data(hrv)
sleep_clean  <- clean_study_data(overallsleep)


# ==============================================================================
# Figure A: Session HRV Mean Difference (SKY − Control) MINUTE-WISE HRV
# ==============================================================================
head(mean_hrv)
table(mean_hrv$grp)

colnames(mean_hrv)[4] <- "p_val"
sig_times <- mean_hrv %>%
  filter(grp == "SKY-Control" & p_val < 0.05) %>%
  pull(time)

FigA <- mean_hrv %>%
  filter(grp == "SKY-Control") %>%
  mutate(is_sig = ifelse(p_val < 0.05, "Significant", "Not Significant")) %>%
  ggplot(aes(x = time, y = Estimate)) +
  # A. Vertical Significance Lines
  geom_vline(data = . %>% filter(p_val < 0.05),
             aes(xintercept = time),
             linetype = "dotted", color = "gray80", linewidth = 0.5) +
  # B. Horizontal Reference Line at Zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  # C. Difference Ribbon and Line (Changed to black per your snippet)
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.15) +
  geom_line(color = "black", linewidth = 1.2) +
  # D. Significance Asterisks at y = 18
  geom_text(data = . %>% filter(p_val < 0.05),
            aes(label = "*", y = -8), 
            color = "black", size = 6) +
  # E. Formatting & Scales
  theme_classic() +
  scale_x_continuous(breaks = c(0,20, 40, 49,60, 75)) +
  # USE coord_cartesian TO ZOOM WITHOUT REMOVING DATA
  coord_cartesian(ylim = c(-10, 25)) + 
  # F. Labels and Styling
  labs(title = "Session HRV Mean Difference (SKY − Control)",
       subtitle = "Values > 0 indicate SKY > Control | * p < 0.05",
       x = "Time Since Start of Treatment (minutes)", 
       y = "HRV Estimate") +
  theme(
    legend.position = c(0.10, 0.55),
    plot.title = element_text(face = "bold", size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    # Optional: If you wanted to change the font size of the axis labels:
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

print(FigA)

ggsave("FigA_Session_HRV.pdf", FigA, width = 6, height = 4)


# ==============================================================================
# Figure B: Sleep HRV Mean Difference (SKY − Control)
# ==============================================================================
time_levels_weekly <- c("Baseline", "Week 1", "Week 2", "Week 3", 
                        "Week 4", "Week 5", "Week 6", "Week 7", "Week 8")

hrv_diff_clean <- hrv_clean %>%
  # Filter for difference group and exclude 'Day' rows
  filter(treatment == "Sky-Control") %>%
  filter(time %in% time_levels_weekly) %>%
  # Set factor levels and convert to numeric 0-8
  mutate(time = factor(time, levels = time_levels_weekly),
         time_num = as.numeric(time) - 1) %>%
  rename(p_val = p) 

# IDENTIFY SIGNIFICANT POINTS
sig_times <- hrv_diff_clean %>%
  filter(p_val < 0.05) %>%
  pull(time_num)

# CREATE THE DIFFERENCE PLOT
FigB <- hrv_diff_clean %>%
  ggplot(aes(x = time_num, y = Estimate)) +
  
  # A. Vertical Significance Lines
  geom_vline(xintercept = sig_times,
             linetype = "dotted", color = "gray80", linewidth = 0.5) +
  
  # B. Horizontal Reference Line at Zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  
  # C. Difference Ribbon and Line
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.15) +
  geom_line(color = "black", linewidth = 1.2) +
  
  # D. ADD BLACK DOTS AT EACH WEEK
  geom_point(color = "black", size = 2) +
  
  # E. Significance Asterisks 
  geom_text(data = . %>% filter(p_val < 0.05),
            aes(label = "*", y = -8), 
            color = "black", size = 8) +
  
  # F. Formatting & Scales
  theme_classic() +
  scale_x_continuous(breaks = 0:8, labels = time_levels_weekly) +
  coord_cartesian(ylim = c(-15, 20)) + 
  
  # G. Labels and Styling
  labs(title = "Sleep HRV Mean Difference (SKY − Control)",
       # subtitle = "Weekly Analysis with Observed Data Points",
       x = "Study Time in Weeks", 
       y = "HRV Estimate") +
  theme(plot.title = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        # Optional: If you wanted to change the font size of the axis labels:
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12))

print(FigB)

ggsave("Fig4B_Sleep_HRV.pdf", FigB, width = 6, height = 4)




# ==============================================================================
# Figure C: Sleep HRV Mean Difference (SKY − Control)
# Elevated baseline inflammation subgroup
# ==============================================================================
time_map <- c("Baseline"="Baseline", "FU: 1"="Week 1", "FU: 2"="Week 2", "FU: 3"="Week 3", 
              "FU: 4"="Week 4", "FU: 5"="Week 5", "FU: 6"="Week 6", "FU: 7"="Week 7", "FU: 8"="Week 8")

df_elevated <- sleep_elevate %>%
  filter(treatment == "SKY-Control" & inflam_grp == "Elevated") %>%
  filter(time %in% names(time_map)) %>%
  mutate(time_label = factor(time_map[time], levels = unname(time_map)),
         time_num = as.numeric(time_label) - 1)

sig_points <- df_elevated %>% filter(p < 0.05)

Fig4C <- ggplot(df_elevated, aes(x = time_num, y = estimate)) +
  # A. Reference Lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  
  # B. Ribbon and Trend Line
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.15) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_vline(data = . %>% filter(p < 0.05),
             aes(xintercept = time_num),
             linetype = "dotted", color = "gray80", linewidth = 0.5) +
  # C. ADD BLACK DOTS AT EACH WEEK
  geom_point(color = "black", size = 2) +
  geom_text(data = . %>% filter(p < 0.05),
            aes(x = time_num, y = -28, label = "*"), 
            color = "black", size = 10, fontface = "bold") +
  
  # E. Formatting & Scales
  scale_x_continuous(breaks = 0:8, labels = unname(time_map)) +
  coord_cartesian(ylim = c(-30, 30)) +
  theme_classic() +
  labs(title = "Sleep HRV Mean Difference (SKY − Control) \nElevated baseline inflammation subgroup",
       subtitle = "Difference (SKY - Control) | * p < 0.05",
       x = "Study Time in Weeks", 
       y = "HRV Estimate") +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        # Optional: If you wanted to change the font size of the axis labels:
        axis.text.x = element_text(size = 12, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 12))
Fig4C

ggsave("Fig4C_Elevated_inflammation_group.pdf", Fig4C, width = 6, height = 4)



# ==============================================================================
# Figure D: Sleep HRV Mean Difference (SKY − Control)
# Low baseline inflammation subgroup
# ==============================================================================
time_map <- c("Baseline"="Baseline", "FU: 1"="Week 1", "FU: 2"="Week 2", "FU: 3"="Week 3", 
              "FU: 4"="Week 4", "FU: 5"="Week 5", "FU: 6"="Week 6", "FU: 7"="Week 7", "FU: 8"="Week 8")

df_low <- sleep_elevate %>%
  filter(treatment == "SKY-Control" & inflam_grp == "Low") %>%
  filter(time %in% names(time_map)) %>%
  mutate(time_label = factor(time_map[time], levels = unname(time_map)),
         time_num = as.numeric(time_label) - 1)

FigD<-ggplot(df_low, aes(x = time_num, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.15) +
  geom_line(color = "black", linewidth = 1.2) +
  scale_x_continuous(breaks = 0:8, labels = unname(time_map)) +
  geom_point(color = "black", size = 2) +
  coord_cartesian(ylim = c(-35, 30)) +
  theme_classic() +
  labs(title = "Mean HRV Difference  (SKY − Control): \nLow baseline inflammation subgroup",
       # subtitle = "Difference (SKY - Control) | No Significant Timepoints",
       x = "Study Time in Weeks", y = "HRV Estimate") +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        # Optional: If you wanted to change the font size of the axis labels:
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12))

Fig4D

ggsave("Fig4D_Low_inflamarion_group.pdf", 
       Fig4D, width = 6, height = 4)
