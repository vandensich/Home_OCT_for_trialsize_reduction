try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

library(dplyr)
library(ggplot2)

# Example reading CSV files (adjust paths as needed)
df_soc_OS   <- read.csv("Data_patients_OnSite_more_SoC.csv", stringsAsFactors = FALSE)
df_soc_OS$Setting <- "OnSite Monitoring"
df_trial_OS <- read.csv("Data_patients_OnSite_more_Trialdrug.csv", stringsAsFactors = FALSE)
df_trial_OS$Setting <- "OnSite Monitoring"

df_soc_NO   <- read.csv("Data_patients_Notal_SoC.csv", stringsAsFactors = FALSE)
df_soc_NO$Setting <- "Home monitoring"
df_trial_NO <- read.csv("Data_patients_Notal_Trialdrug.csv", stringsAsFactors = FALSE)
df_trial_NO$Setting <- "Home monitoring"

df_soc <- rbind(df_soc_NO, df_soc_OS)
df_trial <- rbind(df_trial_NO, df_trial_OS)


df_soc_sub   <- df_soc   %>% subset(name=="CST_obs"& time>=140)
df_trial_sub <- df_trial %>% subset(name=="CST_obs"& time>=140)


# Tag each subset with group labels so we can color them
df_soc_sub$drug   <- "Eylea"
df_trial_sub$drug <- "Trialdrug"

# Combine for convenience in plotting lines per patient
df_combined <- bind_rows(df_soc_sub, df_trial_sub)
df_combined$time <- df_combined$time - 140



# Plot
a <- ggplot(df_combined, aes(x = time, y = value)) +
  facet_wrap(.~Setting, ncol=2)+
  # Individual patient lines (lighter + transparent)
  # geom_line(aes(color = group), alpha = 0.3) +
  # geom_point(aes(color = drug), alpha = 0.2) +
  geom_line(aes(color = drug, group = interaction(drug,ID)), alpha = 0.1) +
  # Overlay the group means (solid, darker)
  stat_summary(
    data = subset(df_combined, drug == "Eylea"),
    fun = mean, 
    aes(x = time, y = value),
    geom = "line",
    size = 1.2, 
    color = "blue"      # a clear, darker red for SoC
  ) +
  stat_summary(
    data = subset(df_combined, drug == "Trialdrug"),
    fun = mean, 
    aes(x = time, y = value),
    geom = "line",
    size = 1.2, 
    color = "Orange"     # a clear, darker blue for Trialdrug
  ) +
  
  # Labels and theme
  labs(
    x = "Time",
    y = "CST in \u00B5m"
  ) +
  theme_minimal()+
  scale_color_manual(
    values = c("Eylea" = "blue", "Trialdrug" = "orange"),
    name   = "Drug"
  ) +
  scale_x_continuous(name = "Time [weeks]", breaks=c(0,56,112,168, 224),
                     labels=c("20", "28", "36","44", "52"))+
  theme(
    text = element_text(size = 20),
    
    legend.position = "none",
    strip.placement = "outside", # Ensure strips appear outside the panel
    plot.background = element_blank(),
    panel.background = element_blank(),
    strip.background =  element_blank(),
    strip.text.x = element_text(size = 20)
    
  )+
  geom_vline(xintercept = 168, color="green")+
  scale_y_log10(limits=c(90,1100))
ggsave("plot_all.jpg", a, width= 11.7, height= 8.3)
