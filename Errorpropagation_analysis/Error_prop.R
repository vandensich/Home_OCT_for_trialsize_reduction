try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
source("functions_and_model.R")

cov_files <- list.files(pattern = "\\_matrix.csv$", full.names = FALSE)

# Remove the .cov extension from the filenames
file_names <- sub("\\_covariance_matrix.csv$", "", cov_files)

error_prop <- do.call(rbind,lapply(file_names, get_errors_from_file))
error_prop <- subset(error_prop, name == "CST_obs")


error_prop_MC <- do.call(rbind,lapply(file_names, get_errors_via_MC))
error_prop_MC_SoC <- error_prop_MC
error_prop$Method <- paste0(error_prop$Setting, " Delta Method")
error_prop_MC$Method <- paste0(error_prop_MC$Setting, " Monte Carlo Method")
save(error_prop_MC_SoC, file = "errorprop_MC_SoC.RData")

tot_error <- rbind(error_prop[,-6], error_prop_MC[,-c(1,6,7,8)])
tot_error_bettertrial <- tot_error
save(tot_error_bettertrial, file="tot_error_bettertrial.RData")

tot_error$UNI <- paste0(tot_error$Method,"_",tot_error$Drug)

library(splines)

spline_models_param <- tot_error %>%
  group_by(UNI) %>%
  do({
    model_se <- lm(Absolute_Error ~ bs(Number_of_Patients, df = 3), data = .)
    model_lower <- lm(Lower_CI ~ bs(Number_of_Patients, df = 3), data = .)
    model_upper <- lm(Upper_CI ~ bs(Number_of_Patients, df = 3), data = .)
    tibble(model_se = list(model_se), model_lower = list(model_lower), model_upper = list(model_upper))

  })

predict_for_setting <- function(setting_name, patient_values, model_type) {
  model <- spline_models_param %>%
    filter(UNI == setting_name) %>%
    pull(!!sym(model_type)) %>%
    .[[1]]
  predict(model, data.frame(Number_of_Patients = patient_values))
}

patient_values <- seq(10, 55, by = 1)

# Get all unique settings
all_settings <- unique(spline_models_param$UNI)

# Generate predictions for each setting
prediction_data <- lapply(all_settings, function(sett) {
  data.frame(
    Uni = sett,
    Number_of_Patients = patient_values,
    Absolute_Error_Pred = predict_for_setting(sett, patient_values, "model_se"),
    Lower_CI_Pred = predict_for_setting(sett, patient_values, "model_lower"),
    Upper_CI_Pred = predict_for_setting(sett, patient_values, "model_upper")
  )
}) %>%
  bind_rows()
prediction_data$Drug <- sub("^[^_]*_", "", prediction_data$Uni)
save(prediction_data, file = "error_pred.RData")

zscore_fit <- data.frame (
  Number_of_Patients = subset(prediction_data, Drug == "Eylea")$Number_of_Patients,
  Absolute_Error_Eylea = subset(prediction_data, Drug == "Eylea")$Absolute_Error_Pred,
  Absolute_Error_Trial = subset(prediction_data, Drug == "Trialdrug")$Absolute_Error_Pred,
  Lower_CI_Eylea = subset(prediction_data, Drug == "Eylea")$Lower_CI_Pred,
  Upper_CI_Eylea = subset(prediction_data, Drug == "Eylea")$Upper_CI_Pred,
  Lower_CI_Trial = subset(prediction_data, Drug == "Trialdrug")$Lower_CI_Pred,
  Upper_CI_Trial = subset(prediction_data, Drug == "Trialdrug")$Upper_CI_Pred,
  zscore = 47.9 / sqrt(subset(prediction_data, Drug == "Eylea")$Absolute_Error_Pred^2 + subset(prediction_data, Drug == "Trialdrug")$Absolute_Error_Pred^2),
  zscore_Lower = 47.9 / sqrt(subset(prediction_data, Drug == "Eylea")$Upper_CI_Pred^2 + subset(prediction_data, Drug == "Trialdrug")$Upper_CI_Pred^2),
  zscore_Upper = 47.9 / sqrt(subset(prediction_data, Drug == "Eylea")$Lower_CI_Pred^2 + subset(prediction_data, Drug == "Trialdrug")$Lower_CI_Pred^2),
  Setting = sub(" .*", "",subset(prediction_data, Drug=="Eylea")$Uni),
  Method = sub("^[^ ]* ", "", sub("_.*", "",subset(prediction_data, Drug=="Eylea")$Uni))
  
)

load(file= "zscore_data.RData")

zscore_fit <- subset(zscore_fit, !(Setting == "Notal" & Number_of_Patients > 40))

zscore <- rbind(zscore_fit, zscore_data)

zscore[zscore$Setting=="Notal",]$Setting <- "Home Monitoring"
zscore[zscore$Setting=="OnSite",]$Setting <- "OnSite Monitoring"


plot_data <- zscore %>%
  pivot_longer(
    cols = c(Absolute_Error_Eylea, Absolute_Error_Trial, zscore),
    names_to = "name",
    values_to = "value"
  ) %>%
  mutate(
    lower_CI = case_when(
      name == "Absolute_Error_Eylea" ~ Lower_CI_Eylea,
      name == "Absolute_Error_Trial" ~ Lower_CI_Trial,
      name == "zscore" ~ zscore_Lower
    ),
    upper_CI = case_when(
      name == "Absolute_Error_Eylea" ~ Upper_CI_Eylea,
      name == "Absolute_Error_Trial" ~ Upper_CI_Trial,
      name == "zscore" ~ zscore_Upper
    )
  )
significance <- ggplot(subset(plot_data, Number_of_Patients <= 55 & (Method== "Data" | Method == "Monte Carlo Method") ),
                       aes(x = Number_of_Patients, y = value, color = Setting)) +
  geom_line(aes(group = Setting), size = 1) +
  geom_ribbon(aes(ymin = lower_CI, ymax=upper_CI, color=Setting, fill=Setting), alpha=0.2)+
  # Add a horizontal line at y=1.96 only in the zscore panels
  geom_hline(aes(yintercept = ifelse(name == "zscore", 1.96, NA)),
             color = "orange")+
  geom_vline(aes(xintercept = ifelse(Setting == "Data", 47, NA)),
             color = "red", linetype = "dashed")+
  geom_vline(aes(xintercept = ifelse(Setting == "Data", 41, NA)),
             color = "lightcoral", linetype = "dashed")+
  geom_vline(aes(xintercept = ifelse(Setting == "Data", 54, NA)),
             color = "lightcoral", linetype = "dashed")+
  geom_vline(aes(xintercept = ifelse(Setting == "Home Monitoring", 34, NA)),
             color = "green", linetype = "dashed")+
  geom_vline(aes(xintercept = ifelse(Setting == "Home Monitoring", 33, NA)),
             color = "lightgreen", linetype = "dashed")+
  geom_vline(aes(xintercept = ifelse(Setting == "Home Monitoring", 35, NA)),
             color = "lightgreen", linetype = "dashed")+
  geom_vline(aes(xintercept = ifelse(Setting == "OnSite Monitoring", 52, NA)),
             color = "blue", linetype = "dashed")+
  geom_vline(aes(xintercept = ifelse(Setting == "OnSite Monitoring", 50, NA)),
             color = "lightblue", linetype = "dashed")+
  geom_vline(aes(xintercept = ifelse(Setting == "OnSite Monitoring", 54, NA)),
             color = "lightblue", linetype = "dashed")+
  labs(
    x = "Number of Patients",
    y = NULL,          # Remove the default "value" y-axis label
    color = "Setting"
  ) +
  theme_dMod() +
  
  facet_grid(
    name ~ Setting,
    scales = "free_y",
    switch = "y",               # Move facet strips to the left side
    labeller = labeller(
      name = function(x) {
        # Replace underscores with spaces
        x <- gsub("_", " ", x)
        
        # Rename "zscore" to "z-score"
        x[x == "zscore"] <- "z-score"
        
        # Add "[um]" to the first two labels (assumed to be the absolute error metrics)
        x[x == "Absolute Error Eylea"] <- "Absolute Error on \n Endpoint for Eylea [um]"
        x[x == "Absolute Error Trial"] <- "Absolute Error on \n Endpoint for Trialdrug [um]"
        
        x
      },
    )
  ) +
  
  theme(
    text = element_text(size = 16),
    legend.position = "bottom",
    strip.placement = "outside", # Ensure strips appear outside the panel
    plot.background = element_blank(),
    panel.background = element_blank(),
    strip.background =  element_blank()
    
  )
ggsave("zscore_lines.pdf",significance, height= 10.5, width= 12 )

