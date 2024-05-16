library(survival)

##checking position of first and last proteins 
names(inflam_data)[(ncol(inflam_data)-4):ncol(inflam_data)]
which(colnames(inflam_data) == "A1BG") ##col 17
which(colnames(inflam_data) =="ZNF174") ##col 744

##double checking there is no NA in time col 
if (any(is.na(inflam_data$time_dem))) {
  print("There are NA values in the time_dem column.")
} else {
  print("There are no NA values in the time_dem column.")
} ##"There are no NA values in the time_dem column." therefore all the censored events have been set to max time

##categorizing ethnicity 
inflam_data$ethn <- as.character(inflam_data$ethn)
inflam_data$ethn_group <- ifelse(inflam_data$ethn %in% c("White", "S.Asian", "Black"), inflam_data$ethn, "Other") %>%
  as.factor()

# dementia cox ------------------------------------------------------------
##prepring for cox without regularisation 
prot <- colnames(inflam_data[, 17:744])
matrix_res <- matrix(data = NA, nrow = length(prot), ncol = 4)
rownames(matrix_res) <- prot
colnames(matrix_res) <- paste(c(paste("coef"), paste("lower .95"), paste("upper .95"), paste("Pr(>|z|)")), sep = "")


## no protein adjustments to loook at apoe effect 
cox_model <- coxph(Surv(time_dem, incident_dem) ~ age + factor(inflam_data$sex) + apoe + PLAUR ,data = inflam_data, na.action = na.exclude)
cox_zph <- cox.zph(cox_model) ## to check redisuals 
cox_zph

##dementia cox only adjusted for sex and age 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_dem, incident_dem) ~ inflam_data[, i] + age + factor(inflam_data$sex), data = inflam_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

write.csv(matrix_res, "dem_inflam_cox_basic.csv", row.names = TRUE)

write.csv(matrix_res, paste("dem_inflam_cox_basic.csv"), row.names = TRUE) # Save tables


##dem cox for sex age apoe 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_dem, incident_dem) ~ inflam_data[, i] + age + factor(inflam_data$sex) + apoe, data = inflam_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

write.csv(matrix_res, "dem_inflam_cox_apoe.csv", row.names = TRUE)
write.csv(matrix_res, paste("dem_inflam_cox_apoe.csv"), row.names = TRUE) # Save tables

##dem cox full adjustment 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_dem, incident_dem) ~ inflam_data[, i] + age + factor(inflam_data$sex) + apoe + 
                       bmi + edu + diabetes  +smoking + hypertension + alcohol + egfr + ethn_group, data = inflam_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

write.csv(matrix_res, "dem_inflam_cox_full.csv", row.names = TRUE)
write.csv(matrix_res, paste("dem_inflam_cox_full.csv"), row.names = TRUE) # Save tables

## full adjustment without apoe 
##dem cox full adjustment 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_dem, incident_dem) ~ inflam_data[, i] + age + factor(inflam_data$sex)  + 
                       bmi + edu + diabetes  +smoking + hypertension + alcohol + egfr + ethn_group, data = inflam_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

data_res <- as.data.frame(matrix_res)
data_res$fdr <- p.adjust(data_res$`Pr(>|z|)`, method = "fdr")
write.csv(data_res, "dem_inflam_cox_full_wo_apoe.csv", row.names = TRUE)


# AD cox ------------------------------------------------------------------
##AD cox only adjusted for sex and age 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_AD, incident_AD) ~ inflam_data[, i] + age + factor(inflam_data$sex), data = inflam_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

write.csv(matrix_res, "AD_inflam_cox_basic.csv", row.names = TRUE)
write.csv(matrix_res, paste("AD_inflam_cox_basic.csv"), row.names = TRUE) # Save tables

##AD cox for sex age apoe 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_AD, incident_AD) ~ inflam_data[, i] + age + factor(inflam_data$sex) + apoe, data = inflam_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

write.csv(matrix_res, "AD_inflam_cox_apoe.csv", row.names = TRUE)
write.csv(matrix_res, paste("AD_inflam_cox_apoe.csv"), row.names = TRUE) # Save tables


##AD cox only adjusted for sex and age 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_AD, incident_AD) ~ inflam_data[, i] + age + factor(inflam_data$sex), data = inflam_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

write.csv(matrix_res, "AD_inflam_cox_basic.csv", row.names = TRUE)
write.csv(matrix_res, paste("AD_inflam_cox_basic.csv"), row.names = TRUE) # Save tables

##AD cox full adjustment 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_AD, incident_AD) ~ inflam_data[, i] + age + factor(inflam_data$sex) + apoe + 
                       bmi + edu + diabetes  +smoking + hypertension + alcohol + egfr + ethn_group, data = inflam_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

write.csv(matrix_res, "AD_inflam_cox_full.csv", row.names = TRUE)
write.csv(matrix_res, paste("AD_inflam_cox_full.csv"), row.names = TRUE) # Save tables


## full AD adjustment without apoe 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_AD, incident_AD) ~ inflam_data[, i] + age + factor(inflam_data$sex)  + 
                       bmi + edu + diabetes  +smoking + hypertension + alcohol + egfr + ethn_group, data = inflam_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

data_res <- as.data.frame(matrix_res)
data_res$fdr <- p.adjust(data_res$`Pr(>|z|)`, method = "fdr")
write.csv(data_res, "AD_inflam_cox_full_wo_apoe.csv", row.names = TRUE)






# applying FDR adjustments  -----------------------------------------------
library(readxl)
df <- read_excel("cox_for_fdr.xlsx")

# Get the column names
col_names <- names(df)

# Iterate over column names
for (col in col_names) {
  # Check if the column name ends with "_p"
  if (endsWith(col, "_p")) {
    # Perform FDR correction and add a new column for corrected p-values
    df[paste0(col, "_fdr")] <- p.adjust(df[[col]], method = "fdr")
  }
}

str(df)
write.csv(df, paste("cox_fdr.csv"), row.names = TRUE)

significant_data <- df[df$basic_dem_p_fdr < 0.05, ]
write.csv(significant_data, paste("sig_cox_fdr.csv"), row.names = TRUE)

# volcano plot ------------------------------------------------------------
# Load the required libraries
library(ggplot2)
library(ggrepel)

# Sort the dataframe by FDR values in ascending order
df <- df[order(df$basic_dem_p_fdr), ]

# Negative logarithm (base 10) of FDR values
df$log_fdr <- -log10(df$basic_dem_p_fdr)

# Filter top significant variables
top_significant <- head(df[df$basic_dem_p_fdr < 0.05, ], 10)


p <- ggplot(df, aes(x = basic_dem_coef, y = log_fdr, color = basic_dem_p_fdr < 0.05)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("red", "black"), name = "FDR < 0.05") +  # Change legend title
  geom_text_repel(data = top_significant, aes(label = prot), size = 5, show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue3") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "blue3") +  # Add vertical line at x = 1
  labs(x = "Hazard Ratio", y = "-log10(FDR)",  # Adjust axis labels
       title = "Model 1") +
  theme_minimal() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),  # Adjust plot margins
        axis.title = element_text(size = 16),  # Increase axis label size
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 18)) +  # Increase axis text size
  coord_cartesian(ylim = c(min(df$log_fdr), max(df$log_fdr) + 1))  # Expand y-axis
# Save plot with increased dimensions
ggsave("volcano_plot.png", plot = p, width = 10, height = 8, units = "in")

# volcano plot for apoe adjusted ------------------------------------------

# Sort the dataframe by FDR values in ascending order
df <- df[order(df$apoe_dem_p_fdr), ]

# Negative logarithm (base 10) of FDR values
df$log_fdr <- -log10(df$apoe_dem_p_fdr)

# Filter top significant variables
top_significant <- head(df[df$apoe_dem_p_fdr < 0.05, ], 10)



p <- ggplot(df, aes(x = apoe_dem_coef, y = log_fdr, color = apoe_dem_p_fdr < 0.05)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("red", "black"), name = "FDR < 0.05") +  # Change legend title
  geom_text_repel(data = top_significant, aes(label = prot), size = 5, show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue3") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "blue3") +  # Add vertical line at x = 1
  labs(x = "Hazard Ratio", y = "-log10(FDR)",  # Adjust axis labels
       title = "Model 2") +
  theme_minimal() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),  # Adjust plot margins
        axis.title = element_text(size = 16),  # Increase axis label size
        axis.text = element_text(size = 16), 
        plot.title = element_text(size = 18)) +  # Increase axis text size
  coord_cartesian(ylim = c(min(df$log_fdr), max(df$log_fdr) + 1))  # Expand y-axis
ggsave("volcano_plot_apoe.png", plot = p, width = 10, height = 8, units = "in")



# sensitivity analysis  ---------------------------------------------------
sum(is.na(inflam_data$apoe)) ##113 are missing apoe genotype 
sum(is.na(inflam_data$apoe) & inflam_data$incident_dem == 1) ## 12 of them are those who have dementia 

sens_data <- inflam_data[!is.na(inflam_data$apoe), ]


##dem model 1 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_dem, incident_dem) ~ sens_data[, i] + age + factor(sens_data$sex), data = sens_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

data_res <- as.data.frame(matrix_res)
data_res$fdr <- p.adjust(data_res$`Pr(>|z|)`, method = "fdr")
write.csv(data_res, "sens_dem_inflam_cox_basic.csv", row.names = TRUE)


##dem model 2 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_dem, incident_dem) ~ sens_data[, i] + age + factor(sens_data$sex) + apoe, data = sens_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

data_res <- as.data.frame(matrix_res)
data_res$fdr <- p.adjust(data_res$`Pr(>|z|)`, method = "fdr")
write.csv(data_res, "sens_dem_inflam_cox_apoe.csv", row.names = TRUE)

##dem model 3 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_dem, incident_dem) ~ sens_data[, i] + age + factor(sens_data$sex)  + apoe +
                       bmi + edu + diabetes  +smoking + hypertension + alcohol + egfr + ethn_group
                     , data = sens_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

data_res <- as.data.frame(matrix_res)
data_res$fdr <- p.adjust(data_res$`Pr(>|z|)`, method = "fdr")
write.csv(data_res, "sens_dem_inflam_cox_full.csv", row.names = TRUE)

##model 4 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_dem, incident_dem) ~ sens_data[, i] + age + factor(sens_data$sex) +
                       bmi + edu + diabetes  +smoking + hypertension + alcohol + egfr + ethn_group
                     , data = sens_data, na.action = na.exclude)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 2]  # Choose coef you want
  L95_gwa <- exp(conf_gwa[1, 1])
  U95_gwa <- exp(conf_gwa[1, 2])
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

data_res <- as.data.frame(matrix_res)
data_res$fdr <- p.adjust(data_res$`Pr(>|z|)`, method = "fdr")
write.csv(data_res, "sens_dem_inflam_cox_full_wo_apoe.csv", row.names = TRUE)s