##cox en on all control with crossvalidation training for lambda 
library(caret)
library(glmnet)

cox_en_data <- inflam_data

cox_en_data$sex <- as.character(cox_en_data$sex)
cox_en_data$sex[cox_en_data$sex == "Female"] <- "1" #assigning these levels would have generated NAs if f and m were not character vectors 
cox_en_data$sex[cox_en_data$sex == "Male"] <- "2"

set.seed(123)
index <- createDataPartition(cox_en_data$incident_dem, p = 0.7, list = FALSE)
train_data <- cox_en_data[index, ]
test_data <- cox_en_data[-index, ]

##training data 
y_train <- Surv(train_data$time_dem , train_data$incident_dem) 
var_train <- data.frame(age = as.numeric(train_data$age), sex = as.numeric(train_data$sex))
##checking for NAs in metadata 
sapply(var_train, function(x) sum(is.na(x)))
prot_train <- train_data[, 17:728] 
x_train <- cbind(var_train, prot_train) %>%
  as.matrix()

##test data 
y_test <- Surv(test_data$time_dem, test_data$incident_dem) 
var_test <- data.frame(age = as.numeric(test_data$age), sex = as.numeric(test_data$sex))
prot_test <- test_data[, 17:728] 
x_test <- cbind(var_test, prot_test) %>%
  as.matrix() 


##crossvalidation for lambda 
set.seed(123)
en.model.cv <- cv.glmnet(x_train, y_train, family = 'cox', alpha = 0.5, type.measure = 'C')
plot(en.model.cv) ##plotting crossvalidation for lambda and C-index for predictive capacity 

print(en.model.cv$lambda.1se) ## 0.01854537
print(en.model.cv$lambda.min) ##0.006072773

# Apply model to test data 
pred.surv <- predict(en.model.cv ,type="response", 
                     newx = x_test, s = 'lambda.1se')

#predictive accuracy 
c_index_1 <- Cindex(pred.surv, y_test) ## 0.6764283
library(survcomp)
cindex.comp(c_index_1, c_index_2)
library(survcomp)

# Compare C-index values between two models
comparison <- compare(c_index_1, c_index_2)


# plot and verify survival predictions  -------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("survcomp")

set.seed(123)
en.model <- glmnet(x_train, y_train, alpha = 0.5, family = "cox" )


lp <- predict(en.model,
              newx= x_test,
              s= en.model.cv$lambda.1se,
              type ="link")


y_test <- test_data[, c("time_dem", "incident_dem")]
dat.test <- data.frame(y_test)
dat.test$pred_status <- ifelse(lp>0,"dementia","control")
fit.surv.all <- survfit(Surv(time_dem, incident_dem) ~ pred_status, 
                        data = dat.test)

ggsurvplot(fit.surv.all,conf.int = TRUE,
           ylim = c(0.7, 1))

## fold difference 
summary_fit <- summary(fit.surv.all)
# Calculate fold difference in survival probability
control_surv <- summary_fit$surv[summary_fit$strata == "pred_status=control"]
dementia_surv <- summary_fit$surv[summary_fit$strata == "pred_status=dementia"]

# Compute fold difference
fold_difference <- dementia_surv / control_surv
average_fold_difference <- mean(fold_difference)

# Print the average fold difference
print(paste("Average Fold Difference (Dementia / Control):", average_fold_difference))
##0.0.960759440477081

# exponentiated
exponentiated_average_fold_difference <- exp(average_fold_difference)
print(paste("Exponentiated Average Fold Difference (Dementia / Control):", exponentiated_average_fold_difference))
##2.61368065477734 fold difference 


##Nelson Aalen plot 
dat.test$pred_status <- as.vector(dat.test$pred_status)
summary_fit <- summary(fit.surv.all)


nelson_aalen_data <- data.frame(
  time = summary_fit$time,
  strata = summary_fit$strata,
  cumulative_hazard = summary_fit$cumhaz
)

# Create the plot
nelson_aalen_plot <- ggplot(data = nelson_aalen_data, aes(x = time, y = cumulative_hazard, color = strata)) +
  geom_step() +
  labs(
    title = "Model 1",
    x = "Time (years)",
    y = "Cumulative Hazard"
  ) +
  theme_bw() +
  scale_color_manual(values = c("pred_status=control" = "blue", "pred_status=dementia" = "red")) +
  ylim(0, 0.1379218)  # Adjust the y-axis limits

# Print the plot
print(nelson_aalen_plot)

##apply model for all inflam proteins and obtain scores 
y <- Surv(cox_en_data$time_dem , cox_en_data$incident_dem) 
var <- data.frame(age = as.numeric(cox_en_data$age), sex = as.numeric(cox_en_data$sex))
prot <- cox_en_data[, 17:728] 
x <- cbind(var, prot) %>%
  as.matrix()

set.seed(123)
en.model.all <- glmnet( x, y, alpha = 0.5, lambda = en.model.cv$lambda.1se, family = "cox")

cox_en_coef <- coef(en.model.all) 
data_frame <- as.data.frame(as.matrix(cox_en_coef))
write.csv(data_frame, "cox_en_coef.csv", row.names = TRUE)


# model 2 -----------------------------------------------------------------

complete_apoe <- inflam_data[!is.na(inflam_data$apoe), ]
cox_en_data <- complete_apoe

cox_en_data$sex <- as.character(cox_en_data$sex)
cox_en_data$sex[cox_en_data$sex == "Female"] <- "1" #assigning these levels would have generated NAs if f and m were not character vectors 
cox_en_data$sex[cox_en_data$sex == "Male"] <- "2"

set.seed(123)
index <- createDataPartition(cox_en_data$incident_dem, p = 0.7, list = FALSE)
train_data <- cox_en_data[index, ]
test_data <- cox_en_data[-index, ]

##training data 
y_train <- Surv(train_data$time_dem , train_data$incident_dem) 
var_train <- data.frame(age = as.numeric(train_data$age), sex = as.numeric(train_data$sex), apoe = as.numeric(train_data$apoe))
##checking for NAs in metadata 
sapply(var_train, function(x) sum(is.na(x)))
prot_train <- train_data[, 17:728] 
x_train <- cbind(var_train, prot_train) %>%
  as.matrix()

##test data 
y_test <- Surv(test_data$time_dem, test_data$incident_dem) 
var_test <- data.frame(age = as.numeric(test_data$age), sex = as.numeric(test_data$sex), apoe = as.numeric(test_data$apoe))
prot_test <- test_data[, 17:728] 
x_test <- cbind(var_test, prot_test) %>%
  as.matrix() 


##crossvalidation for lambda 
set.seed(123)
en.model.cv <- cv.glmnet(x_train, y_train, family = 'cox', alpha = 0.5, type.measure = 'C')
plot(en.model.cv) ##plotting crossvalidation for lambda and C-index for predictive capacity 

print(en.model.cv$lambda.1se) ## 0.02098809
print(en.model.cv$lambda.min) ##0.009971043

# Apply model to test data 
pred.surv.apoe <- predict(en.model.cv ,type="response", 
                          newx = x_test, s = 'lambda.1se')

#predictive accuracy 
c_index_2 <- Cindex(pred.surv.apoe, y_test) ## 0.6973593

# plot and verify survival predictions  -------------------------------
set.seed(123)
en.model.apoe <- glmnet(x_train, y_train, alpha = 0.5, family = "cox")


lp <- predict(en.model.apoe,
              newx= x_test,
              s= en.model.cv$lambda.1se,
              type ="link")


y_test <- test_data[, c("time_dem", "incident_dem")]
dat.test <- data.frame(y_test)
dat.test$pred_status <- ifelse(lp>0,"dementia","control")
fit.surv.apoe <- survfit(Surv(time_dem, incident_dem) ~ pred_status, 
                         data = dat.test)

ggsurvplot(fit.surv.apoe,conf.int = TRUE,
           ylim = c(0.7, 1))

## fold difference 
# Extract summary of survival probabilities at the maximum time point
summary_fit <- summary(fit.surv.apoe)
# Calculate fold difference in survival probability
control_surv <- summary_fit$surv[summary_fit$strata == "pred_status=control"]
dementia_surv <- summary_fit$surv[summary_fit$strata == "pred_status=dementia"]

# Compute fold difference
fold_difference <- dementia_surv / control_surv
average_fold_difference <- mean(fold_difference)

# Print the average fold difference
print(paste("Average Fold Difference (Dementia / Control):", average_fold_difference))
##0.961953143869645

# exponentiated
exponentiated_average_fold_difference <- exp(average_fold_difference)
print(paste("Exponentiated Average Fold Difference (Dementia / Control):", exponentiated_average_fold_difference))
##2.6168024771363 fold difference 



##Nelson Aalen plot 
dat.test$pred_status <- as.vector(dat.test$pred_status)
summary_fit <- summary(fit.surv.apoe)

nelson_aalen_data <- data.frame(
  time = summary_fit$time,
  strata = summary_fit$strata,
  cumulative_hazard = summary_fit$cumhaz
)

# Create the plot
nelson_aalen_plot <- ggplot(data = nelson_aalen_data, aes(x = time, y = cumulative_hazard, color = strata)) +
  geom_step() +
  labs(
    title = "Model 2",
    x = "Time (years)",
    y = "Cumulative Hazard"
  ) +
  theme_bw() +
  scale_color_manual(values = c("pred_status=control" = "blue", "pred_status=dementia" = "red")) +
  ylim(0, 0.1379218)  # Adjust the y-axis limits

# Print the plot
print(nelson_aalen_plot)


##determining the fold difference between the two 
summary_fit <- summary(fit.surv)
surv_prob_group1 <- summary(fit.surv)$surv[1, nrow(summary(fit.surv)$surv)]
surv_prob_group2 <- summary(fit.surv)$surv[2, nrow(summary(fit.surv)$surv)]

# Calculate the fold difference
fold_difference <- surv_prob_group1 / surv_prob_group2

# Print the fold difference
print(paste("Fold Difference:", round(fold_difference, digits = 2)))


##apply model for all inflam proteins and obtain scores 
y <- Surv(cox_en_data$time_dem , cox_en_data$incident_dem) 
var <- data.frame(age = as.numeric(cox_en_data$age), sex = as.numeric(cox_en_data$sex), apoe = as.numeric(cox_en_data$apoe))
prot <- cox_en_data[, 17:728] 
x <- cbind(var, prot) %>%
  as.matrix()

set.seed(123)
en.model.all <- glmnet( x, y, alpha = 0.5, lambda = en.model.cv$lambda.1se, family = "cox")

cox_en_coef <- coef(en.model.all) 
data_frame <- as.data.frame(as.matrix(cox_en_coef))
write.csv(data_frame, "cox_en_coef_model2.csv", row.names = TRUE)


# Control model  ----------------------------------------------------------------
complete_apoe <- inflam_data[!is.na(inflam_data$apoe), ]
cox_en_data <- complete_apoe

cox_en_data$sex <- as.character(cox_en_data$sex)
cox_en_data$sex[cox_en_data$sex == "Female"] <- "1" #assigning these levels would have generated NAs if f and m were not character vectors 
cox_en_data$sex[cox_en_data$sex == "Male"] <- "2"

set.seed(123)
index <- createDataPartition(cox_en_data$incident_dem, p = 0.7, list = FALSE)
train_data <- cox_en_data[index, ]
test_data <- cox_en_data[-index, ]

##training data 
y_train <- Surv(train_data$time_dem , train_data$incident_dem) 
var_train <- data.frame(age = as.numeric(train_data$age), sex = as.numeric(train_data$sex), apoe = as.numeric(train_data$apoe))
##checking for NAs in metadata 
x_train <- var_train %>%
  as.matrix()

##test data 
y_test <- Surv(test_data$time_dem, test_data$incident_dem) 
var_test <- data.frame(age = as.numeric(test_data$age), sex = as.numeric(test_data$sex), apoe = as.numeric(test_data$apoe))
x_test <- var_test %>%
  as.matrix() 


##crossvalidation for lambda 
set.seed(123)
en.model.cv <- cv.glmnet(x_train, y_train, family = 'cox', alpha = 0.5, type.measure = 'C')
plot(en.model.cv) ##plotting crossvalidation for lambda and C-index for predictive capacity 

print(en.model.cv$lambda.1se) ## 0.1229277
print(en.model.cv$lambda.min) ##0.1229277

# Apply model to test data 
pred.surv <- predict(en.model.cv ,type="response", 
                     newx = x_test, s = 'lambda.1se')

#predictive accuracy 
c_index_3 <- Cindex(pred.surv, y_test) ## 0.6649445
c_index_3

# plot and verify survival predictions  -------------------------------
set.seed(123)
en.model <- glmnet(x_train, y_train, alpha = 0.5, family = "cox" )


lp <- predict(en.model,
              newx= x_test,
              s= en.model.cv$lambda.1se,
              type ="link")


y_test <- test_data[, c("time_dem", "incident_dem")]
dat.test <- data.frame(y_test)
dat.test$pred_status <- ifelse(lp>0,"dementia","control")
fit.surv <- survfit(Surv(time_dem, incident_dem) ~ pred_status, 
                    data = dat.test)

ggsurvplot(fit.surv,conf.int = TRUE,
           ylim = c(0.7, 1))

##determining the fold difference between the two 
summary_fit <- summary(fit.surv)
surv_prob_group1 <- summary(fit.surv)$surv[1, nrow(summary(fit.surv)$surv)]
surv_prob_group2 <- summary(fit.surv)$surv[2, nrow(summary(fit.surv)$surv)]

# Calculate the fold difference
fold_difference <- surv_prob_group1 / surv_prob_group2

# Print the fold difference
print(paste("Fold Difference:", round(fold_difference, digits = 2)))


##apply model for all inflam proteins and obtain scores 
y <- Surv(cox_en_data$time_dem , cox_en_data$incident_dem) 
var <- data.frame(age = as.numeric(cox_en_data$age), sex = as.numeric(cox_en_data$sex), apoe = as.numeric(cox_en_data$apoe))
prot <- cox_en_data[, 17:728] 
x <- cbind(var, prot) %>%
  as.matrix()

set.seed(123)
en.model.all <- glmnet( x, y, alpha = 0.5, lambda = en.model.cv$lambda.1se, family = "cox")

cox_en_coef <- coef(en.model.all) 
data_frame <- as.data.frame(as.matrix(cox_en_coef))
write.csv(data_frame, "cox_en_coef_model2.csv", row.names = TRUE)