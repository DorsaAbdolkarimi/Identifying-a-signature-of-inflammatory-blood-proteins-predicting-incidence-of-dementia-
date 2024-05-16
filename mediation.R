en_prot <-c("age", "sex", "apoe","time_dem", "incident_dem",
            "PLAUR",
            "CD276",
            "C7",
            "TCN1",
            "FN1", 
            "SFRP4",
            "ITGA11",
            "SERPINA5",
            "LCAT",
            "APOE", 
            "CLEC3B",
            "TNFSF10",
            "PEPD", 
            "MENT",
            "APOA2",
            "TNFRSF11B",
            "TFF2",
            "DDX39A",
            "EDAR",
            "SERPINI1",
            "KLKB1",
            "ITGA2",
            "FN1",
            "CLEC3B",
            "SERPINA5")

data_med <- complete_apoe[, names(complete_apoe) %in% en_prot]
complete.cases(data_med)
if (any(is.na(data_med))) {
  print("There are NA values in data_med")
} else {
  print("There are no NA values in data_med")
}
#create lists that hold metabolites and snps
prot <- colnames(data_med[, 6:27])
matrix_res <- matrix(data = NA, nrow = length(prot), ncol = 4)
rownames(matrix_res)<-prot
colnames(matrix_res) <- paste(c(paste("coef"), paste("lower.95"), paste("upper.95"), paste("Pr(>|z|)")), sep = "")


for (i in prot){
  str(i)
  linear_gwa<-lm(data_med[,i]~data_med$apoe, data=data_med)
  sum_gwa<-summary(linear_gwa) #save summaries table
  coef_gwa<-sum_gwa$coefficients #save coef table
  conf_gwa<-confint(linear_gwa, level=0.95) #save conf int table
  b_gwa<-coef_gwa[2,1]  #choose coef you want
  L95_gwa<-conf_gwa[2,1]
  U95_gwa<-conf_gwa[2,2]
  p_gwa<-coef_gwa[2,4]
  #add them in table
  matrix_res[i,1]<-b_gwa
  matrix_res[i,2]<-L95_gwa
  matrix_res[i,3]<-U95_gwa
  matrix_res[i,4]<-p_gwa
}

write.csv(matrix_res, paste("lm_prot_apoe_test.csv"), row.names=T)	#save tables
lm_test <- read_csv("lm_prot_apoe_test.csv") %>%
  as.data.frame()
str(lm_test)

prot <- colnames(data_med[, 6:27])
matrix_res <- matrix(data = NA, nrow = length(prot), ncol = 4)
rownames(matrix_res)<-prot
colnames(matrix_res) <- paste(c(paste("coef"), paste("lower.95"), paste("upper.95"), paste("Pr(>|z|)")), sep = "")



##dementia cox only adjusted for sex and age 
for (i in prot) {
  str(i)
  cox_model <- coxph(Surv(time_dem, incident_dem) ~ data_med[, i] + apoe, data = data_med)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[2, 1]  # Choose coef you want
  L95_gwa <- conf_gwa[2, 1]
  U95_gwa <- conf_gwa[2, 2]
  p_gwa <- coef_gwa[2, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}

write.csv(matrix_res, "med_step3_apoe_non exp.csv", row.names = TRUE)

Step1 <- coxph(Surv(time_dem, incident_dem) ~ apoe, data = data_med)

for (i in prot) {
  str(i)
  cox_model_apoe <- coxph(Surv(time_dem, incident_dem) ~ apoe, data = data_med)
  sum_gwa <- summary(cox_model) # Save summaries table
  coef_gwa <- sum_gwa$coefficients # Save coef table
  conf_gwa <- confint(cox_model, level = 0.95) # Save conf int table
  b_gwa <- coef_gwa[1, 1]  # Choose coef you want
  L95_gwa <- conf_gwa[1, 1]
  U95_gwa <- conf_gwa[1, 2]
  p_gwa <- coef_gwa[1, 5]
  # Add them in table
  matrix_res[i, 1] <- b_gwa
  matrix_res[i, 2] <- L95_gwa
  matrix_res[i, 3] <- U95_gwa
  matrix_res[i, 4] <- p_gwa
}
write.csv(matrix_res, "med_step1.csv", row.names = TRUE)


# bootstrapping -----------------------------------------------------------
library(boot)
#APOA2
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(APOA2 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + APOA2, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_apoa2 <- boot(data_med, mediation_function, R = 5000)
boot_mediation_apoa2
#Bootstrap Statistics :
#      original        bias    std. error
#t1* 0.01815311 5.430018e-05 0.005664328

boot.ci(boot.out = boot_mediation_apoa2, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   ( 0.0070,  0.0292 )   ( 0.0072,  0.0296 ) )  
#calculations and Intervals on Original Scale

plot(boot_mediation)
sd(boot_mediation$t) #0.005571648 which is the same as the function gave us earlier 

bias= mean(boot_mediation$t) - boot_mediation$t0
print(bias)

print(boot_mediation$t0-bias + c(-1, 1)*1.96*sd(boot_mediation$t))
#0.007268262 0.029109123 matched the result 


#TNFSF10
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(TNFSF10 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + TNFSF10, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_TNFSF10 <- boot(data_med, mediation_function, R = 5000)
boot_mediation_TNFSF10
#Bootstrap Statistics :
#      original        bias    std. error
#t1* 0.008750731 9.976226e-06 0.003105713

boot.ci(boot.out = boot_mediation_TNFSF10, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   ( 0.0027,  0.0148 )   ( 0.0035,  0.0156 ) 
#calculations and Intervals on Original Scale

plot(boot_mediation)
sd(boot_mediation$t) #0.005571648 which is the same as the function gave us earlier 

bias= mean(boot_mediation$t) - boot_mediation$t0
print(bias)

print(boot_mediation$t0-bias + c(-1, 1)*1.96*sd(boot_mediation$t))
#0.007268262 0.029109123 matched the result 


##APOE
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(APOE~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + APOE, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_APOE<- boot(data_med, mediation_function, R = 5000)
boot_mediation_APOE
#Bootstrap Statistics :
#      original        bias    std. error
#t1* 0.07448602 0.003638411  0.03910112

boot.ci(boot.out = boot_mediation_APOE, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   (-0.0058,  0.1475 )   ( 0.0093,  0.1623 )   
#calculations and Intervals on Original Scale


mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
for (i in en_prot){
  formula= as.formula(paste0(i," ~ apoe"))
  print(formula)}
result_a <- lm(formula, data= data_med)
a <- result_a$coefficients[2]
formula_2 = as.formula(paste0("Surv(time_dem, incident_dem) ~", i, "+ apoe"))
result_b <- coxph(formula_2, data = data_med)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_APOE <- boot(data_med, mediation_function, R = 5000)
boot_mediation_APOE
#Bootstrap Statistics :
#      original        bias    std. error
#t1* 0.07448602 0.003638411  0.03910112

boot.ci(boot.out = boot_mediation_APOE, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   (-0.0037,  0.1466 )   ( 0.0080,  0.1578 )   
#calculations and Intervals on Original Scale
plot(boot_mediation_APOE)

#CD276
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(CD276 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + CD276, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_CD276 <- boot(data_med, mediation_function, R = 5000)
boot_mediation_CD276
#Bootstrap Statistics :
#      original        bias    std. error
#0.006218289 -2.805695e-06  0.00271796

boot.ci(boot.out = boot_mediation_CD276, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   ( 0.0009,  0.0115 )   ( 0.0016,  0.0122 )     
#calculations and Intervals on Original Scale
plot(boot_mediation_CD276)

#FN1
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(FN1 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + FN1, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_FN1 <- boot(data_med, mediation_function, R = 5000)
boot_mediation_FN1
#Bootstrap Statistics :
#      original        bias    std. error
#t1* 0.006029991 7.049512e-06 0.002487566

boot.ci(boot.out = boot_mediation_FN1, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   ( 0.0011,  0.0109 )   ( 0.0018,  0.0116 )  
#calculations and Intervals on Original Scale


#SERPINA5
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(SERPINA5 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + SERPINA5, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_SERPINA5<- boot(data_med, mediation_function, R = 5000)
boot_mediation_SERPINA5

#Level      Normal             Percentile     
#t1* 0.005672572 4.292114e-06 0.002686823

boot.ci(boot.out = boot_mediation_SERPINA5, type = c("norm", "perc"))
#Level      Normal             Percentile     
#95%   ( 0.0004,  0.0109 )   ( 0.0012,  0.0116 )  

#TCN1
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(TCN1 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + TCN1, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_TCN1<- boot(data_med, mediation_function, R = 5000)
boot_mediation_TCN1
#Bootstrap Statistics :
#      original        bias    std. error
#t1* 0.003171416 6.373721e-06 0.001725149

boot.ci(boot.out = boot_mediation_TCN1, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   (-0.0002,  0.0065 )   ( 0.0000,  0.0068 )  
#calculations and Intervals on Original Scale


#SFRP4
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(SFRP4 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + SFRP4, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_SFRP4<- boot(data_med, mediation_function, R = 5000)
boot_mediation_SFRP4
#Bootstrap Statistics :
#      original        bias    std. error
#t1* 0.005264741 1.859105e-05 0.002325289

boot.ci(boot.out = boot_mediation_SFRP4, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   ( 0.0007,  0.0098 )   ( 0.0014,  0.0105 )  
#calculations and Intervals on Original Scale

#SERPINI1
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(SERPINI1 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + SERPINI1, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_SERPINI1<- boot(data_med, mediation_function, R = 5000)
boot_mediation_SERPINI1
#Bootstrap Statistics :
#      original        bias    std. error
#t1* 0.003650535 4.087116e-05  0.00203393

boot.ci(boot.out = boot_mediation_SERPINI1, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%    (-0.0004,  0.0076 )   ( 0.0004,  0.0082 )   
#calculations and Intervals on Original Scale

#KLKB1
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(KLKB1 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + KLKB1, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_KLKB1<- boot(data_med, mediation_function, R = 5000)
boot_mediation_KLKB1
#Bootstrap Statistics :
#      original        bias    std. error
#t1* 0.004303806 -2.159077e-05 0.002286915

boot.ci(boot.out = boot_mediation_KLKB1, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   (-0.0002,  0.0088 )   ( 0.0004,  0.0093 ) 
#calculations and Intervals on Original Scale

#ITGA2
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(ITGA2 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + ITGA2, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_ITGA2<- boot(data_med, mediation_function, R = 5000)
boot_mediation_ITGA2
#Bootstrap Statistics :
#      original        bias    std. error
#t1* -0.004345355 -7.327737e-06 0.002323541

boot.ci(boot.out = boot_mediation_ITGA2, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   (-0.0089,  0.0002 )   (-0.0093, -0.0002 )   
#calculations and Intervals on Original Scale

#ITGA11
set.seed(123)
mediation_function <- function(data_used, i)
{data_temp = data_med[i, ]
result_a <- lm(ITGA11 ~ apoe, data=data_temp)
a <- result_a$coefficients[2]
result_b <- coxph(Surv(time_dem, incident_dem) ~ apoe + ITGA11, data = data_temp)
b <- result_b$coefficients[2]
indirect_effect <- a*b
return(indirect_effect)}

boot_mediation_ITGA11<- boot(data_med, mediation_function, R = 5000)
boot_mediation_ITGA11
#Bootstrap Statistics :
#      original        bias    std. error
#t1* -0.005472432 4.320632e-05 0.002442138

boot.ci(boot.out = boot_mediation_ITGA11, type = c("norm", "perc"))
#Intervals :
#Level      Normal             Percentile
#95%   (-0.0103, -0.0007 )   (-0.0107, -0.0013 )  
#calculations and Intervals on Original Scale


# plots -------------------------------------------------------------------

library(ggplot2)
path_a <- ggplot(data=med, aes(x=prot, y=coef_a, ymin=lower_a, ymax=upper_a)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Protein") + ylab("Coefficient (95% CI)") +
  theme_bw() +  # use a white background
  ggtitle("Path a") +  # Add your title here
  theme(axis.title.x = element_text(size = 14),  # adjust x-axis label size
        axis.title.y = element_text(size = 14),  # adjust y-axis label size
        axis.text.x = element_text(size = 12),   # adjust x-axis text size
        axis.text.y = element_text(size = 12),   # adjust y-axis text size
        plot.title = element_text(size = 16))  # adjust title size
print(path_a)



path_b <- ggplot(data=med, aes(x=prot, y=coef_b, ymin=lower_b, ymax=upper_b)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Protein") + ylab("Coefficient (95% CI)") +
  theme_bw() +
  ggtitle("Path b") +  # Add your title here
  theme(axis.title.x = element_text(size = 14),  # adjust x-axis label size
        axis.title.y = element_text(size = 14),  # adjust y-axis label size
        axis.text.x = element_text(size = 12),   # adjust x-axis text size
        axis.text.y = element_text(size = 12),   # adjust y-axis text size
        plot.title = element_text(size = 16))  # adjust title size
print(path_b)

path_ct <- ggplot(data=med, aes(x=prot, y=coef_ct, ymin=lower_ct, ymax=upper_ct)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Protein") + ylab("Coefficient (95% CI)") +
  theme_bw() +
  ggtitle("Path c'") +  # Add your title here
  theme(axis.title.x = element_text(size = 14),  # adjust x-axis label size
        axis.title.y = element_text(size = 14),  # adjust y-axis label size
        axis.text.x = element_text(size = 12),   # adjust x-axis text size
        axis.text.y = element_text(size = 12),   # adjust y-axis text size
        plot.title = element_text(size = 16))  # adjust title size
print(path_ct)

path_ab <- ggplot(data=med, aes(x=prot, y=coef_ab, ymin=lower_ab, ymax=upper_ab)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Protein") + ylab("Coefficient (95% CI)") +
  theme_bw() +
  ggtitle("Path a x b") +  # Add your title here
  theme(axis.title.x = element_text(size = 14),  # adjust x-axis label size
        axis.title.y = element_text(size = 14),  # adjust y-axis label size
        axis.text.x = element_text(size = 12),   # adjust x-axis text size
        axis.text.y = element_text(size = 12),   # adjust y-axis text size
        plot.title = element_text(size = 16))  # adjust title size
print(path_ab)s