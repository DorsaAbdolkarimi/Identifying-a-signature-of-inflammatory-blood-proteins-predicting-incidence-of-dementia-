my_ukb_data<-readRDS("UKBDorsa_78867r674582_140224.rds")
data <- my_ukb_data
apoe <- data$APOE4_alleles 
age <- data$R_Age_f21003_0
sex <- data$R_Sex_f31
incident_dem <- data$R_ML_C42C240Xf41270f20002_Dementia_Incident
incident_AD <- data$R_ML_C42C240Xf41270_Alzheimers_Incident
time_dem <- data$ML_C42C240Xf41270f20002_DementiaDiagyr
time_AD <- data$ML_C42C240Xf41270_AlzheimersDiagyr
bmi <- data$R_BMI_f21001_0
edu <- data$R_ML_Qualification_CompositeFinal_SecNonSec_f845f613
ethn <- data$R_Ethnicity_f21000
diabetes <- data$ML_C42C240Xf41270f20002_DiabetesType2
smoking <- data$ever_smoked_f20160_0_0 ##ever smoked yes or no 
hypertension <- data$ML_C42C240Xf41270f20002_Hypertension
alcohol <- data$R_AlcStatus_f20117_0
egfr <- data$creatinine_f30700_0_0
center <- data$uk_biobank_assessment_centre_f54_0_0

which(colnames(data) == "ML_PROT_f30900_A1BG_1_0") ##[1] 1617  is the first protein col
proteins <- data[, 1617:4535]

meta.data <- data.frame(age = age, sex = sex, apoe = apoe, incident_AD = incident_AD, time_AD = time_AD, time_dem = time_dem,
                        incident_dem = incident_dem, bmi = bmi, edu = edu, ethn = ethn, diabetes = diabetes, smoking = smoking, 
                        hypertension = hypertension, alcohol = alcohol, egfr = egfr, center = center)

data <- cbind(meta.data, proteins)

table(data$incident_dem)

##organising time and incident for dem 
data$incident_dem[is.na(data$incident_dem)] <- 2 ##setting incident == 2 for controls to distinguish them from 0 == prevelent cases
##Removing prevelant incidents 
data <- subset(data, incident_dem != 0)

##now assigning zero to no dementia 
data$incident_dem[data$incident_dem== 2] <- 0 ## from here dementia has incident == 1 and no dementia incident == 0

##setting up the value for those who never developed dementia as the maximum time 
censored_indicents <- is.na(data$time_dem)
data$time_dem[censored_indicents] <- max(data$time_dem, na.rm = TRUE)
str(data$time_dem)

##organising time and incident for AD 
table(data$incident_AD)
data$incident_AD[is.na(data$incident_AD)] <- 0 ## Setting incident == 2 for controls to distinguish them from 0 == prevalent cases


## Setting up the value for those who never developed dementia as the maximum time 
censored_indices <- is.na(data$time_AD)
data$time_AD[censored_indices] <- max(data$time_AD, na.rm = TRUE)
str(data$time_AD)


# for the demographics table ---------------------------------------------
##step 1 and do the code below 
data_all  <- data
##step 2 
data_all <- data_proteomoics
##step 3 
data_all <- imp_data 
##step 4 
data_all <- matched_data


##organising time and incident for dem 
data_all$incident_dem[is.na(data_all$incident_dem)] <- 2 ##setting incident == 2 for controls to distinguish them from 0 == prevelent cases
data_all$incident_AD[is.na(data_all$incident_AD)] <- 2 ##setting incident == 2 for controls to distinguish them from 0 == prevelent cases

##age
mean(data_all$age[data_all$incident_dem== 1], na.rm = TRUE)
mean(data_all$age)

##sex
table(data_all$sex)
table(data_all$sex[data_all$incident_dem== 1])

##apoe
table(data_all$apoe)
table(data_all$apoe[data_all$incident_dem== 1])

##ethn 
table(data_all$ethn[data_all$incident_dem== 1])
table(data_all$ethn[data_all$incident_AD== 1])
table(data_all$ethn)

##BMI
mean(data_all$bmi[data_all$incident_dem== 1], na.rm = TRUE)
mean(data_all$bmi, na.rm = TRUE)

##smoking
if (!("2" %in% levels(data_all$smoking))) {
  data_all$smoking <- factor(data_all$smoking, levels = c(levels(data_all$smoking), "2"))
}
data_all$smoking[is.nan(data_all$incident_dem)] <- 2
table(data_all$smoking[data_all$incident_dem== 1])
table(data_all$smoking)

##alcohol 
table(data_all$alcohol[data_all$incident_dem== 1])
table(data_all$alcohol)

##hypertension 
table(data_all$hypertension[data_all$incident_dem== 1])
table(data_all$hypertension[data_all$incident_AD== 1])
table(data_all$hypertension)

##diabetes
table(data_all$diabetes)
table(data_all$diabetes[data_all$incident_dem== 1])

##education 
table(data_all$edu)
table(data_all$edu[data_all$incident_dem== 1])

##centre 
table(data_all$center)
table(data_all$center[data_all$incident_dem== 1])
table(data_all$center[data_all$incident_AD== 1])
