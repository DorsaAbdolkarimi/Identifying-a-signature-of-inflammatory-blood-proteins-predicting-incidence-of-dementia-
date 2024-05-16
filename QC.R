
which(colnames(data) == "ML_PROT_f30900_A1BG_1_0") 

complete_prot <-  colSums(!is.na(data[,17:2935])) ##this is how many measurement each protein has 
write.csv(complete_prot,"test.csv")

missing_rows <- rowSums(is.na(data[,17:2935])) ##there is 2919 proteins
write.csv(missing_rows,"test_rows.csv")
sum(missing_rows == 2919) ##449128 people have 2919 proteins missing = this is all proteins so means they don't have proteomics data 
502356 - 449128

##eliminating individuals wo proteomics data 
data_proteomoics <- data[missing_rows != 2919, ]

##now out of these poeple with proteomics data, see how many NAs each person has 
missing_rows <- rowSums(is.na( data_proteomoics[,17:2935]))
0.80*2919

##how many individuals with more than 10% missing 
missing_10 <- sum( missing_rows < 291.9) ## 41364 out of 52,164 have measurements for at least 90% of the proteins 
missing_20 <- sum( missing_rows < 583.8) ## 44150 out of 52,164 have measurements for at least 80% of the proteins 

52164*0.9
52164*0.8
prot_10 <- sum(complete_prot >= 46947.6) ##1457 out of 2919 proteins have less than 10% of people missing 
prot_20 <- sum(complete_prot >= 41731.2) ##2913 out of 2919 proteins have less than 20% of people missing 

##Exclude individuals with missing data in more than 20% of  proteins -!!! 
data_80 <-data[which(rowMeans(is.na(data[,17:2935])) <0.2),] ##only kept the 44190 people 
dim(data_80)

data_90 <- data[which(rowMeans(is.na(data[,17:2935])) <0.1),] ##only kept the 41400 people 
dim(data_90)

##Exclude metabolites with missing data in more than 20% of metabolites####
cols_to_remove <- c("ML_PROT_f30900_GLIPR1_1173_0",
                    "ML_PROT_f30900_NPM1_1889_0",
                    "ML_PROT_f30900_PCOLCE_1991_0",
                    "ML_PROT_f30900_CST1_709_0",
                    "ML_PROT_f30900_CTSS_732_0",
                    "ML_PROT_f30900_ENDOU_923_0")

# Remove the specified columns from the dataframe
data_80_80 <- data_80[, !names(data_80) %in% cols_to_remove]

data_80_80 <-data_80[,which(colMeans(is.na(data_80[, 17:2935])) < 0.2)] ##80% rows 80% cols ##thsi code hadn't deletecd the col 
dim(data_80_80)

data_90_80 <- data_90[,which(colMeans(is.na(data_90[, 17:2935])) < 0.2)]



# imputation  -------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
library(impute)
which(colnames(data_80_80) == "ML_PROT_f30900_A1BG_1_0")
i_prot_80 <- data_80_80[, 17:2929]
i_prot_80m <- as.matrix(i_prot_80)
meta_data_80 <- data_80_80[, 1:16]


# Perform KNN imputation
set.seed(2)
KNN_r <- impute.knn(i_prot_80m, k = 10)
KNN_r <- as.data.frame(KNN_r$data)
str(KNN_r)
imp_data <- cbind(meta_data_80, KNN_r) ##now there are no NAs in the data 
str(imp_data)


# selecting control cases -------------------------------------------------
library(MatchIt)

set.seed(2)
DementiaFMatches10 <-  matchit(incident_dem ~ sex + age + ethn + center,  data = imp_data, method = "nearest", ratio = 10, na.rm = TRUE)
summary(DementiaFMatches10)
matched_data <- match.data(DementiaFMatches10)
str(matched_data)
matched_data <- as.data.frame(matched_data)
##now the data is cleared from variables/obs with more than 20% NAs, imputaed for NAs and 1:10 case controls are selected based on similarities



# extracting only the inflammatory proteins data  -------------------------
##1. extract only protein names 
inflam_data <- matched_data 
colnames(inflam_data)[17:2929] <- str_extract(colnames(inflam_data)[17:2929], "(?<=ML_PROT_f30900_)[^_]+")
str(inflam_data)
which(colnames(inflam_data) =="ZPR1") ##last protein is on col 2929

##2. match with data from sun et al 
sun_data <- read_csv("sun_data.csv")
common_protein_names <- intersect(colnames(data_en)[17:2929], sun_data$protein) 

inflam_data <- inflam_data %>% select(all_of(common_protein_names)) ##only extracted proteins in common
var <-matched_data[, 1:16] ##extracting metadata 
inflam_data <- cbind(var, inflam_data)
str(inflam_data) ##dataframe now contains only inflammatory proteins and the metadata 


