library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(openxlsx)
#The following codes contain two parts, analysis on 1)all records 2) BCG-removed records.




###Analysis on all records
setwd("/Users/lidiayung/Downloads/Fwd_ MULTIR_ BC peptidomics datasets")
clin_raw <- read.xlsx("BC_CE-MS peptidomics_Metadata.xlsx")
peptides_raw <-read.xlsx("TransBioBC_08.04.2024_MosaID_1_7_5_MFinder_vs_MV_HybridSolution_v5_ML1_seq_md.xlsx",startRow=2)
# check if object has zero times; our min>0 
min(clin_raw$Follow.up.till.death)>0

##Sub-grouping
#Two groups of stages(MIBC:T2+= or NMIBC:T1)
clin_raw$Stage_baseline_MN <- ifelse(grepl("^pT2",clin_raw$Stage_grade_baseline),"MIBC","NMIBC")
#Two groups of stages(MIBC:T2+= or NMIBC:T1)
clin_raw$Stage_followup_MN <- ifelse(grepl("^pT2|^pT3|^pT4",clin_raw$Stage_grade_followup),"MIBC","NMIBC")
#Treatment groups, BCG, epirubicin, MMC
clin_raw$Treatment_groups <- ifelse(grepl("MMC", clin_raw$Treatment) & grepl("epirubicin", clin_raw$Treatment), "e-M Combined",
                                    ifelse(grepl("BCG", clin_raw$Treatment), "BCG",
                                           ifelse(grepl("MMC", clin_raw$Treatment), "MMC",
                                                  ifelse(grepl("epirubicin", clin_raw$Treatment), "epirubicin", "NA/none"))))
clin_raw$localisation_groups <- ifelse(grepl("\\bwall\\b", clin_raw$localisation), "bladder walls",
                                       ifelse(grepl("0", clin_raw$localisation), "0",
                                              ifelse(grepl("unknown", clin_raw$localisation), "unknown", "bladder neck, ostium and urethra")))

clin_raw$cytology_group <- ifelse(grepl("high", clin_raw$Cytology), "high_(grade 3)","low/others")

surv_obj_raw <- Surv(time = clin_raw$Follow.up.till.death, 
                     event = clin_raw$`Death.(1/.0)`=="1")

fit <- survfit(surv_obj_raw ~ 1, data = clin_raw)
ggsurvplot(fit, data = clin_raw, xlab = "Day", ylab = "Overall survival",risk.table = TRUE)
# Fit the survival curve
fit <- survfit(surv_obj_raw ~ 1, data = clin_raw)

# Plot the survival curve
surv_plot <- ggsurvplot(fit, data = clin_raw, xlab = "Time", ylab = "Overall survival", risk.table = TRUE)

# Calculate median survival time
median_time <- median(fit$time)

# Add a vertical line at the median time point
surv_plot <- surv_plot + geom_vline(xintercept = median_time, linetype = "dashed")
surv_plot
# Print the plot
print(surv_plot)

# Check the class and structure of median_time
class(median_time)
str(median_time)

# Check the class and structure of surv_plot
class(surv_plot)
str(surv_plot)



#KM curve ~ Variable of interest(Treatment_groups, gender, tumor_size_group, etc)
fit <- survfit(Surv(time = clin_raw$Follow.up.till.death, 
                    event = clin_raw$`Death.(1/.0)`=="1") ~ Gender, data = clin_raw)
ggsurvplot(fit, data = clin_raw, pval = TRUE,xlab = "Time", ylab = "Overall survival",risk.table = TRUE)

#Hazard ratio
fit.coxph<- coxph(surv_obj_raw ~ Gender+Stage_baseline_MN+Stage_followup_MN+cytology_group, data = clin_raw)
ggforest(fit.coxph, data = clin_raw)

library("glmpath")
library("glmnet")
library("penalized")

#Retrieve coefficients from the whole dataset
fit_glm <- glmnet(peptides_raw,surv_obj_raw,family="cox")# standardize = TRUE, maxit = 1000)
print(fit_glm)

cfs = coef(fit_glm,s=0.017530) #selecting 5 peptides 0.085230; 22 0.017530

meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
meaning_vals
meaning_coefs
coef_data <-data.frame(peptides=meaning_coefs,coefficient=meaning_vals)
#write.csv(coef_data,"/Users/lidiayung/Downloads/coefs_5peptides.csv")
#
coef_data <- coef_data[!(row.names(coef_data) %in% c("x99904775", "x99907555", "x99909218", "x99910097", "x99919757")), ]
sorted_coef_abs <- coef_data[order(-abs(coef_data$coefficient)), ]

length(meaning_vals)
top <- head(sorted_coef_abs,length(meaning_vals))
list(top$variable)

# Assuming top is a data frame containing the top coefficients
paste(top$variable, collapse = ",")

# Create bar plot
ggplot(coef_data, aes(x = reorder(peptides, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Peptides", y = "Coefficient") +
  theme_minimal()

##DPD hazard ratio 
rownames(peptides_raw) <- peptides_raw$`CE-MS.Analysis.ID`
peptides_subset <- peptides_raw[,coef_data$peptides]
length(peptides_subset)
#dot product
row.names(peptides_subset)
colnames(coef_data)

coef_data_subset <-coef_data[,-1]
coef_data_subset
result <- as.matrix(peptides_subset) %*% as.matrix(coef_data_subset)

#assign result to DPD_prognosis
clin_raw$DPD_prognosis <- result[row.names(result) %in% clin_raw$`CE-MS.Analysis.ID`, ]

clin_raw$outcome <- ifelse(clin_raw$DPD_prognosis > 0, "positive", "negative;0")
                           #ifelse(clin_raw$DPD_prognosis == 0, "0", "negative"))

#KM curve ~ Variable of interest(DPD outcome)
fit <- survfit(Surv(time = clin_raw$Follow.up.till.death, 
                    event = clin_raw$`Death.(1/.0)`=="1") ~ outcome, data = clin_raw)


ggsurvplot(fit, data = clin_raw, pval = TRUE,xlab = "Time", ylab = "Overall survival",risk.table = TRUE)

#HR calculation
fit.coxph<- coxph(surv_obj_raw ~ Gender+Stage_baseline_MN+Stage_followup_MN+cytology_group+outcome, data = clin_raw,iter.max=30)
ggforest(fit.coxph, data = clin_raw)

fit.coxph
###Removed BCG patients(3 records removed)
clin <-clin_raw[!grepl('^BCG', clin_raw$Treatment), ]

peptides <-peptides_raw[peptides_raw$`CE-MS.Analysis.ID` %in% clin$`CE-MS.Analysis.ID`, ]

##Sub-grouping
#Two groups of stages(MIBC:T2+= or NMIBC:T1)
clin$Stage_baseline_MN <- ifelse(grepl("^pT2",clin$Stage_grade_baseline),"MIBC","NMIBC")

#Two groups of stages(MIBC:T2+= or NMIBC:T1)
clin$Stage_followup_MN <- ifelse(grepl("^pT2|^pT3|^pT4",clin$Stage_grade_followup),"MIBC","NMIBC")

#Treatment groups, BCG, epirubicin, MMC
clin$Treatment_groups <- ifelse(grepl("MMC", clin$Treatment) & grepl("epirubicin", clin$Treatment), "e-M Combined",
                                ifelse(grepl("BCG", clin$Treatment), "BCG",
                                       ifelse(grepl("MMC", clin$Treatment), "MMC",
                                              ifelse(grepl("epirubicin", clin$Treatment), "epirubicin", "NA/none"))))
#localisation
clin$localisation_groups <- ifelse(grepl("\\bwall\\b", clin$localisation), "bladder walls",
                                   ifelse(grepl("0", clin$localisation), "0",
                                          ifelse(grepl("unknown", clin$localisation), "unknown", "bladder neck, ostium and urethra")))
#cytology
clin$cytology_group <- ifelse(grepl("high", clin$Cytology), "high_(grade 3)","low/others")
#ifelse(grepl("low", clin$Cytology), "low_(grade 2/1)",
#       ifelse(grepl("suspect (grade 1)", clin$Cytology), "low (grade 2/1)", "0/no cytology/no abnormalities/unknown")))

#tumor size
clin$tumor_size_group <- ifelse(grepl("0", clin$tumor.sizegroups), "0",
                                ifelse(grepl("unknown", clin$tumor.sizegroups), "unknown",
                                       ifelse(grepl("1-3",clin$tumor.sizegroups),"medium",
                                              ifelse(grepl("<1",clin$tumor.sizegroups),"small","large"))))

##KM curves 
#Overall survival curve
surv_obj <- Surv(time = clin$Follow.up.till.death, 
                 event = clin$`Death.(1/.0)`=="1")
fit <- survfit(surv_obj ~ 1, data = clin)
ggsurvplot(fit, data = clin, xlab = "Time", ylab = "Overall survival",risk.table = TRUE)

#KM curve surv_obj~ Variable of interest(Treatment_groups, gender, tumor_size_group, etc)
kmfit_filtered <- survfit(surv_obj ~ clin$Treatment_groups)
ggsurvplot(kmfit_filtered, data = clin, pval = TRUE,risk.table = TRUE)

#Treatment groups epirubicin vs the rest
#subset_group <- clin[clin$Treatment_groups=="epirubicin" |clin$Treatment_groups=="NA/none",]
#surv_obj_sub <- Surv(time = subset_group$Follow.up.till.death, 
                     #event = subset_group$`Death.(1/.0)`=="1")
#kmfit2_filtered <- survfit(surv_obj_sub ~ subset_group$Treatment_groups)
#ggsurvplot(kmfit2_filtered, data = subset_group, pval = TRUE,risk.table = TRUE)

#Hazard ratio
fit.coxph<- coxph(surv_obj ~(Gender+Stage_baseline_MN+Stage_followup_MN+cytology_group), data = clin)
ggforest(fit.coxph, data = clin)

# glmnet
#Retrieve coefficients from the dataset(without BCG)
fit_glm <- glmnet(peptides,surv_obj,family="cox")# standardize = TRUE, maxit = 1000)
print(fit_glm)

cfs = coef(fit_glm,s=0.087130) 
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
meaning_vals
meaning_coefs
coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)

sorted_coef_abs <- coef_data[order(-abs(coef_data$coefficient)), ]

length(meaning_vals)
top <- head(sorted_coef_abs,length(meaning_vals))
list(top$variable)

# Assuming top is a data frame containing the top coefficients
paste(top$variable, collapse = ",")


# Create bar plot
ggplot(coef_data, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Genes", y = "Coefficient") +
  theme_minimal()


##DPD hazard ratio 
rownames(peptides) <- peptides$`CE-MS.Analysis.ID`
peptides_subset <- peptides[,coef_data$variable]
length(peptides_subset)
#dot product
row.names(peptides_subset)
colnames(coef_data)

coef_data_subset <-coef_data[,-1]
coef_data_subset
result <- as.matrix(peptides_subset) %*% as.matrix(coef_data_subset)

#assign result to DPD_prognosis
clin$DPD_prognosis <- result[row.names(result) %in% clin$`CE-MS.Analysis.ID`, ]

clin$outcome <- ifelse(clin$DPD_prognosis > 0, "positive",
                       ifelse(clin$DPD_prognosis == 0, "0", "negative"))

#KM curve ~ Variable of interest(DPD outcome)
fit <- survfit(Surv(time = clin$Follow.up.till.death, 
                    event = clin$`Death.(1/.0)`=="1") ~ outcome, data = clin)
ggsurvplot(fit, data = clin, pval = TRUE,xlab = "Time", ylab = "Overall survival",risk.table = TRUE)

#HR calculation
fit.coxph<- coxph(surv_obj ~ Gender+Stage_baseline_MN+Stage_followup_MN+cytology_group+outcome, data = clin,iter.max=30)
ggforest(fit.coxph, data = clin)
