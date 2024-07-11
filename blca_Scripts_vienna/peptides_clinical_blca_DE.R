library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)
library(openxlsx)


#### Read files ####
clin_raw <- read.xlsx("Fwd_ MULTIR_ BC peptidomics datasets/BC_CE-MS peptidomics_Metadata.xlsx")
peptides <-read.xlsx("Fwd_ MULTIR_ BC peptidomics datasets/TransBioBC_08.04.2024_MosaID_1_7_5_MFinder_vs_MV_HybridSolution_v5_ML1_seq_md.xlsx",startRow=2)

clin <- clin_raw[clin_raw$`CE-MS.Analysis.ID` %in% peptides$`CE-MS.Analysis.ID`, ]



#### Understanding and refining groups ####
table(clin$Stage_grade_followup)
#Two groups of stages(MIBC:T2+= or NMIBC:T1)
clin$Stage_baseline_MN <- ifelse(grepl("^pT2",clin$Stage_grade_baseline),"MIBC","NMIBC")

#Two groups of stages(MIBC:T2+= or NMIBC:T1)
clin$Stage_followup_MN <- ifelse(grepl("^pT2|^pT3|^pT4",clin$Stage_grade_followup),"MIBC","NMIBC")

#tables
table(clin$Stage_grade_baseline)
table(clin$Stage_followup_MN)

#Treatment groups, BCG, epirubicin, MMC
clin$Treatment_groups <- ifelse(grepl("MMC", clin$Treatment) & grepl("epirubicin", clin$Treatment), "e-M Combined",
                                ifelse(grepl("BCG", clin$Treatment), "BCG",
                                       ifelse(grepl("MMC", clin$Treatment), "MMC",
                                              ifelse(grepl("epirubicin", clin$Treatment), "epirubicin", "NA/none"))))
table(clin$Treatment_groups)                                            

subset_group <- clin[clin$Treatment_groups=="epirubicin" |clin$Treatment_groups=="NA/none",]
surv_obj_sub <- Surv(time = subset_group$Follow.up.till.death, 
                 event = subset_group$`Death.(1/.0)`=="1")
kmfit2_filtered <- survfit(surv_obj_sub ~ subset_group$Treatment_groups)
ggsurvplot(kmfit2_filtered, data = subset_group, pval = TRUE,risk.table = TRUE)

#localisation
clin$localisation_groups <- ifelse(grepl("\\bwall\\b", clin$localisation), "bladder walls",
                                           ifelse(grepl("0", clin$localisation), "0",
                                                  ifelse(grepl("unknown", clin$localisation), "unknown", "bladder neck, ostium and urethra")))

table(clin$localisation)

#cytology groups
unique(clin$Cytology)

table(clin$Cytology)

# Assign cytology_group based on the conditions
clin$cytology_group <- case_when(
  grepl("high", clin$Cytology, ignore.case = TRUE) ~ "high_(grade 3)",
  grepl("low", clin$Cytology, ignore.case = TRUE) ~ "low_(grade 2/1)",
  grepl("suspect \\(grade 1\\)", clin$Cytology, ignore.case = TRUE) ~ "low (grade 2/1)",
  grepl("no abnormalities", clin$Cytology, ignore.case = TRUE) ~ "0/no cytology/no abnormalities/unknown",
  grepl("no cytology", clin$Cytology, ignore.case = TRUE) ~ "0/no cytology/no abnormalities/unknown",
  grepl("unknown", clin$Cytology, ignore.case = TRUE) ~ "0/no cytology/no abnormalities/unknown",
  TRUE ~ "0/no cytology/no abnormalities/unknown"
)

# Check the result
table(clin$cytology_group)

table(clin$`Death.(1/.0)`)
table(clin$Treatment)
#Filtered Treatment#Filtered by stages of interest if needed

#### KM curves and HZ ration ####
# create a survival object consisting of times & censoring
surv_obj <- Surv(time = clin$Follow.up.till.death, 
                 event = clin$`Death.(1/.0)`=="1")

#surv_obj 
fit <- survfit(surv_obj ~ 1, data = clin)
ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival",risk.table = TRUE)


# Fit the survival model for gender and plot
kmfit2_filtered <- survfit(surv_obj ~ clin$Gender)
ggsurvplot(kmfit2_filtered, data = clin, pval = TRUE,risk.table = TRUE)

# fit multivariate model (COX proportional hazard) 
fit.coxph<- coxph(surv_obj ~(Gender+Stage_baseline_MN+Stage_followup_MN), data = clin)
ggforest(fit.coxph, data = clin)


# now creating object without zero times (here all >0)
clin_filt <- clin[clin$Follow.up.till.death > 0,]


# create a survival object consisting of times & censoring
surv_filt <- Surv(time = clin_filt$Follow.up.till.death, 
                 event = clin$`Death.(1/.0)`=="1")

fit <- survfit(surv_filt ~ 1, data = clin_filt)
ggsurvplot(fit, data = clin_filt, xlab = "Month", ylab = "Overall survival")



#! surv_filt might contain null values because clin_flit contains
which(is.na(surv_filt))
#Remove rows that have null values

#### GLM, coefficients retrieval and dot product ####
# glmnet
library("glmpath")
library("glmnet")
library("penalized")
length(surv_filt)

fit_glm <- glmnet(peptides,surv_filt,family="cox")# standardize = TRUE, maxit = 1000)
print(fit_glm)

#glm_data <- data.frame(lambda = lambda_values, coef_glm)
# Write the data frame to a CSV file
#write.csv(glm_data, "glm_peptides_lambda.csv", row.names = FALSE)

# analysing results
cfs = coef(fit_glm,s=0.085230) #works 0.085230 #worked 0.016730
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
meaning_vals
meaning_coefs
coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)

#write.csv(coef_data,"coefs_peptides.csv")

sorted_coef_abs <- coef_data[order(-abs(coef_data$coefficient)), ]

top_5 <- head(sorted_coef_abs,5)
list(top_5$variable)

# Assuming top_100 is a data frame containing the top 100 coefficients
paste(top_5$variable, collapse = ",")


# Create bar plot
ggplot(coef_data, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Genes", y = "Coefficient") +
  theme_minimal()

print(paste(meaning_coefs,collapse=" + "))
# cutting to find most important
ncut = 5
vals_surv = sort(abs(meaning_vals),decreasing = TRUE)[1:ncut]
print(paste(names(vals_surv),collapse=" + "))
# if we want to run coxph we have to copy-paste the string
fit.coxph <- coxph(surv_filt ~ x99900895 + x99910097 + x99904775 + x99905720 + x99902448,data = peptides)
ggforest(fit.coxph, data = peptides)

#DPD hazard ratio 
peptides_subset <- peptides[,coef_data$variable]
peptides_subset
#dot product

colnames(coef_data)

coef_data_subset <-coef_data[,-1]
result <- as.matrix(peptides_subset) %*% as.matrix(coef_data_subset)
hist(result)
clin$DPD_prognosis <- result[str_sub(row.names(result)) %in% row.names(clin), ]

#write.csv(clin$DPD_prognosis, file = "DPD_prognosis.csv", row.names = FALSE)

clin$outcome <- ifelse(clin$DPD_prognosis > 0, "positive",
                       ifelse(clin$DPD_prognosis == 0, "0", "negative"))




#### KM curve and HZ ratio combining results above #####

# create a survival object consisting of times & censoring
surv_filt <- Surv(time = clin$Follow.up.till.death, 
                  event = clin$`Death.(1/.0)`=="1")

fit <- survfit(surv_obj ~ 1, data = clin)
ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival",risk.table = TRUE)


# Fit the survival model for gender and plot
kmfit2_filtered <- survfit(surv_obj ~ clin$outcome)
ggsurvplot(kmfit2_filtered, data = clin, pval = TRUE,risk.table = TRUE)
kmfit2_filtered

fit <- survfit(Surv(time = clin$Follow.up.till.death, 
                    event = clin$`Death.(1/.0)`=="1") ~ outcome, data = clin)

ggsurvplot(fit, data = clin, pval = TRUE,xlab = "Month", ylab = "Overall survival",risk.table = TRUE)


#Complete separation leads to problems from here on
fit.coxph<- coxph(surv_obj ~ Gender+Stage_baseline_MN+Stage_followup_MN+cytology_group, data = clin,iter.max=30)
ggforest(fit.coxph, data = clin)

summary(fit.coxph)

