library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("/Users/lidiayung/project/resource/perturbations/skcm_tcga_pan_can_atlas_2018")

clin_raw <- read.delim("data_clinical_patient.txt", sep = '\t',skip = 4)
rownames(clin_raw) <- clin_raw$PATIENT_ID

RNA_raw <- read.delim("data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",check.names = FALSE)
RNA_raw[is.na(RNA_raw)] <- 0
RNA_raw <- RNA_raw[RNA_raw$Hugo_Symbol!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Hugo_Symbol),]
rownames(RNA_raw) <- RNA_raw$Hugo_Symbol
RNA <- as.data.frame(t(RNA_raw[-1:-2]))
# Align clinical data:
clin <- clin_raw[str_sub(row.names(RNA), end = -4),]


# Return count of AJCC_PATHOLOGIC_TUMOR_STAGE per subgroup
stage_counts <- table(clin$AJCC_PATHOLOGIC_TUMOR_STAGE)
print(stage_counts)



###
surv_obj <- Surv(time = clin$OS_MONTHS, 
                 event = clin$OS_STATUS=="1:DECEASED")

survfit(formula = surv_obj ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clin)

clin_filtered <- clin[!is.na(clin$AJCC_PATHOLOGIC_TUMOR_STAGE), ]
survfit(formula = surv_obj ~ AJCC_PATHOLOGIC_TUMOR_STAGE, data = clin_filtered)

                        
#RNA selection

###
#1-year-survival function
calculate_survival_probability <- function(fit, time_point = 12, alpha = 0.05) {
  # Extract the 1-year survival probability and its standard error
  surv_prob_1_year <- summary(fit, times = time_point)$surv
  std_err_1_year <- summary(fit, times = time_point)$std.err
  
  # Calculate the z-score for the specified alpha level
  z_score <- qnorm(1 - alpha/2)
  
  # Calculate the lower and upper bounds of the confidence interval
  lower_bound <- surv_prob_1_year - z_score * std_err_1_year
  upper_bound <- surv_prob_1_year + z_score * std_err_1_year
  
  # Return the results as a named list
  result <- list(
    "1-year Survival Probability" = round(surv_prob_1_year, 3),
    "95% CI for 1-year Survival Probability" = c(round(lower_bound, 3), round(upper_bound, 3))
  )
  return(result)
}

#Filter Stage III and IV

clin_filtered <- clin[clin$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE III",  "STAGE IIIA",
                                                              "STAGE IIIB","STAGE IIIC",
                                                              "STAGE IV"), ]

#RNA selection
indices <- row.names(clin_filtered)

# Filter RNA based on the indices obtained
RNA_skmc <- RNA[indices,]



# Check if any row has all NA values
any_row_all_na <- apply(RNA_skmc, 1, function(row) all(is.na(row)))

# Check if there is at least one row with all NA values
if (any(any_row_all_na)) {
  print("There is at least one row with all NA values.")
} else {
  print("There are no rows with all NA values.")
}

# Assuming RNA_skmc is your data frame
# Check if any row has all NA values
rows_all_na <- apply(is.na(surv_filt), 1, all)

# Get the indices of rows where all values are NA
indices_all_na <- which(rows_all_na)
indices_all_na
print(paste("Indices of rows with all NA values:", indices_all_na))

RNA_skmc_nonna<- RNA_skmc[-indices_all_na,]
clin_filtered<-clin_filtered[-indices_all_na,]
# Create the survival object with Stage III and IV
my.surv_filtered <- Surv(clin_filtered$OS_MONTHS, clin_filtered$OS_STATUS == '1:DECEASED',)


# Fit the survival model for gender and plot
kmfit2_filtered <- survfit(my.surv_filtered ~ clin_filtered$SEX)
fit<- survfit(Surv(clin_filtered$OS_MONTHS, clin_filtered$OS_STATUS == '1:DECEASED') ~ SEX, data = clin_filtered)
ggsurvplot(fit, data = clin_filtered,
           risk.table = TRUE,                  # Add No at risk table
           cumevents = TRUE,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE    ,pval = TRUE)

#Overall
fit_skmc <- survfit(my.surv_filtered ~ 1, data = clin_filtered)

ggsurvplot(fit_skmc, data = clin, xlab = "Month", ylab = "Overall survival")
fit_skmc


###BLCA genes of interest
#TERT+FGFR3+KRAS+HRAS+PIK3CA+KDM6A+TSC1+CDKN2A+TP53+ERCC2+
#+ATM +ATR +BRCA1 +BRCA2 +POLE +FANCA+STAG2+PTEN+FOXA1


fit.coxph <- coxph(my.surv_filtered ~ EGFR +ATM +ATR +BRCA1 +BRCA2 , 
                   data = RNA_skmc)
ggforest(fit.coxph, data = RNA)
###

zscores.isna()

library(progeny)
zscores = as.matrix(t(RNA_skmc_nonna))
pathways <- progeny(zscores, scale=TRUE, organism="Human")  #, top = 100, perm = 1)
path_df = as.data.frame(pathways)
# fit multivariate model (COX proportional hazard) 
fit.coxph <- coxph(my.surv_filtered ~ MAPK + PI3K  +  p53 + TGFb  + `JAK-STAT` + TNFa, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)

fit.coxph <- coxph(my.surv_filtered ~ NFkB + Androgen +  VEGF + `JAK-STAT` + TGFb, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)

fit.coxph <- coxph(my.surv_filtered ~ MAPK + PI3K  +  p53 + TGFb  + `JAK-STAT` + TNFa + 
                     WNT + Androgen +  Estrogen + Hypoxia +  Trail + VEGF, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)



library(corrplot)
corrplot(cor(path_df), type = "upper", order = "hclust",tl.col = "black", tl.srt = 45)


# analysing on what MAPK activity depend
MAPK_df = cbind(MAPK = path_df$MAPK, PI3K = path_df$PI3K,
                KSR1 = RNA$KSR1, KSR2 = RNA$KSR2, IQGAP1 = RNA$IQGAP1, 
                IQGAP2 = RNA$IQGAP2, IQGAP3 = RNA$IQGAP3, GAB1 = RNA$GAB1, GAB2 = RNA$GAB2,
                KRAS = RNA$KRAS, NRAS = RNA$NRAS, HRAS = RNA$HRAS, BRAF=RNA$BRAF, 
                CRAF=RNA$RAF1, ARAF=RNA$ARAF)
corrplot(cor(MAPK_df), type = "upper", order = "hclust",tl.col = "black", tl.srt = 45)



# now creating object without zero times
clin_filtered <- clin_filtered[clin_filtered$OS_MONTHS > 0,]

RNA_filt <- RNA_skmc[clin_filtered$OS_MONTHS > 0,]
path_filt <- path_df[clin_filtered$OS_MONTHS > 0,]
# create a survival object consisting of times & censoring
surv_filt <- Surv(time = clin_filtered$OS_MONTHS, 
                  event = clin_filtered$OS_STATUS=="1:DECEASED")
fit <- survfit(surv_filt ~ 1, data = clin_filtered)
ggsurvplot(fit, data = clin_filtered, xlab = "Month", ylab = "Overall survival")


# glmnet
library("glmpath")
library("glmnet")
library("penalized")

fit_glm <- glmnet(RNA_filt,surv_filt,family="cox") #, alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)

# analysing results
cfs = coef(fit_glm,s=0.18)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
write.csv(meaning_vals,file='survival_coeffs_glmnet.csv')
print(paste(meaning_coefs,collapse=" + "))

fit_glm <- glmnet(RNA_skmc_nonna,surv_filt,family="cox") #, alpha = 1, standardize = TRUE, maxit = 1000

cfs = coef(fit_glm,s=0.18)

meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_coefs
meaning_vals = cfs[cfs[,1]!=0,]
write.csv(meaning_vals,file='survival_coeffs_glmnet.csv')
print(paste(meaning_coefs,collapse=" + "))

coef_data <-data.frame(gene=meaning_coefs,coefficient=meaning_vals)
# Create bar plot
ggplot(coef_data, aes(x = reorder(gene, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "gene", y = "Coefficient") +
  theme_minimal()



# cutting to find most important
ncut = 5
vals_surv = sort(abs(meaning_vals),decreasing = TRUE)[1:ncut]
print(paste(names(vals_surv),collapse=" + "))
# if we want to run coxph we have to copy-paste the string

fit.coxph <- coxph(surv_filt ~ SH2D2A, data = RNA_filt)
ggforest(fit.coxph, data = RNA_filt)


# glmnet over progeny
fit_glm <- glmnet(path_filt[-171,],surv_filt,family="cox") # , alpha = 1, standardize = TRUE, maxit = 1000
print(fit_glm)


# analysing results
cfs = coef(fit_glm,s=0.02)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
write.csv(meaning_vals,file='survival_coeffs_progeny_glmnet.csv')
print(paste(meaning_coefs,collapse=" + "))


coef_data <-data.frame(variable=meaning_coefs,coefficient=meaning_vals)
# Create bar plot
ggplot(coef_data, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Pathways", y = "Coefficient") +
  theme_minimal()

# cutting to find most important
ncut = 5
vals_surv = sort(abs(meaning_vals),decreasing = TRUE)[1:ncut]
print(paste(names(vals_surv),collapse=" + "))
# if we want to run coxph we have to copy-paste the string
fit.coxph <- coxph(surv_filt ~  NFkB + VEGF + `JAK-STAT` + Androgen + TGFb, data = path_filt)
ggforest(fit.coxph, data = RNA_filt)
colnames(path_filt)

# plotting Kaplan-Mayer curves
pathway = 'MAPK'
pathway_data = path_filt$NFkB
# sort age 
uni_path = sort(unique(pathway_data))
# store results
results_path = matrix(1,length(uni_path))
# do a for loop for every unique value of age mat
for (i in 2:(length(uni_path)-1)){ # Starting from 2 because the first element would yield no elements higher than.
  path_i = 1*(pathway_data>uni_path[i])
  # survdiff is the function from the survival package 
  logrank = survdiff(surv_filt ~ path_i)
  # store in results_age
  results_path[i] = logrank$pvalue
}
# Plot unique elements of age against p-value
plot(uni_path, results_path, log = 'y')
# Select minimum P-value
min_p_path = which.min(results_path)
# here are 2 good thresholds, -1 and 1
opt_thr = uni_path[min_p_path]
#opt_JAK = opt_thr
#opt_thr = 0.0
# I recalculated the P-value as a sanity check
pval = survdiff(surv_filt ~ pathway_data>opt_thr)$pvalue
nplus = sum(pathway_data>opt_thr)   # how many patients we have in high group
nminus = sum(pathway_data<opt_thr)   # how many patients we have in low group
# fit Kaplan Meier model
path_filt <- path_filt %>% mutate(MAPK_group = ifelse(MAPK >= 1, "high", ifelse(MAPK < 1.001*opt_thr,"low","intermediate")))
nhigh = sum(pathway_data>1)
ninter = sum((pathway_data<1) & (pathway_data > 1.001*opt_thr))
nlow = sum(pathway_data<1.001*opt_thr)
#KM = survfit(surv_filt ~ pathway_data>opt_thr,data = path_filt)
KM = survfit(surv_filt ~ MAPK_group,data = path_filt)
# Plot Kaplan Meier 
#plot(KM, lwd = 3, col = c(1,2), cex.axis = 1.5, xlab = 'Months', ylab = 'Survival Probability' , cex.lab = 1.5)
#ggsurvplot(KM, data = path_filt,pval = TRUE,xlab = 'Overall survival time, months',
#           legend.labs=c(paste('Low ',pathway,' activity, ',nminus,' patient(s)',sep = ''),paste('High ',pathway,' activity, ',nplus,' patient(s)',sep = '')),
#           palette = c('blue','red'),legend.title="")
p <- ggsurvplot(KM, data = path_filt,pval = TRUE,xlab = 'Overall survival time, months',
                legend.labs=c(paste("High MAPK activity,\n",nhigh," patients",sep=""),paste("Intermediate MAPK activity,\n",ninter," patients",sep=""),paste("Low MAPK activity,\n",nlow," patients",sep="")),
                legend.title=""
)

ggpar(p, 
      font.main = c(13, "bold"),
      font.x = c(16, "bold"),
      font.y = c(16, "bold"),
      font.caption = c(16, "bold"), 
      font.legend = c(13, "bold"), 
      font.tickslab = c(14, "bold"))


# RF
library("randomForestSRC")
# building RF model
B <- 1000
# Building a RSF
status_surv <- clin_filt$OS_STATUS=="1:DECEASED"
time_surv <- clin_filt$OS_MONTHS
dataSetRF <- cbind(time_surv,status_surv,RNA_filt)

names(dataSetRF)<-make.names(names(dataSetRF))

RF_obj <- rfsrc(Surv(time_surv,status_surv)~., dataSetRF,  ntree = B,  membership = TRUE, importance=TRUE)
# Printing the RF object  
print(RF_obj)

# Vadiable importance
jk.obj <- subsample(RF_obj)
#pdf("VIMPsur.pdf", width = 15, height = 20)
#par(oma = c(0.5, 10, 0.5, 0.5))
#par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
#plot(jk.obj, xlab = "Variable Importance (x 100)", cex = 1.2)
plot(jk.obj)
#dev.off()

fit.coxph <- coxph(surv_filt ~ TRIM67 + RASA1 + RHOF + NCAPD3 + NEU2 + ARNTL2 + CDK6, data = RNA_filt)
ggforest(fit.coxph, data = RNA_filt)