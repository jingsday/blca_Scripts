library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(stringr)

setwd("/Users/lidiayung/project/resource/perturbations/blca_tcga_pan_can_atlas_2018")



clin_raw <- read.delim("data_clinical_patient.txt", sep = '\t',skip = 4)

table(clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE)
#Filtered by stages of interest

clin_raw <- clin_raw[clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE II","STAGE III", "STAGE IV"), ]
#double check
table(clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE)

rownames(clin_raw) <- clin_raw$PATIENT_ID
length(clin_raw$PATIENT_ID)

RNA_raw <- read.delim("data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",check.names = FALSE)
RNA_raw[is.na(RNA_raw)] <- 0
RNA_raw <- RNA_raw[RNA_raw$Hugo_Symbol!='',]
RNA_raw <- RNA_raw[!duplicated(RNA_raw$Hugo_Symbol),]
rownames(RNA_raw) <- RNA_raw$Hugo_Symbol
RNA <- as.data.frame(t(RNA_raw[-1:-2]))

#retrieve RNAs of interest
RNA <- RNA[str_sub(row.names(RNA), end = -4) %in% row.names(clin_raw), ]


clin <- clin_raw[str_sub(row.names(RNA), end = -4),]

# create a survival object consisting of times & censoring
surv_obj <- Surv(time = clin$OS_MONTHS, 
                 event = clin$OS_STATUS=="1:DECEASED")

#surv_obj 

fit <- survfit(surv_obj ~ 1, data = clin)
ggsurvplot(fit, data = clin, xlab = "Month", ylab = "Overall survival",surv.median.line = "hv")


# fit multivariate model (COX proportional hazard) 
fit.coxph <- coxph(surv_obj ~ EGFR + ERBB2 + ERBB3 + ERBB4 + NRG1 + NRG2 + NRG3 + NRG4, 
                   data = RNA)
ggforest(fit.coxph, data = RNA)
#new.coxph <- predict(fit.coxph, newdata = RNA[,c("KRAS", "ERBB2")])
#new.coxph
fit.coxph <- coxph(surv_obj ~ KSR1 + KSR2 + IQGAP1 + GAB1 + GAB2, 
                   data = RNA)
ggforest(fit.coxph, data = RNA)




library(progeny)
zscores = as.matrix(t(RNA))
pathways <- progeny(zscores, scale=TRUE, organism="Human")  #, top = 100, perm = 1)
path_df = as.data.frame(pathways)
# fit multivariate model (COX proportional hazard) 
fit.coxph <- coxph(surv_obj ~ MAPK + PI3K  +  p53 + TGFb  + `JAK-STAT` + TNFa, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)

fit.coxph <- coxph(surv_obj ~ NFkB + Androgen +  VEGF + `JAK-STAT` + TGFb, 
                   data = path_df)
ggforest(fit.coxph, data = path_df)

fit.coxph <- coxph(surv_obj ~ MAPK + PI3K  +  p53 + TGFb  + `JAK-STAT` + TNFa + 
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

cor(MAPK_df)



# now creating object without zero times
clin_filt <- clin[clin$OS_MONTHS > 0,]
RNA_filt <- RNA[clin$OS_MONTHS > 0,]
path_filt <- path_df[clin$OS_MONTHS > 0,]
# create a survival object consisting of times & censoring
surv_filt <- Surv(time = clin_filt$OS_MONTHS, 
                  event = clin_filt$OS_STATUS=="1:DECEASED")
fit <- survfit(surv_filt ~ 1, data = clin_filt)
ggsurvplot(fit, data = clin_filt, xlab = "Month", ylab = "Overall survival")

#removal NAs
which(is.na(surv_filt))
#Remove rows that have null values
surv_filt<- surv_filt[-261,]
RNA_filt<- RNA_filt[-261,]

# glmnet
library("glmpath")
library("glmnet")
library("penalized")
fit_glm <- glmnet(RNA_filt,surv_filt,family="cox") # , alpha = 1, standardize = TRUE, maxit = 1000
plot(fit_glm)
print(fit_glm)
# analysing results

cfs = coef(fit_glm,s=0.18)
meaning_coefs = rownames(cfs)[cfs[,1]!= 0]
meaning_vals = cfs[cfs[,1]!=0,]
write.csv(meaning_vals,file='survival_coeffs_glmnet.csv')
print(paste(meaning_coefs,collapse=" + "))
# cutting to find most important
ncut = 15
vals_surv = sort(abs(meaning_vals),decreasing = TRUE)[1:ncut]
print(paste(names(vals_surv),collapse=" + "))
# if we want to run coxph we have to copy-paste the string
fit.coxph <- coxph(surv_filt ~ ACTL6A + CASKIN2 + KIAA0195 + LY6D + MET + MRPL3 + PWWP3A + SEC61A2 + SFXN5 + THSD1P1 + TRIM67 + UCA1 + USP20, data = RNA_filt)
ggforest(fit.coxph, data = RNA_filt)


# CV glmnet
cvfit <- cv.glmnet(data.matrix(RNA_filt),surv_filt,family="cox",type.measure = "C")
plot(cvfit)


ggplot(coef_data, aes(x = reorder(variable, coefficient), y = coefficient)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotate axis labels
  labs(title = "Coefficients from glmnet Model", 
       x = "Genes", y = "Coefficient") +
  theme_minimal()
