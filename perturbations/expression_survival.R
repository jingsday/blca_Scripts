library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)




setwd("/Users/lidiayung/project/resource/perturbations/lusc_tcga_pan_can_atlas_2018")

clin_raw <- read.delim("data_clinical_patient.txt", sep = '\t',skip = 4)
#Filtered by stages of interest
table(clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE)

clin_raw <- clin_raw[clin_raw$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE III","STAGE IIIA", "STAGE IIIB",
                                                                 "STAGE IIIC","STAGE IV"), ]
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



any(colnames(clin) %in% c("OS_STATUS", "DAYS_LAST_FOLLOWUP", "OS_MONTHS"))
which(colnames(clin) %in% c("OS_STATUS", "DAYS_LAST_FOLLOWUP", "OS_MONTHS"))


#RNA<- RNA[row.names(RNA) %in% clin$col_rna,]

# days_to_death, that is the number of days passed from the initial diagnosis to the patient’s death (clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.


# change certain values the way they are encoded
clin$deceased <- ifelse(clin$OS_STATUS == "0:LIVING", FALSE, TRUE)


# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clin$overall_survival <- ifelse(is.na(clin$OS_STATUS == "0:LIVING"),
                                clin$DAYS_LAST_FOLLOWUP,
                                clin$OS_MONTHS)
clin$overall_survival[clin$OS_MONTHS == 0] <- NA

table(clin$OS_STATUS)

### Categorize into dysregulated and intact groups
dys_rna <- t(apply(RNA, 1, function(x) ifelse(abs(x) > 1.96,"dysregulated","intact")))




###

fin_dat <- data.frame(gene = dys_rna[, colnames(dys_rna) == "C14orf4"])

fin_dat$PATIENT_ID<- str_sub(row.names(fin_dat), end = -4)


# Merge the data frames based on the ID column
fin_dat <- merge(fin_dat, clin, by = "PATIENT_ID", all.x = TRUE)


# fitting model
fit1 <- survfit(Surv(overall_survival, deceased) ~ gene, data = fin_dat)
print(fit1)


# calculating pvalue
fit2 <- survdiff(Surv(overall_survival, deceased) ~ gene, data = fin_dat)
pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
print(pv)


# KM plot,Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit1,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(),
           palette = c("#990000", "#000099"))



survdiff.fit <- function(y, x, strat, rho=0) {
  #
  # This routine is almost always called from survdiff
  #  If called directly, remember that it does no error checking
  #
  n <- length(x)
  if (ncol(y) !=2) stop ("Invalid y matrix")
  if (nrow(y) !=n | length(x) !=n) stop("Data length mismatch")
  
  ngroup <- length(unique(x))
  if (ngroup <2) stop ("There is only 1 group")
  if (inherits(x, "factor")) x <- as.numeric(x)
  else x <- match(x, unique(x))
  
  if (missing(strat)) strat <- rep(1,n)
  else strat <- as.numeric(as.factor(strat))
  nstrat <- length(unique(strat))
  if (length(strat) !=n) stop("Data length mismatch")
  
  ord <- order(strat, y[,1], -y[,2])
  strat2 <- c(1*(diff(strat[ord])!=0), 1)
  
  xx <- .C(Csurvdiff2, as.integer(n),
           as.integer(ngroup),
           as.integer(nstrat),
           as.double(rho),
           as.double(y[ord,1]),
           as.integer(y[ord,2]),
           as.integer(x[ord]),
           as.integer(strat2),
           observed = double(ngroup*nstrat),
           expected = double(ngroup*nstrat),
           var.e    = double(ngroup * ngroup),
           double(ngroup), double(n))
  
  if (nstrat==1)  list(expected = xx$expected,
                       observed = xx$observed,
                       var      = matrix(xx$var.e, ngroup, ngroup))
  else            list(expected = matrix(xx$expected, ngroup),
                       observed = matrix(xx$observed, ngroup),
                       var      = matrix(xx$var.e, ngroup, ngroup))
}



#Performing survival analysis for all genes
all_gene <- colnames(dys_rna)
result = data.frame( gene=character(0), pval=numeric(0), dysregulated=numeric(0), intact=numeric(0))

for (i in all_gene){
  fin_dat <- data.frame(gene = dys_rna[, colnames(dys_rna) == i])
  fin_dat$PATIENT_ID<- str_sub(row.names(fin_dat), end = -4)
  
  fin_dat <- merge( fin_dat, clin, by = "PATIENT_ID",all.x = TRUE)
  for (i in all_gene) {
    fin_dat <- data.frame(gene = dys_rna[, colnames(dys_rna) == i])
    fin_dat$PATIENT_ID <- str_sub(row.names(fin_dat), end = -4)
    
    fin_dat <- merge(fin_dat, clin, by = "PATIENT_ID", all.x = TRUE)
    
    if (dim(table(fin_dat$gene)) > 1 | any(colSums(fin_dat > 2) > 0)) {
      # Use tryCatch to handle the error and continue
      fit2 <- tryCatch({
        survdiff(Surv(overall_survival, deceased) ~ gene, data = fin_dat)
      }, error = function(e) {
        return(NULL)  # Return NULL if an error occurs
      })
      
      if (!is.null(fit2)) {
        pv <- ifelse(is.na(fit2), next, (round(1 - pchisq(fit2$chisq, length(fit2$n) - 1), 3)))[[1]]
        
        gene <- i
        dysregulated <- table(fin_dat$gene)[1]
        intact <- table(fin_dat$gene)[2]
        pval <- pv
        result[i, ] <- c(gene, pval, dysregulated, intact)
      }
    }
  }
  


# retrieve info
gene_names <- c("EMP1", "FGFR1", "TPM1", "NRP2", "LATS2")
genes_subset <- result[result$gene %in% gene_names, ]

result_sig <- result[result$pval<0.5, ]

###


################ Finding genes which has values under the three categories High, Norm, Low
#mapping high, low, medium
#low (z-score <= -1.96), normal (-1.96 < z-score > 1.96), and high (z-score >= 1.96)
dys_rna <- t(apply(RNA, 1, function(x) {
  ifelse(x <= -1.96, "Low", 
         ifelse(x >= 1.96, "High", "Normal"))
}))

gene_table <- data.frame( gene=character(0), High=numeric(0), Norm = numeric(0), Low=numeric(0), lab = character(0))
gene <- colnames(dys_rna)
for (i in gene) {
  df <- data.frame(gene = dys_rna[,colnames(dys_rna) == i ])
  tb <- table(df)
  if (dim(tb) == 3){
    gene <- i
    High <- tb[1]
    Low <- tb[2]
    Norm <- tb[3]
    lab <- paste(names(tb)[1],names(tb)[2],names(tb)[3], sep = "_" )
    gene_table[i, ] = c(gene, High, Norm, Low, lab)
  }
  
}

############## finding genes which has data only in two of the three states (High,Norm,Low)
gene_table_2x2 <- data.frame( gene=character(0), state1=numeric(0), state2= numeric(0), lab1 = character(0), lab2 = character(0))

gene <- colnames(dys_rna)
for (i in gene) {
  df <- data.frame(gene = dys_rna[,colnames(dys_rna) == i ])
  tb <- table(df)
  if (dim(tb) == 2){
    gene <- i
    state1 <- tb[1]
    state2 <- tb[2]
    lab1 <- paste(names(tb)[1])
    lab2 <- paste(names(tb)[2])
    gene_table_2x2[i, ] = c(gene, state1, state2, lab1, lab2)
  }
  
}
#h.table
h.tab <- gene_table_2x2[gene_table_2x2$lab1 == "High", ]
names(h.tab)[2] <- names(table(h.tab$lab1))[1]
names(h.tab)[3] <- "Normal"
h.tab$Low <- NA
#h.tab <- h.tab[, colnames(gene_table)[1:4]]

#l.table
l.tab <- gene_table_2x2[gene_table_2x2$lab1 == "Low", ]
names(l.tab)[2] <- names(table(l.tab$lab1))[1]
names(l.tab)[3] <- names(table(l.tab$lab2))[1]
l.tab$High <- NA
#l.tab <- l.tab[, colnames(gene_table)[1:4]]
#
gene_table_2x2 <- rbind(h.tab, l.tab)
gene_table_2x2$lab <- NA
#

###Make sure columns counts the same
to_amend <- gene_table_2x2[c("gene", "High", "Normal", "Low", "lab")]
names(to_amend)[3] <- "Norm"

gene_table <- rbind(gene_table, to_amend)



#gene_table <- gene_table[,-5]

# geting the result table from the previous analysis:
gene_table[,2:4] <- sapply(gene_table[, 2:4], as.numeric) # setting type of columns for numbers as numeric
gene_table[is.na(gene_table)] <- 0
###


# defining classes for each gene
gene_table$state <- ifelse(gene_table$Norm < 15 & gene_table$High >= 15 & gene_table$Low >= 15, "HL", 
                           ifelse(gene_table$Low < 15 & gene_table$High >= 15 & gene_table$Norm >= 15, "HN", 
                                  ifelse(gene_table$High < 15 & gene_table$Low >= 15 & gene_table$Norm >= 15, "LN",
                                         ifelse(gene_table$High >= 15 & gene_table$Low >= 15 & gene_table$Norm >= 15, "HNL", "flag"))))
# see what we get
table(gene_table$state)
##################Here
#hnl
length(gene_table[gene_table$state == "HNL", ]$gene)

hnl.dys_rna <- t(apply(RNA[,colnames(RNA) %in% gene_table[gene_table$state == "HNL", ]$gene ], 1, function(x) ifelse(x >= 1.96,"High",ifelse(x <= -1.96, "Low", "Norm"))))

gene_t <- colnames(hnl.dys_rna)

hnl.result = data.frame( gene=character(0), pval=numeric(0), High=numeric(0), Norm=numeric(0), Low=numeric(0))

for (i in gene){
  fin_dat <- data.frame(gene = dys_rna[,colnames(dys_rna) == i ])
  fin_dat$PATIENT_ID<-clin$PATIENT_ID
  
  fin_dat <- merge( fin_dat, clin, by = "PATIENT_ID")
  if (dim(table(fin_dat$gene)) > 2){
    fit2 <- survdiff(Surv(overall_survival, deceased) ~ gene, data = fin_dat)
    
    pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
    
    gene <- i
    High <- table(fin_dat$gene)[1]
    Norm <- table(fin_dat$gene)[3]
    Low <- table(fin_dat$gene)[2]
    pval = pv
    hnl.result[i, ] = c(gene, pval, High,Norm, Low)
  }
}

#hn
hn.dys_rna <- t(apply(RNA[,colnames(RNA) %in% gene_table[gene_table$state == "HN", ]$gene ], 
                      1, function(x) ifelse(x >= 1.96,"High",ifelse(x <= -1.96, "Low", "Norm"))))
gene <- gene_table[gene_table$newcol == "HN", ]$gene
hn.result = data.frame( gene=character(0), pval=numeric(0), High=numeric(0), Norm=numeric(0))

for (i in all_gene){

  fin_dat <- data.frame(gene = hn.dys_rna[,colnames(hn.dys_rna) == i ])
  fin_dat$PATIENT_ID<-clin$PATIENT_ID
  
  fin_dat <- merge( fin_dat, clin, by = "PATIENT_ID")
  fin_dat <- fin_dat[-which(fin_dat$gene == "Low"), ]
  if (dim(table(fin_dat$gene)) > 1){

    fit2 <- survdiff(Surv(overall_survival, deceased) ~ gene, data = fin_dat)

    pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
    #
    gene <- i
    High <- table(fin_dat$gene)[1]
    Norm <- table(fin_dat$gene)[2]
    pval = pv
    hn.result[i, ] = c(gene, pval, High, Norm)
  }
}

#ln

ln.dys_rna <- t(apply(RNA[,colnames(RNA) %in% gene_table[gene_table$state == "LN", ]$gene ],1, function(x) ifelse(x >= 1.96,"High",ifelse(x <= -1.96, "Low", "Norm"))))

gene <- gene_table[gene_table$newcol == "LN", ]$gene
ln.result = data.frame( gene=character(0), pval=numeric(0), Low=numeric(0), Norm=numeric(0))

for (i in all_gene){
  fin_dat <- data.frame(gene = ln.dys_rna[,colnames(ln.dys_rna) == i ])
  fin_dat$PATIENT_ID<-clin$PATIENT_ID
  
  fin_dat <- merge( fin_dat, clin, by = "PATIENT_ID")
  fin_dat <- fin_dat[-which(fin_dat$gene == "High"), ]
  if (dim(table(fin_dat$gene)) > 1){
    fit2 <- survdiff(Surv(overall_survival, deceased) ~ gene, data = fin_dat)
    pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
    #
    gene <- i
    Low <- table(fin_dat$gene)[1]
    Norm <- table(fin_dat$gene)[2]
    pval = pv
    ln.result[i, ] = c(gene, pval, High, Norm)
  }
}

# combining results into one dataset
hn.result$Low <- NA
ln.result$High <- NA
hn.result <- hn.result[, colnames(hnl.result)]
ln.result <- ln.result[, colnames(hnl.result)]
#
result <- rbind(hnl.result, hn.result, ln.result)
# selecting significantly associated genes
result$pval <- as.numeric(result$pval)
sig.result <- result[result$pval <= 0.05 ,]



write.table(sig.result, file = "stages_sig_blca.surv.associated.gene.csv", row.names = T, quote = F)


#Plot centain genes
fin_dat <- data.frame(gene = dys_rna[, colnames(dys_rna) == "BDP1"])

fin_dat <- data.frame(gene = dys_rna[, colnames(dys_rna) == "ABCC1"], PATIENT_ID = row.names(RNA))
fin_dat$PATIENT_ID<- str_sub(fin_dat$PATIENT_ID, end = -4)

# Merge the data frames based on the ID column
fin_dat <- merge(fin_dat, clin, by = "PATIENT_ID", all.x = TRUE)


# KM plot,Change color, linetype by strata, risk.table color by strata


# fitting model
fit1 <- survfit(Surv(overall_survival, deceased) ~ gene, data = fin_dat)
print(fit1)


# calculating pvalue
fit2 <- survdiff(Surv(overall_survival, deceased) ~ gene, data = fin_dat)
pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
print(pv)

ggsurvplot(fit1,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(),)+
  ggtitle("Survival analysis with ABCC1 expression")

