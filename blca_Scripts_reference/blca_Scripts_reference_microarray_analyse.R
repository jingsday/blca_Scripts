# we have mouse4302cdf
# cdf=Mouse430_2 (45101 affyids)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("affy")
#BiocManager::install("gcrma")
#BiocManager::install("mouse4302.db")

library(affy)
library(gcrma)
library(limma)

library(annotate)
library(mouse4302cdf)
library(mouse4302probe)
library("mouse4302.db")

cels = list.files("all_data/", pattern = "cel")
setwd('all_data')
raw_data=ReadAffy(verbose=TRUE, filenames=cels, cdfname = "mouse4302") #From bioconductor

# normalizing the data
data.gcrma = gcrma(raw_data)
# alternative variant
#my.affinity.info = compute.affinities.local(raw_data)
#data.gcrma = gcrma(raw_data,affinity.info=my.affinity.info)

# just stubs
#ID = featureNames(raw_data)
#gNames = getSYMBOL(ID,"mouse4302.db")
#cdfName(raw_data)
#featureNames(raw_data)
#probeNames(raw_data)
#raw_data@experimentData
#raw_data@featureData

ph = raw_data@phenoData

# 1 variable classification
ph@data[ ,2] = c("B_I0_M0","B_I0_M0","B_I0_M0","B_I0_M30","B_I0_M30","B_I0_M30","B_I0_M6","B_I0_M6","B_I0_M6","B_I0_M78","B_I0_M78","B_I0_M78","B_I24_M0","B_I24_M0","B_I24_M0","B_I24_M30","B_I24_M30","B_I24_M30","B_I24_M6","B_I24_M6","B_I24_M6","B_I24_M78","B_I24_M78","B_I24_M78",
                 "B1H_I0_M0","B1H_I0_M0","B1H_I0_M0","B1H_I0_M30","B1H_I0_M30","B1H_I0_M30","B1H_I0_M6","B1H_I0_M6","B1H_I0_M6","B1H_I0_M78","B1H_I0_M78","B1H_I0_M78","B1H_I24_M0","B1H_I24_M0","B1H_I24_M0","B1H_I24_M30","B1H_I24_M30","B1H_I24_M30","B1H_I24_M6","B1H_I24_M6","B1H_I24_M6","B1H_I24_M78","B1H_I24_M78","B1H_I24_M78",
                 "H_I0_M0","H_I0_M0","H_I0_M0","H_I0_M30","H_I0_M30","H_I0_M30","H_I0_M6","H_I0_M6","H_I0_M6","H_I0_M78","H_I0_M78","H_I0_M78","H_I24_M0","H_I24_M0","H_I24_M0","H_I24_M30","H_I24_M30","H_I24_M30","H_I24_M6","H_I24_M6","H_I24_M6","H_I24_M78","H_I24_M78","H_I24_M78",
                 "H1B_I0_M0","H1B_I0_M0","H1B_I0_M0","H1B_I0_M30","H1B_I0_M30","H1B_I0_M30","H1B_I0_M6","H1B_I0_M6","H1B_I0_M6","H1B_I0_M78","H1B_I0_M78","H1B_I0_M78","H1B_I24_M0","H1B_I24_M0","H1B_I24_M0","H1B_I24_M30","H1B_I24_M30","H1B_I24_M30","H1B_I24_M6","H1B_I24_M6","H1B_I24_M6","H1B_I24_M78","H1B_I24_M78","H1B_I24_M78")
colnames(ph@data)[2]="Sample_name"
groups = ph@data$Sample_name
f = factor(groups,levels = c("B_I0_M0","B_I0_M30","B_I0_M6","B_I0_M78","B_I24_M0","B_I24_M30","B_I24_M6","B_I24_M78",
                             "B1H_I0_M0","B1H_I0_M30","B1H_I0_M6","B1H_I0_M78","B1H_I24_M0","B1H_I24_M30","B1H_I24_M6","B1H_I24_M78",
                             "H_I0_M0","H_I0_M30","H_I0_M6","H_I0_M78","H_I24_M0","H_I24_M30","H_I24_M6","H_I24_M78",
                             "H1B_I0_M0","H1B_I0_M30","H1B_I0_M6","H1B_I0_M78","H1B_I24_M0","H1B_I24_M30","H1B_I24_M6","H1B_I24_M78"))
design = model.matrix(~ 0 + f)
colnames(design) = levels(f)
data.fit = lmFit(data.gcrma, design)

contrast.matrix = makeContrasts(B_I0_M0-B_I0_M0,B_I0_M6-B_I0_M0,B_I0_M30-B_I0_M0,B_I0_M78-B_I0_M0,
                                B_I24_M0-B_I0_M0,B_I24_M6-B_I0_M0,B_I24_M30-B_I0_M0,B_I24_M78-B_I0_M0,
                                B1H_I0_M0-B_I0_M0,B1H_I0_M6-B_I0_M0,B1H_I0_M30-B_I0_M0,B1H_I0_M78-B_I0_M0,
                                B1H_I24_M0-B_I0_M0,B1H_I24_M6-B_I0_M0,B1H_I24_M30-B_I0_M0,B1H_I24_M78-B_I0_M0,
                                H_I0_M0-B_I0_M0,H_I0_M6-B_I0_M0,H_I0_M30-B_I0_M0,H_I0_M78-B_I0_M0,
                                H_I24_M0-B_I0_M0,H_I24_M6-B_I0_M0,H_I24_M30-B_I0_M0,H_I24_M78-B_I0_M0,
                                H1B_I0_M0-B_I0_M0,H1B_I0_M6-B_I0_M0,H1B_I0_M30-B_I0_M0,H1B_I0_M78-B_I0_M0,
                                H1B_I24_M0-B_I0_M0,H1B_I24_M6-B_I0_M0,H1B_I24_M30-B_I0_M0,H1B_I24_M78-B_I0_M0,
                                levels=design)
data.fit.con = contrasts.fit(data.fit,contrast.matrix)

# 2 variable classification
#ph@data[ ,2] = c("I0_M0","I0_M0","I0_M0","I0_M30","I0_M30","I0_M30","I0_M6","I0_M6","I0_M6","I0_M78","I0_M78","I0_M78","I24_M0","I24_M0","I24_M0","I24_M30","I24_M30","I24_M30","I24_M6","I24_M6","I24_M6","I24_M78","I24_M78","I24_M78",
#                 "I0_M0","I0_M0","I0_M0","I0_M30","I0_M30","I0_M30","I0_M6","I0_M6","I0_M6","I0_M78","I0_M78","I0_M78","I24_M0","I24_M0","I24_M0","I24_M30","I24_M30","I24_M30","I24_M6","I24_M6","I24_M6","I24_M78","I24_M78","I24_M78",
#                 "I0_M0","I0_M0","I0_M0","I0_M30","I0_M30","I0_M30","I0_M6","I0_M6","I0_M6","I0_M78","I0_M78","I0_M78","I24_M0","I24_M0","I24_M0","I24_M30","I24_M30","I24_M30","I24_M6","I24_M6","I24_M6","I24_M78","I24_M78","I24_M78",
#                 "I0_M0","I0_M0","I0_M0","I0_M30","I0_M30","I0_M30","I0_M6","I0_M6","I0_M6","I0_M78","I0_M78","I0_M78","I24_M0","I24_M0","I24_M0","I24_M30","I24_M30","I24_M30","I24_M6","I24_M6","I24_M6","I24_M78","I24_M78","I24_M78")
#colnames(ph@data)[2]="Treatment"
#ph@data[ ,3] = c("B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B","B",
#                 "B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H","B1H",
#                 "H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H",
#                 "H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B","H1B")
#colnames(ph@data)[3]="Cell_line"
#groupsC = ph@data$Cell_line 
#groupsT = ph@data$Treatment
#fc = factor(groupsC,levels=c("B","B1H","H","H1B"))
#ft = factor(groupsT,levels=c("I0_M0","I0_M30","I0_M6","I0_M78","I24_M0","I24_M30","I24_M6","I24_M78"))
#design = model.matrix(~ fc + ft)
#data.fit = lmFit(data.gcrma,design)


data.fit.eb = eBayes(data.fit.con)
#data.fit.eb$coefficients

setwd('..')

# saving results
res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=1)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B_I0_M0_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=2)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B_I0_M6_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=3)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B_I0_M30_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=4)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B_I0_M78_vs_B_I0_M0_DEGs.csv",row.names = FALSE)


res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=5)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B_I24_M0_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=6)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B_I24_M6_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=7)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B_I24_M30_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=8)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B_I24_M78_vs_B_I0_M0_DEGs.csv",row.names = FALSE)


res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=9)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B1H_I0_M0_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=10)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B1H_I0_M6_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=11)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B1H_I0_M30_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=12)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B1H_I0_M78_vs_B_I0_M0_DEGs.csv",row.names = FALSE)


res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=13)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B1H_I24_M0_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=14)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B1H_I24_M6_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=15)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B1H_I24_M30_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=16)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"B1H_I24_M78_vs_B_I0_M0_DEGs.csv",row.names = FALSE)


res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=17)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H_I0_M0_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=18)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H_I0_M6_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=19)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H_I0_M30_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=20)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H_I0_M78_vs_B_I0_M0_DEGs.csv",row.names = FALSE)


res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=21)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H_I24_M0_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=22)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H_I24_M6_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=23)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H_I24_M30_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=24)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H_I24_M78_vs_B_I0_M0_DEGs.csv",row.names = FALSE)


res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=25)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H1B_I0_M0_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=26)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H1B_I0_M6_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=27)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H1B_I0_M30_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=28)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H1B_I0_M78_vs_B_I0_M0_DEGs.csv",row.names = FALSE)


res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=29)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H1B_I24_M0_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=30)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H1B_I24_M6_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=31)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H1B_I24_M30_vs_B_I0_M0_DEGs.csv",row.names = FALSE)

res<-topTable(data.fit.eb, number=Inf, adjust.method="BH", coef=32)
gprobes=row.names(res)
Symbols = unlist(mget(gprobes, mouse4302SYMBOL, ifnotfound=NA))
EntrezIDs = unlist(mget(gprobes, mouse4302ENTREZID, ifnotfound=NA))
res.names=cbind(gprobes,Symbols,EntrezIDs,res)
write.csv(res.names,"H1B_I24_M78_vs_B_I0_M0_DEGs.csv",row.names = FALSE)


