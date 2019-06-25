# This script is for generating sig.score in TNBC/Basal samples
# Jiayi Tan
# 06/24/2019

rm(list = ls(all.names = TRUE))
cat("Loading genefu library...\n")
library(genefu)
cat("Loading METABRIC expression data...\n")
st = Sys.time()
meta_expr<-readRDS("metabric_expr_ilid.RDS")
print(Sys.time()-st)
data.meta<-t(meta_expr)

cat("Loading clinical data...\n")
library("xlsx")
st = Sys.time()
clin_info<-read.xlsx("merge_clin_info_manual_checked.xlsx",1,header = TRUE)
print(Sys.time()-st)

cat("Select TNBC/Basal samples from clinical data frame\n")
Basal_info<-subset(clin_info, pam50_claudin =="Basal")
TNBC_info<-subset(clin_info, er_status =="-"& her2_status =="-"& pr_status =="-")
TNBC_Basal<-rbind(Basal_info,TNBC_info)
TNBC_Basal<-unique(TNBC_Basal)

cat("Generate expression data from TNBC/Basal samples\n")
data.meta_TB<-subset(data.meta,rownames(data.meta) %in% TNBC_Basal$pid)

cat("Loading METABRIC annotation data...\n")
library(readxl)
annot.meta<-read_excel("HumanHT_12_v30_R3_cleaned.xlsx", sheet = 1)
annot.meta<-annot.meta[c(5,9,14)]
colnames(annot.meta)[1]="NCBI.gene.symbol"
colnames(annot.meta)[2]="EntrezGene.ID"
colnames(annot.meta)[3]="probe"

cat("Load TRM signature gene expression data\n")
TRM_sig<-read_excel("TRM_sig1.xlsx", sheet = 2)
TRM_sig<-TRM_sig[1:2]
colnames(TRM_sig)[1]="NCBI.gene.symbol"
colnames(TRM_sig)[2]="coefficient"
TRM_sig$EntrezGene.ID<-annot.meta$EntrezGene.ID[match(TRM_sig$NCBI.gene.symbol, annot.meta$NCBI.gene.symbol)]
TRM_sig$probe<-annot.meta$probe[match(TRM_sig$NCBI.gene.symbol, annot.meta$NCBI.gene.symbol)]

cat("Run TRM gene sig.score in TNBC/Basal samples\n")
TRM.score <- sig.score(x=TRM_sig, data=data.meta_TB, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)