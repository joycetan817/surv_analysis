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
st = Sys.time()
clin_info<-readRDS("merge_clin_info_manual_checked.RDS")
print(Sys.time()-st)

cat("Select TNBC/Basal samples from clinical data frame\n")
Basal_info<-subset(clin_info, Pam50Subtype =="Basal")
TNBC_info<-subset(clin_info, er_status =="-"& her2_status =="-"& pr_status =="-")
TNBC_Basal<-rbind(Basal_info,TNBC_info)
TNBC_Basal<-unique(TNBC_Basal)

cat("Generate expression data from TNBC/Basal samples\n")

data.meta_T_B<-subset(data.meta,rownames(data.meta) %in% TNBC_Basal$pid)


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
TRM.score <- sig.score(x=TRM_sig, data=data.meta_T_B, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)



TRM_score<- as.data.frame (TRM.score$score)
library(tibble)
TRM_score<-rownames_to_column(TRM_score,var ="patient_id")
colnames(TRM_score)[2]="sig.score"
TRM_score<-TRM_score[order(-TRM_score$sig.score),]
library(dplyr)
top_25<-TRM_score %>% top_n(99)
bot_75<-TRM_score %>% top_n(-296)
top_25$score="25%"
bot_75$score="75%"
TRM.T_B<-rbind(top_25,bot_75)

cat("Extract survival information from clinical data")
TRM.T_B$OS<-TNBC_Basal$T[match(TNBC_Basal$pid, TRM.T_B$patient_id)]
TRM.T_B$os_status<-TNBC_Basal$DeathBreast[match(TNBC_Basal$pid, TRM.T_B$patient_id)]
TNBC_Basal$DFS<-ifelse(TNBC_Basal$TLR>TNBC_Basal$TDR,TNBC_Basal$TDR,TNBC_Basal$TLR)
TNBC_Basal$dfs_status<-ifelse(c(TNBC_Basal$LR==0 & TNBC_Basal$DR==0),0,1) 
TRM.T_B$DFS<-TNBC_Basal$DFS[match(TRM.T_B$patient_id,TNBC_Basal$pid)]
TRM.T_B$dfs_status<-TNBC_Basal$dfs_status[match(TRM.T_B$patient_id,TNBC_Basal$pid)]

cat("Generate survival plot\n")

library("survival")
library("survminer")

TRM.T_B$os_status<-as.numeric(TRM.T_B$os_status)


os_fit <- survfit(Surv(OS, os_status) ~ score, data = TRM.T_B)
dfs_fit <- survfit(Surv(DFS, dfs_status) ~ score, data = TRM.T_B)

os_plot<-ggsurvplot(os_fit,
          pval = TRUE, conf.int = FALSE, ylab = "Overall survival",
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          palette = c("#E7B800", "#2E9FDF"))

dfs_plot <-ggsurvplot(dfs_fit,
          pval = TRUE, conf.int = FALSE, ylab = "Disease-free survival",
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          palette = c("#E7B800", "#2E9FDF"))

ggsave(print(os_plot), file= "c:/Users/jitan/Documents/TRM.T_B_os.jpeg")
ggsave(print(dfs_plot), file= "c:/Users/jitan/Documents/TRM.T_B_dfs.jpeg")