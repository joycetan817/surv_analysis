# This script is for survaval analysis of subtype samples with defined function
# Jiayi Tan
# 06/28/2019

rm(list = ls(all.names = TRUE))
sub_clin = function (clin, subtype, coloi) {
	temp_mask = clin[,coloi] == subtype
  IDC_info<-subset(clin_info, oncotree_code =="IDC")
	clin_oi = intersect(IDC_info, clin[temp_mask,])
}


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



subtype_clin=sub_clin(clin = clin_info, subtype = "LumA", coloi = "Pam50Subtype")
subtype = "LumA"


cat("Loading METABRIC annotation data...\n")
library(readxl)
annot.meta<-read_excel("HumanHT_12_v30_R3_cleaned.xlsx", sheet = 1)
annot.meta<-annot.meta[c(5,9,14)]
colnames(annot.meta)[1]="NCBI.gene.symbol"
colnames(annot.meta)[2]="EntrezGene.ID"
colnames(annot.meta)[3]="probe"

cat("Load exhuasted T cell signature data\n")
library(xlsx)
ET<-read.xlsx(file="trm_tex_brtissue_only_all_markers.xlsx",4, header = TRUE)
ET<-ET[c(1,3)]
colnames(ET)[1]="NCBI.gene.symbol"
colnames(ET)[2]="coefficient"
ET$EntrezGene.ID<-annot.meta$EntrezGene.ID[match(ET$NCBI.gene.symbol, annot.meta$NCBI.gene.symbol)]
ET$probe<-annot.meta$probe[match(ET$NCBI.gene.symbol, annot.meta$NCBI.gene.symbol)]

cat("Run ET gene sig.score in BC samples\n")
ET.score <- sig.score(x=ET, data=data.meta, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)

cat("Generate subtype signature score and survival data\n")
ET_score<- as.data.frame (ET.score$score)
library(tibble)
ET_score<-rownames_to_column(ET_score,var ="patient_id")
colnames(ET_score)[2]="sig.score"

sub_score<-subset(ET_score, ET_score$patient_id %in% subtype_clin$pid)

sub_score<-sub_score[order(-sub_score$sig.score),]
sub_score$OS<-subtype_clin$T[match(subtype_clin$pid, sub_score$patient_id)]
sub_score$os_status<-subtype_clin$DeathBreast[match(subtype_clin$pid, sub_score$patient_id)]
subtype_clin$DFS<-ifelse(subtype_clin$TLR>subtype_clin$TDR,subtype_clin$TDR,subtype_clin$TLR)
subtype_clin$dfs_status<-ifelse(c(subtype_clin$LR==0 & subtype_clin$DR==0),0,1) 
sub_score$DFS<-subtype_clin$DFS[match(sub_score$patient_id,subtype_clin$pid)]
sub_score$dfs_status<-subtype_clin$dfs_status[match(sub_score$patient_id,subtype_clin$pid)]

print(dim(sub_score))

stop("JUST STOP!!!!!!!")

library(dplyr)
top_25<-subset(sub_score, sig.score > quantile(sig.score, prob = 1 - 25/100))
bot_25<-subset(sub_score, sig.score < quantile(sig.score, prob = 25/100))
top_25$score="t25%"
bot_25$score="b25%"
sub_score<-rbind(top_25,bot_25)





filename<- paste("c:/Users/jitan/Documents/ET_", subtype, ".xlsx", sep = "")
library(xlsx)
write.xlsx(sub_score, filename)


cat("Generate histogram plot of signature score\n")

title<-paste("exhausted T cell sig.score in ", subtype, " samples", sep = "")

library(ggplot2)
hist<-ggplot(sub_score, aes(x=sig.score)) + geom_histogram(color="darkblue", fill="lightblue") +
  labs(title=title, x="sig.score", y = "Count")+
  theme_classic()

filename<- paste("c:/Users/jitan/Documents/hist_", subtype, "_ET.jpeg", sep = "")

ggsave(hist, file= filename)


cat("Generate survival plot\n")
library("survival")
library("survminer")

sub_score$os_status<-as.character(sub_score$os_status)
sub_score$os_status<-as.numeric(sub_score$os_status)
os_fit <- survfit(Surv(OS, os_status) ~ score, data = sub_score)


os_plot<-ggsurvplot(os_fit,
          pval = TRUE, conf.int = FALSE, ylab = "Overall survival",
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          palette = c("#E7B800", "#2E9FDF"))

filename<- paste("c:/Users/jitan/Documents/", subtype, "_ET_os.jpeg", sep = "")

ggsave(file=filename, print(os_plot))


dfs_fit <- survfit(Surv(DFS, dfs_status) ~ score, data = sub_score)
dfs_plot <-ggsurvplot(dfs_fit,
          pval = TRUE, conf.int = FALSE, ylab = "Disease-free survival",
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          palette = c("#E7B800", "#2E9FDF"))

filename<- paste("c:/Users/jitan/Documents/", subtype, "_ET_dfs.jpeg", sep = "")

ggsave(file=filename, print(dfs_plot))

