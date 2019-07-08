# This script is for generating sig.score in TNBC/Basal samples
# Jiayi Tan
# 06/24/2019

#rm(list = ls(all.names = TRUE))

test_func = function (clin, subtype, coloi) {
#	clin_oi = subset(clin, Pam50Subtype == subtype)
	temp_mask = clin[,coloi] == subtype
	print(head(temp_mask))
	clin_oi = clin[temp_mask,]
	print(dim(clin_oi))
	return(clin_oi)
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

used_clin = test_func(clin = clin_info, subtype = "LumA", coloi = "Pam50Subtype")
stop("JUST STOP!!!!!!!")
if (FALSE) {
cat("Select TNBC/Basal samples from clinical data frame\n")
Basal_info<-subset(clin_info, Pam50Subtype =="Basal")
TNBC_info<-subset(clin_info, er_status =="-"& her2_status =="-"& pr_status =="-")
TNBC_Basal<-rbind(Basal_info,TNBC_info)
TNBC_Basal<-unique(TNBC_Basal)
IDC_info<-subset(clin_info, oncotree_code =="IDC")
library(dplyr)
TNBC_Basal.IDC<-intersect(TNBC_Basal,IDC_info)
TNBC.Basal<-intersect(TNBC_info,Basal_info)

cat("Generate expression data from TNBC/Basal/IDC samples\n")
data.meta_TNBC<-subset(data.meta,rownames(data.meta) %in% TNBC_info$pid)
data.meta_Basal<-subset(data.meta,rownames(data.meta) %in% Basal_info$pid)
data.meta_T_B<-subset(data.meta,rownames(data.meta) %in% TNBC_Basal$pid)
data.meta_T.B<-subset(data.meta,rownames(data.meta) %in% TNBC.Basal$pid)
data.meta_T_B.IDC<-subset(data.meta,rownames(data.meta) %in% TNBC_Basal.IDC$pid)


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
TRM_T.score <- sig.score(x=TRM_sig, data=data.meta_TNBC, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)

TRM_B.score <- sig.score(x=TRM_sig, data=data.meta_Basal, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)

TRM_TB.score <- sig.score(x=TRM_sig, data=data.meta_T.B, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)

TRM_IDC.score <- sig.score(x=TRM_sig, data=data.meta_T_B.IDC, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)


TRM_score<- as.data.frame (TRM.score$score)
library(tibble)
TRM_score<-rownames_to_column(TRM_score,var ="patient_id")
colnames(TRM_score)[2]="sig.score"

test<- as.data.frame (TRM.score$score)
library(tibble)
test1<-rownames_to_column(test, var ="patient_id")
test2<-test1[order(-test1$sig.score),]
library(dplyr)
top_25<-test2 %>% top_n(93)
bot_75<-test2 %>% top_n(-279)
top_25$score="25%"
bot_75$score="75%"
TNBC_Basal1<-rbind(top_25,bot_75)
TNBC_Basal1$os<-TNBC_Basal$os_time[match(TNBC_Basal$pid, TNBC_Basal1$patient_id)]
TNBC_Basal1$status<-TNBC_Basal$os[match(TNBC_Basal$pid, TNBC_Basal1$patient_id)]
colnames(TNBC_Basal1)[4]="time"
TNBC_Basal1$status<-c("LIVING"=1,"DECEASED"=2)
fit <- survfit(Surv(time, status) ~ score, data = TNBC_Basal1)
print(fit)
summary(fit)$table


library(ggplot2)
ggplot(TRM.score, aes(x=sig.score)) + geom_histogram()}