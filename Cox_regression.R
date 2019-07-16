## This script is for Cox regression analysis 
# Jiayi Tan
# 06/28/2019

rm(list = ls(all.names = TRUE))
sub_clin = function (clin, subtype, coloi) {
	temp_mask = clin[,coloi] == subtype
	clin_oi = clin[temp_mask,]
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

subtype_clin=sub_clin(clin = clin_info, subtype = "Basal", coloi = "Pam50Subtype")
subtype = "Basal"

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
TRM.score <- sig.score(x=TRM_sig, data=data.meta, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)

cat("Generate subtype signature score and survival data\n")
TRM_score<- as.data.frame (TRM.score$score)
library(tibble)
TRM_score<-rownames_to_column(TRM_score,var ="patient_id")
colnames(TRM_score)[2]="sig.score"

cat("Rank data by 25/75 cutoff of TRM sig.score")
sub_score<-subset(TRM_score, TRM_score$patient_id %in% subtype_clin$pid)
sub_score<-sub_score[order(-sub_score$sig.score),]
library(dplyr)
top_25<-subset(sub_score, sig.score > quantile(sig.score, prob = 1 - 25/100))
bot_75<-subset(sub_score, sig.score < top_25$sig.score)
top_25$score="25%"
bot_75$score="75%"
sub_score<-rbind(top_25,bot_75)

cat("Extract survival information from clinical data to subtype sig.score data frame\n")
sub_score$OS<-subtype_clin$T[match(subtype_clin$pid, sub_score$patient_id)]
sub_score$os_status<-subtype_clin$DeathBreast[match(subtype_clin$pid, sub_score$patient_id)]
subtype_clin$DFS<-ifelse(subtype_clin$TLR>subtype_clin$TDR,subtype_clin$TDR,subtype_clin$TLR)
subtype_clin$dfs_status<-ifelse(c(subtype_clin$LR==0 & subtype_clin$DR==0),0,1) 
sub_score$DFS<-subtype_clin$DFS[match(sub_score$patient_id,subtype_clin$pid)]
sub_score$dfs_status<-subtype_clin$dfs_status[match(sub_score$patient_id,subtype_clin$pid)]
colnames(sub_score)[2]="TRM"
sub_score$age<-clin_info$age_at_diagnosis[match(sub_score$patient_id, clin_info$pid)]
sub_score$grade<-clin_info$grade[match(sub_score$patient_id, clin_info$pid)]
sub_score$size<-clin_info$tumor_size[match(sub_score$patient_id, clin_info$pid)]
sub_score$nodal_status<-clin_info$Lymph.Nodes.Positive[match(sub_score$patient_id, clin_info$pid)]
sub_score$size<-as.numeric(sub_score$size)

cat("Run Cox regression survival analysis by TRM sig.score adjusted for age, nodal, grade and tumor size\n")
res.cox <- coxph(Surv(DFS, dfs_status) ~ TRM + age + grade + size + nodal_status, data =  sub_score)

data.expr<-as.data.frame(data.meta)
library(tibble)
data.expr<-rownames_to_column(data.expr, var ="patient_id")


cat("Run Cox regression survival analysis by other gene expression signature adjusted for age, nodal, grade and tumor size\n")
Cox_control = function (gene) {
	gene_info<-subset(annot.meta, NCBI.gene.symbol == gene)
	gene_probe<-select(data.expr, patient_id, gene_info$probe)
	gene_probe$gene<-rowMeans(gene_probe[2:length(gene_probe)]) # generate mean value of gene expression when more than one probe
	sub_score$gene<-gene_probe$gene[match(sub_score$patient_id, gene_probe$patient_id)]
	colnames(sub_score)[length(sub_score)] <- gene 
	return(sub_score)
}

sub_score<-Cox_control(gene = c("ITGAE", "CD8A", "CD3G", "STAT1"))

res.cox <- coxph(Surv(DFS, dfs_status) ~ gene + age + grade + size + nodal_status, data =  sub_score)
summary(res.cox)

cat("Run survival analysis in subtype samples according to CD8A expression\n")
library("survival")
library("survminer")

CD8A_score<-sub_score[order(-sub_score$CD8A),]
top25<-subset(CD8A_score, CD8A > quantile(CD8A, prob = 1 - 25/100))
bot75<-subset(CD8A_score, CD8A < top25$CD8A)
top25$CD8_score="25%"
bot75$CD8_score="75%"
CD8A_score<-rbind(top25,bot75)

os_fit <- survfit(Surv(OS, os_status) ~ CD8_score, data = CD8A_score)
os_plot<-ggsurvplot(os_fit,
          pval = TRUE, conf.int = FALSE, ylab = "Overall survival",
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          palette = c("#E7B800", "#2E9FDF"))

dfs_fit <- survfit(Surv(DFS, dfs_status) ~ CD8_score, data = CD8A_score)
dfs_plot<-ggsurvplot(dfs_fit,
          pval = TRUE, conf.int = FALSE, ylab = "Overall survival",
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          palette = c("#E7B800", "#2E9FDF"))

cat("Run survival analysis in subtype samples according to TRM and CD8A expression\n")

dfs_fit <- survfit(Surv(DFS, dfs_status) ~ CD8_score + score, data = CD8A_score)
dfs_plot<-ggsurvplot(dfs_fit,
            pval = TRUE, conf.int = FALSE, ylab = "Disease-free survival", legend.labs=c("CD8hiTRMhi", "CD8hiTRMlo","CD8loTRMhi","CD8loTRMlo")
         )

cat("Run correlation analysis between TRM and other genes\n")
sub_cor<-sub_score[c(2,12:15)]
M<-cor(sub_cor)

cat("Compute the p-value of correlations\n")
cor.mtest <- function(mat, ...) {
     mat <- as.matrix(mat)
     n <- ncol(mat)
     p.mat<- matrix(NA, n, n)
     diag(p.mat) <- 0
     for (i in 1:(n - 1)) {
         for (j in (i + 1):n) {
             tmp <- cor.test(mat[, i], mat[, j], ...)
             p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
         }
     }
     colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
     p.mat
}
p.mat <- cor.mtest(sub_cor)
head(p.mat)
library(corrplot)
cor_plot<-corrplot(M, 
         p.mat = p.mat, sig.level = 0.01)

