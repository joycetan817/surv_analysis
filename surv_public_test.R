# This script is for survaval analysis according to any customized signature score of subtype samples
# Weihua Guo, Ph.D.
# 7/17/2019
# Jiayi Tan
# 07/01/2019

rm(list = ls(all.names = TRUE))
sub_clin = function (clin, subtype, coloi) {
	temp_mask = clin[,coloi] == subtype
  IDC_info<-subset(clin_info, oncotree_code =="IDC")
	clin_oi = intersect(IDC_info, clin[temp_mask,])
}

# Please load all the packages at the very begining of each script
cat("Loading genefu library...\n")
suppressMessages(library(genefu))
suppressMessages(library(readxl))
suppressMessages(library(ggplot2))


# Personally prefer to having a independent section for all the parameters or variables
# which may be tuned for different inputs and tasks

# For each variable, please at least leave a simple comments describing what this 
# variable represent for

work_dir = "/home/weihua/mnts/group_plee/Weihua/surv_validation/" # working directory/path for survival validation
db_name = "metabric"
sg_name = "loi_trm"
# data_dir = "/home/weihua/mnts/group_plee/Weihua/metabric_use/" # directory/path for public data
data_dir = paste(work_dir, db_name, "/", sep = "") # generate the directory with all the public data
sign_dir = paste(work_dir, sg_name, "/", sep = "") # generate the directory with signatures and corresponding results

expr_file = "metabric_expr_ilid.RDS" # Expression file
clin_file = "merge_clin_info.xlsx" # clinical information
annot_file = "HumanHT_12_v30_R3_cleaned_v2.xlsx" # Microarray/Genome annotation
# sign_file = "trm_tex_brtissue_only_all_markers.xlsx" # Signature file
sign_file = "loi_trm_signature.txt" # Signature file

pamst = "LumA"
gp_app = "oneqcut"
qcut = 0.25


cat("Loading METABRIC expression data...\n")
st = Sys.time()
## Please use either the full path of the file or change the work directory here
expr = readRDS(paste(data_dir, expr_file, sep = ""))
print(Sys.time()-st)
# print(meta_expr[1:9,1:6]) # Check the input in terminal
# data.meta<-t(meta_expr)

cat("Loading clinical data...\n")
st = Sys.time()
# clin_info<-readRDS("merge_clin_info_manual_checked.RDS") ### I canNOT find this file
clin_info = as.data.frame(read_excel(paste(data_dir, clin_file, sep = ""),1))
# print(head(clin_info))
print(Sys.time()-st)

cat("Start to filter by clinical info...\n")
sub_clin = clin_info
cat("\tOriginal patient number: ", dim(sub_clin)[1], "\n")
# For IDC
sub_clin = sub_clin[sub_clin[,"oncotree_code"] == "IDC",]
sub_clin = sub_clin[complete.cases(sub_clin$pid),]
cat("\tFiltered patient number: ", dim(sub_clin)[1], "\n")
# For LumA
sub_clin = sub_clin[sub_clin[,"Pam50Subtype"] == pamst,]
sub_clin = sub_clin[complete.cases(sub_clin$pid),]
cat("\tFiltered patient number: ", dim(sub_clin)[1], "\n")
## Since my filters are not too many and not complex, I can directly write one line for each filter
# subtype_clin=sub_clin(clin = clin_info, subtype = "LumA", coloi = "Pam50Subtype")
# subtype = "LumA"


cat("Loading METABRIC annotation data...\n")
# library(readxl)
# annot.meta<-read_excel("HumanHT_12_v30_R3_cleaned.xlsx", sheet = 1)
st = Sys.time()
annot = as.data.frame(read_excel(paste(data_dir, annot_file, sep = "")))
# annot = read.table(paste(data_dir, annot_file, sep = ""), header = TRUE)
# print(head(annot))
print(Sys.time()-st)

## Save this for later
if (FALSE) {
annot.meta<-annot.meta[c(5,9,14)]
colnames(annot.meta)[1]="NCBI.gene.symbol"
colnames(annot.meta)[2]="EntrezGene.ID"
colnames(annot.meta)[3]="probe"}
cat("Load gene signature...\n")
# library(xlsx)
sign = read.table(paste(sign_dir, sign_file, sep = ""), header = TRUE, row.names = 1)
# ET<-read.xlsx(file="trm_tex_brtissue_only_all_markers.xlsx",4, header = TRUE)
# print(head(sign))

## PREP for sig.score
# data = [r:sample, c:gene/probe]
# annot = [r:gene/probe, c: genome annotation] 
# x = [r: genes in signature, c: probe+EntrezGene.ID, coefficient]
ssdata = t(expr) # sig.score data
cat("Expression data dimension, sample: probe ", dim(ssdata), "\n")
# print(ssdata[1:9,1:6])

ssannot = annot # sig.score annotation
rownames(ssannot) = ssannot$Probe_Id
ssannot$Probe_Id = NULL
ssannot[,"RealEntrezGene.ID"] = ssannot[,"EntrezGene.ID"]
ssannot[,"EntrezGene.ID"] = ssannot[,"ILMN_Gene"]
cat("Annotation data dimension, probe: annot ", dim(ssannot), "\n")
# print(all(ssannot$probe == rownames(ssannot))) # Check the probe is exactly equal to rownames
# print(head(ssannot))

sign_gene = toupper(rownames(sign))
cat("Signature gene number: ", length(sign_gene),"\n")
anno_gene = annot$ILMN_Gene
ol_gene = intersect(sign_gene, anno_gene)
unava_gene = setdiff(sign_gene, ol_gene)
cat("Unavailable",length(unava_gene) ," genes: ", unava_gene, "\n")
sub_annot = ssannot[ssannot$ILMN_Gene %in% sign_gene,]
cat("Available probe: ", dim(sub_annot)[1], "\n")
sssign = as.data.frame(matrix(ncol = 3, nrow = dim(sub_annot)[1])) # Initialize sig.score x
colnames(sssign) = c("probe", "EntrezGene.ID", "coefficient")
rownames(sssign) = rownames(sub_annot)
sssign$probe = rownames(sub_annot)
sssign$EntrezGene.ID = sub_annot$ILMN_Gene
# sssign$EntrezGene.ID = sub_annot$EntrezGene.ID

sssign$coefficient = sign[sssign$EntrezGene.ID,"logfc"]
# sssign$coefficient = sign[sub_annot$ILMN_Gene,"logfc"]

# print(head(sssign))

ssin = list("data" = ssdata, "annot" = ssannot, "x" = sssign)
# saveRDS(ssin, file = paste(sign_dir, "sig_score_inputs.RDS", sep = ""))

## RUN sig.score
cat("Run sig.score in all patient samples\n")
sig_score <- sig.score(x=sssign, data=ssdata, annot=ssannot, do.mapping=TRUE, signed=TRUE, verbose=TRUE)
# print(head(sig_score$score))

# saveRDS(sig_score, file = paste(sign_dir, "sig_score_outputs.RDS", sep = ""))
cat("Generate subtype signature score and survival data\n")
sc_res = as.data.frame(sig_score$score)
colnames(sc_res) = "sig_score"
sc_res$pid = rownames(sc_res)
print(head(sc_res))

print(head(sub_clin))
print(dim(sub_clin))


# print(sub_clin$pid)
print(head(sc_res))

sub_scres = sc_res[sc_res$pid %in% sub_clin$pid,]

## Histogram for sig.score
cat("Generate histogram plot of signature score\n")
title = paste(db_name, sg_name, pamst, "IDC", sep = "  ")
sc_hist = ggplot(sub_score, aes(x=sig_score)) + 
	geom_histogram(color="darkblue", fill="lightblue") +
	labs(title=title, x="sig.score", y = "Count") + 
	theme_classic()

hist_tif = paste(sign_dir, db_name, sg_name, pamst, "IDC.tiff", sep = "_")

ggsave(sc_hist, file = hist_tif)

print(dim(sub_scres))
q(save = "no")

# ET<-ET[c(1,3)]
# colnames(ET)[1]="NCBI.gene.symbol"
# colnames(ET)[2]="coefficient"
# ET$EntrezGene.ID<-annot.meta$EntrezGene.ID[match(ET$NCBI.gene.symbol, annot.meta$NCBI.gene.symbol)]
# ET$probe<-annot.meta$probe[match(ET$NCBI.gene.symbol, annot.meta$NCBI.gene.symbol)]

# cat("Run ET gene sig.score in all patient samples\n")
# ET.score <- sig.score(x=ET, data=data.meta, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)

# cat("Generate subtype signature score and survival data\n")
# ET_score<- as.data.frame (ET.score$score)
# library(tibble)
# ET_score<-rownames_to_column(ET_score,var ="patient_id")
# colnames(ET_score)[2]="sig.score"

cat("Generate exhuasted T cell sig. score in subtype\n")
sub_score<-subset(ET_score, ET_score$patient_id %in% subtype_clin$pid)

cat("Extract survival information from clinical data to subtype sig.score data frame")
sub_score<-sub_score[order(-sub_score$sig.score),]
sub_score$OS<-subtype_clin$T[match(subtype_clin$pid, sub_score$patient_id)]
sub_score$os_status<-subtype_clin$DeathBreast[match(subtype_clin$pid, sub_score$patient_id)]
subtype_clin$DFS<-ifelse(subtype_clin$TLR>subtype_clin$TDR,subtype_clin$TDR,subtype_clin$TLR)
subtype_clin$dfs_status<-ifelse(c(subtype_clin$LR==0 & subtype_clin$DR==0),0,1) 
sub_score$DFS<-subtype_clin$DFS[match(sub_score$patient_id,subtype_clin$pid)]
sub_score$dfs_status<-subtype_clin$dfs_status[match(sub_score$patient_id,subtype_clin$pid)]

stop("JUST STOP!!!!!!!")

cat("Select patient samples by 25/25 cutoff of exhuasted T sig.score")

library(dplyr)
top_25<-subset(sub_score, sig.score > quantile(sig.score, prob = 1 - 25/100))
bot_25<-subset(sub_score, sig.score < quantile(sig.score, prob = 25/100))
top_25$score="t25%"
bot_25$score="b25%"
sub_score<-rbind(top_25,bot_25)



filename<- paste("c:/Users/jitan/Documents/ET_", subtype, ".xlsx", sep = "")
library(xlsx)
write.xlsx(sub_score, filename)




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


sub_score$grade<-clin_info$grade[match(sub_score$patient_id, clin_info$pid)]
sub_score$grade<-as.character(sub_score$grade)
sub_score$grade<-as.numeric(sub_score$grade)
sub_score<-na.omit(sub_score)

sub_score$grade<-as.factor(sub_score$grade)

cat("Perform ANOVA test by BC grade in BC subtype based on exhuasted T cell sig.score\n")
p <- ggplot(sub_score, aes(x = grade, y = sig.score, color = grade))
p + geom_point()+stat_compare_means(method = "anova") + theme_bw()

cat("Perform t-test between each two groups based on BC grade\n")
my_comparisons <- list( c("1", "2"), c("1", "3"), c("2", "3") )
ggplot(sub_score, aes(x=grade, y=sig.score, color = grade)) + geom_point() + 
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = 10) + theme_bw()
