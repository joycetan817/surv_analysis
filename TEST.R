## This script is for Cox regression analysis 
# Jiayi Tan
# 06/28/2019










subgroup_score=function(subgroup){
    sub_clin<-if (subgroup=="Basal"){
        subset(clin_info, Pam50Subtype =="Basal")
    } else if (subgroup =="TNBC"){
        subset(clin_info, er_status =="-"& Nor_status =="-"& pr_status =="-")
    }
    subgroup_score<-subset(TRM_score,TRM_score$patient_id %in% sub_clin$pid)
}





p <- ggplot(sub_score, aes(grade, ET.sig))
p + geom_point()+stat_compare_means(method = "anova")




IDC_info<-subset(clin_info, oncotree_code =="IDC")
library(dplyr)
TNBC_Basal.IDC<-intersect(TNBC_Basal,IDC_info)
TNBC.Basal<-intersect(TNBC_info,Basal_info)

data.meta_TNBC<-subset(data.meta,rownames(data.meta) %in% TNBC_info$pid)
data.meta_Basal<-subset(data.meta,rownames(data.meta) %in% Basal_info$pid)

data.meta_T.B<-subset(data.meta,rownames(data.meta) %in% TNBC.Basal$pid)
data.meta_T_B.IDC<-subset(data.meta,rownames(data.meta) %in% TNBC_Basal.IDC$pid)

TRM_T.score <- sig.score(x=TRM_sig, data=data.meta_TNBC, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)

TRM_B.score <- sig.score(x=TRM_sig, data=data.meta_Basal, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)

TRM_TB.score <- sig.score(x=TRM_sig, data=data.meta_T.B, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)

TRM_IDC.score <- sig.score(x=TRM_sig, data=data.meta_T_B.IDC, annot=annot.meta, do.mapping=FALSE, signed=TRUE, verbose=TRUE)




stop("JUST STOP!!!!!!!")

print(dim(sub_score))


singene_expr = function (gene) {
+     gene_info<-subset(annot.meta, NCBI.gene.symbol == gene) # extract gene probe ID
+     gene_probe<-select(data.expr, patient_id, gene_info$probe) # extract gene expression data in patient sample
+     gene_probe$gene<-rowMeans(gene_probe[2:length(gene_probe)]) # generate mean value of gene expression when there is more than one probe
+     sub_score$gene<-gene_probe$gene[match(sub_score$patient_id, gene_probe$patient_id)] # extract gene expression to subtype sig.score data frame
+     colnames(sub_score)[length(sub_score)] <- gene # assign the column names with gene name
+     cat("Run Cox regression survival analysis by other gene expression signature adjusted for age, nodal, grade and tumor size\n")
+     res.cox <- coxph(Surv(DFS, dfs_status) ~ gene + age + grade + size + nodal_status, data =  sub_score)
+     summary(res.cox)
+     return(sub_score)
+ } 
> sub_score<-singene_expr(gene = "CD8A")


cat("Run correlation analysis between TRM and other genes\n")
sub_cor<-sub_score[,c("sig.score", "CD8A")]  
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
jpeg("test.jpeg")
cor_plot<-corrplot(M, type="upper",
         p.mat = p.mat, sig.level = 0.01) ## Specialized the insignificant value according to the significant level
dev.off()

if(cox_reg) {
	res.cox <- coxph(Surv(rfst, rfse) ~ gene+age+grade+tsize +node_stat, data = sub_scres)
	prescox = summary(res.cox) # Print results of COX
		sum_txt = paste(gene, "_summary_cox_results.txt", sep = "")
		sink(sum_txt)
		print(summary(res.cox))
		sink()
		cox_rds <- paste(gene, "_survana_cox_res.rds", sep="")
		cat("\t\tP-value from logRank test: ", (prescox$sctest["pvalue"]), "\n") 
		saveRDS(res.cox, file=cox_rds)
	}

if(cor_plot != "") {
	sub_score<-sub_scres[complete.cases(sub_scres),]
	res.cox <- coxph(Surv(rfst, rfse) ~ gene+age+grade+tsize+node_stat, data = sub_score)
	prescox = summary(res.cox) # Print results of COX
		sum_txt = paste(gene, "_summary_cox_results.txt", sep = "")
		sink(sum_txt)
		print(summary(res.cox))
		sink()
		cox_rds <- paste(gene, "_survana_cox_res.rds", sep="")
		cat("\t\tP-value from logRank test: ", (prescox$sctest["pvalue"]), "\n") 
		saveRDS(res.cox, file=cox_rds)
	}

sub_clin=subset(clin_info, ER.Expr == "-" & PR.Expr == "-" & Her2.Expr == "-")
sum(sub_clin$CT=="AC")
sum(sub_clin$CT=="AC/CMF")
sum(sub_clin$CT=="APD")
 sum(sub_clin$CT=="CAPE")
 sum(sub_clin$CT=="CEF")
 sum(sub_clin$CT=="CMF")
 sum(sub_clin$CT=="ECMF")
 sum(sub_clin$CT=="FAC")
 sum(sub_clin$CT=="HER")
 sum(sub_clin$CT=="PACL")
 sum(sub_clin$CT=="TAXOID")
 sum(sub_clin$CT=="NO/NA")
 sum(sub_clin$CT=="OTHER")
 sum(sub_clin$CT=="CHEMO")
 sum(sub_clin$RT=="CW")
 sum(sub_clin$RT=="CW_NODAL")
 sum(sub_clin$RT=="Ext R/T")
 sum(sub_clin$RT=="Ext  R/T")
 sum(sub_clin$RT=="Unspecified R/T; Ext  R/T")
 sum(sub_clin$RT=="Unspecified R/T; Ext R/T")
 sum(sub_clin$RT=="HD Caesium implant alone; Ext  R/T")
 sum(sub_clin$RT=="HD Caesium implant alone; Ext R/T")
 sum(sub_clin$RT=="NODAL")
 sum(sub_clin$RT=="UNK SITE")
 sum(sub_clin$RT=="Unspecified R/T")
 sum(sub_clin$RT=="Nnne")
 sum(sub_clin$RT=="HD Caesium implant alone")
 sum(sub_clin$RT=="Ext R/T +Iridium implant")
 sum(sub_clin$RT=="NO/NA")
 sum(sub_clin$RT=="NONE RECORDED IN LANTIS")
 sum(sub_clin$HT=="TAM")
 sum(sub_clin$HT=="AI")
 sum(sub_clin$HT=="TAM/AI")
 sum(sub_clin$HT=="OO")
 sum(sub_clin$HT=="GNRHA")
 sum(sub_clin$HT=="TAM/ZOLADEX")
 sum(sub_clin$HT=="ATAC")
 sum(sub_clin$HT=="MEGACE")
 sum(sub_clin$HT=="NO/NA")
 sum(sub_clin$HT=="OTHER")
 sum(sub_clin$HT=="ENDOCRINE")
 sum(sub_clin$HT=="Y")

singene_expr = function (gene, expr, annot, subdf, caltype = "mean") {
	# caltype: calculation type for multiple probes, mean, median, max, min
	gene_info = subset(annot, ILMN_Gene == gene) # ILMN_Gene is the column of gene name which is used to extract gene probe ID
	if (dim(gene_info)[1] == 0) {stop("No MATCHED gene found!!!")}
#	print(gene_info)
	if (sum(gene_info$probe %in% rownames(expr)) == 0) {stop("NO MATCHED probe found!!!")}
	gene_probe = expr[gene_info$probe,]

	if (caltype == "mean") {
		sub_expr = as.data.frame(apply(gene_probe, 2, FUN = mean))
		colnames(sub_expr) = gene 
	}
	if (caltype == "median") {
		sub_expr = as.data.frame(apply(gene_probe, 2, FUN = median))
		colnames(sub_expr) = gene 
	}
	if (caltype == "max") {
		sub_expr = as.data.frame(apply(gene_probe, 2, FUN = max))
		colnames(sub_expr) = gene 
	}
	if (caltype == "min") {
		sub_expr = as.data.frame(apply(gene_probe, 2, FUN = min))
		colnames(sub_expr) = gene 
	}

	sub_expr$pid = rownames(sub_expr)
	sub_res = merge(subdf, sub_expr, by = "pid")
	rownames(sub_res) = sub_res$pid
	if (dim(sub_res)[1] == dim(subdf)[1]) {
		cat("\tAll the patients are matched!\n")
	} else {
		cat("\tWarning: missed patients!!!\n")
	}
	return(sub_res)
}
histype = "IDC" # histology type: IDC/DCIS
pamst = "" # PAM50 status: LumA/LumB/Basal/Normal/Her2
gdoi = 0 #c(1) # Grade of interest: 1/2/3
hrtype = c("N", "N", "N") # N: Negative, P: Positive, "-": DON'T CARE
sig_save = FALSE
gp_app = "symqcut" # oneqcut: one quantile cutoff, symqcut: symmetric quantile cutoff
qcut = 0.25 # This is TOP quantile for oneqcut approach
gp_gene = "" # Group gene used for categorizing the cohort(if run cox regression of single gene)
# Default "": use signature score 
corr_gene = c("CD8A", "CD3G", "ITGAE", "STAT1") # Genes need to be correlated with signature scores
gptype = "Tex sig.score"



sssign = as.data.frame(matrix(ncol = 3, nrow = dim(sub_annot)[1]))

test = as.data.frame(matrix(ncol = 8, nrow = length(unique(clin_info$CT))))
rownames(test)<-c(unique(clin_info$CT))
colnames(test)<-c("LumA", "LumB", "Her2", "Basal", "Normal", "ER+Her2+", "ER+Her2-","TNBC")
"NO/NA"  "OTHER"  "AC"     "CMF"    "ECMF"   "CAPE"  
 [7] "AC/CMF" "PACL"   "FAC"    "CEF"    "APD"    "CHEMO" 
[13] "TAXOID" "HER" 


for (ig in 1:length(corr_gene)) {
	sub_corr = singene_expr(gene = corr_gene[ig], expr = expr, annot = annot, subdf = sub_corr)
}


BC_treatment = function(Treat_type, subtype) {
	test[1,1]=dim(subset(sub_clin, CT ==chemotype[1])
}
chemotype = unique(clin_info$CT)
chemotype = c("NO/NA", "OTHER", "AC", "CMF", "ECMF", "CAPE", "AC/CMF", "PACL", "FAC", "CEF","APD","CHEMO", "TAXOID", "HER")
for (ct in 1:length(chemotype)) {
	T[1]<-as.numeric(dim(subset(sub_clin, CT ==chemotype[ct]))[1])
}
test[1]<-list

for (ct in 1:length(chemotype)) {
list[[ct]]<-dim(subset(sub_clin, CT ==chemotype[ct]))[1]
}
sum(subset(sub_clin, CT =="NO/NA"))

sub_clin<-subset(clin_info, oncotree_code == "IDC" & Pam50Subtype == "Normal")
for (ht in 1:length(endotype)) {
list[[ht]]<-dim(subset(sub_clin, HT == endotype[ht]))[1]
}
HT[5]<-list