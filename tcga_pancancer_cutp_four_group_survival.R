# This script is for survaval analysis according to any customized signature score of subtype samples in pancancer from TCGA
# Jiayi Tan
# 05/09/2020


#rm(list = ls(all.names = TRUE))

survana <- function(data, gptype, plot = "") {
	# NOTE: plot/csv/cox should be a folder directory where to store the results
	# NOTE: gptype should be the definition or basis of the grouping method, default is signature score.

	cat("Analyzing overall survival...\n")
	surv_data = data	
	numh = sum(surv_data$group == "High")
	numl = sum(surv_data$group == "Low")
	if (dim(surv_data)[1] != numh + numl) {stop("More than two groups!!!")}
#	print(head(surv_data))
	fit <- survfit(Surv(ost, ose) ~ group, data=surv_data)


	if (plot != "") {
		max_ost = max(surv_data$ost)	
		xlim_max = ceiling(max_ost/1000)*1000
		
		surv_plot <- paste(plot, gptype, "_survana_curves.tiff", sep="")
		survcurv <- ggsurvplot(fit, data = surv_data,
				       xlim = c(0,xlim_max), ylim = c(0.00, 1.00),
				       pval = TRUE, pval.size = 6, pval.coord = c(xlim_max/4*3, 0.95),
				       conf.int = FALSE, conf.int.alpha = 0.2,
				       xlab = "Time (days)", ylab = "Overall Survival", legend.title = "",
				       legend.labs = c(paste("High, n =", numh), paste("Low, n =", numl)), legend = c(0.15,0.15),
#                                      surv.median.line = "hv",
				       ggtheme = theme_classic(),
				       palette = c("#D15466", "#0A9CC7"), 
				       font.x = c(14, "bold"), font.y = c(18, "bold"), 
				       font.tickslab = c(16, "plain", "black"), font.legend = c(18, "bold"),
				       risk.table = FALSE, ncensor.plot = FALSE)
		#note_on_plot <- paste("nHigh = ", numh, "\t", "nLow = ", numl, "\t")
		#survcurv$plot <- survcurv$plot + annotate("text", x=2000, y=0.1, label=note_on_plot, size = 6)
		ggsave(surv_plot, plot = survcurv$plot, dpi = 120, width = 9, height = 6, units = 'in')
	}
}

# Extract single gene expression from expression data to subtype sig.score data frame
singene_expr = function (gene, expr, annot, subdf, caltype = "mean", map = TRUE) {
	if (map) {
		# caltype: calculation type for multiple probes, mean, median, max, min
		gene_info = subset(annot, ILMN_Gene == gene) # ILMN_Gene is the column of gene name which is used to extract gene probe ID
		if (dim(gene_info)[1] == 0) {stop("No MATCHED gene found!!!")}
#		print(gene_info)
		if (sum(gene_info$Probe_Id %in% rownames(expr)) == 0) {stop("NO MATCHED probe found!!!")}
		gene_probe = expr[gene_info$Probe_Id,]

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
	} else {
			print(expr[1:9,1:6])
			if (sum(rownames(expr) == gene) == 1) {
				sub_expr = as.data.frame(t(expr[gene, 3:dim(expr)[2]]))
			} else {stop("No matched gene")}
		}

	sub_expr$pid = rownames(sub_expr)
	sub_res = merge(subdf, sub_expr, by = "pid")
	rownames(sub_res) = sub_res$pid
	if (dim(sub_res)[1] == dim(subdf)[1]) {
		cat("\tAll the patients are matched!\n")
	} else {
		cat("\tWarning: missed patients!!!\n")
	}
	print(head(sub_res))
	return(sub_res)

}
# Please load all the packages at the very begining of each script
cat("Loading genefu library...\n")
suppressMessages(library(genefu))
suppressMessages(library(readxl))
suppressMessages(library(ggplot2))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(corrplot))
suppressMessages(library(ggpubr))
suppressMessages(library(Hmisc))
suppressMessages(library(survMisc))
# Personally prefer to having a independent section for all the parameters or variables
# which may be tuned for different inputs and tasks

# For each variable, please at least leave a simple comments describing what this 
# variable represent for

#################################################################################
# work_dir = "/home/weihua/mnts/group_plee/Weihua/surv_validation/" # working directory/path for survival validation
work_dir = "//Bri-net/citi/Peter Lee Group/Joyce/TCGA_pancancer/"
organ = "breast"
db_name = paste("tcga_", organ, sep = "")
#sg_name = "loi_trm" # Loi's TRM sig
sg_name = "tex_brtissue" # Colt's Tex sig from breast tissue c2
#sg_name = "mamma" # mamma sig
#expr_type = "" # ilid: raw data from EGA, median: raw median data from cbioportal, medianz: zscore from cbioportal
selfmap = TRUE # NOTE: ilid/tcga requires this as TRUE; median as FALSE
log2_trans = TRUE #use log2 transformation data or not
mutation_corr = FALSE
subtype = ""


data_dir = paste(work_dir, db_name, "/", sep = "") # generate the directory with all the public data
sign_dir = paste(work_dir, sg_name, "/", sep = "") # generate the directory with signatures and corresponding results

mutation_file = paste(work_dir, "mutation/", subtype, "_varscan/", sep = "")

if (organ == "melanoma") {
	if (log2_trans) {
		expr_file = paste(organ, "_Primary_Tumor_Metastatic_FPKM-UQ_merged_tcga_log2_trans_v2.RDS", sep = "")
	} else {expr_file = paste(organ, "_Primary_Tumor_Metastatic_FPKM-UQ_merged_tcga.RDS", sep = "")}
} else {
	if (log2_trans) {
		expr_file = paste(organ, "_Primary_Tumor_FPKM-UQ_merged_tcga_log2_trans_v2.RDS", sep = "")
	} else {expr_file = paste(organ, "_Primary_Tumor_FPKM-UQ_merged_tcga.RDS", sep = "")}
}
	
#expr_file = "tcga_brca_log2trans_fpkm_uq_v2.RDS" # Expression file
clin_csv = "clinical_v2.csv" # clinical information with merged disease-free survival 
annot_file = "gencode.gene.info.v22.xlsx" # Microarray/Genome annotation


#sign_file = "loi_trm_signature.txt" # Signature file Loi's TRM
sign_file = "tex_signature_colt_c2.txt" # Signature file Colt's Tex
#sign_file = "mamma_signature_v1.txt" # Signature file mamma
#sign_file = "tex_signature_tirosh.txt"
#sign_file = "tex_signature_schumacher.txt"
#sign_file = "amit_melanoma_logfc1.txt" # Signature file Colt's Tex
#sign_file = "tex_guo_lung.txt" # Signature file Colt's Tex
#sign_file = "cd26_c1_all_tumor_v1.txt"
#sign_file = "tcf1_fc25_all_tumor_wt377_0420_cbpt.txt"






gdoi = 0 #c(1) # Grade of interest: 1/2/3
stageoi = 0 # Stage of interest: 1/2/3/4
Tstageoi = 0 #c("T3", "T4") # T stage of interest : T1/T2/T3/T4
histype = "IDC"
hrtype = c("N", "N", "N")# N: Negative, P: Positive, "-": DON'T CARE
sig_save = FALSE
gp_app = "symqcut"#"symqcut" # oneqcut: one quantile cutoff (upper percential), symqcut: symmetric quantile cutoff
qcut = 0.25 #0.25 # This is TOP quantile for oneqcut approach
gp_gene = "CD8A" # Group gene used for categorizing the cohort(if run cox regression of single gene)
# Default "": use signature score 
corr_gene = ""#c("CD274", "PDCD1") #c("CD8A", "CD3G", "ITGAE", "STAT1") # Genes need to be correlated with signature scores
gptype = "Tex sig.score"
trt_type = "" #c("ct", "rt", "ht") # check the correlation between sig.score and treatment
cox_reg = FALSE
group_in_priMarker = FALSE ##second group strata within the primary gene marker group
optimal_cutoff = TRUE


#################################################################################
# Work for experiment records
res_folder = "optimal_tex_CD8A_breast_TNBC_log2_tcga_overdrive" # NOTE: Please change this folder name to identify your experiments
res_dir = paste(sign_dir, res_folder, "/", sep ="")
dir.create(file.path(sign_dir, res_folder), showWarnings = FALSE)
# COPY the used script to the result folder for recording what experiment was run
### !!!Please change the script_dir to the folder directory where this script is located
script_dir = "~/GitHub/surv_analysis/"
script_name = "tcga_pancancer_survival_test.R"
file.copy(paste(script_dir, script_name, sep = ""), res_dir)  

#################################################################################
cat("Loading expression data...\n")
st = Sys.time()
## Please use either the full path of the file or change the work directory here
expr = readRDS(paste(data_dir, expr_file, sep = ""))
#expr = readRDS("metabric_expr_ilid.RDS") # When test the script using metabric
#expr = readRDS("tcga_brca_log2trans_fpkm_uq_v2.RDS") # When test the script using tcga
#expr = readRDS("tcga_brca_rsem_expr_log2trans.RDS")
#expr <- readRDS("primary_tumor_cleaned_merged_raw_counts.rds")
#expr = readRDS("data_expression_median.RDS") # When test the script using cBioportal
#expr = readRDS("tcga_portal_data_expr_v3.RDS")
print(Sys.time()-st)
# print(meta_expr[1:9,1:6]) # Check the input in terminal
expr<-as.data.frame(expr)
#expr<-log2(expr)
#expr[expr== -Inf]=0

cat("Loading clinical data...\n")
st = Sys.time()

clin_info = read.csv(paste(data_dir, clin_csv, sep = ""))

print(Sys.time()-st)

cat("Start to filter by clinical info...\n")
sub_clin = clin_info
cat("\tOriginal patient number: ", dim(sub_clin)[1], "\n")


if (gdoi != 0) {
	# print(head(sub_clin))
	sub_clin = sub_clin[sub_clin[,"grade"] %in% gdoi,]
	sub_clin = sub_clin[complete.cases(sub_clin$pid),]
	cat("\tFiltered patient number: ", dim(sub_clin)[1], "\n")
}

if (stageoi != 0) {
	sub_clin = subset(sub_clin, Stage %in% stageoi)
	sub_clin = sub_clin[complete.cases(sub_clin$pid),]
	cat("\tFiltered patient number: ", dim(sub_clin)[1], "\n")

}

if (Tstageoi != 0) {
	sub_clin = subset(sub_clin, Size_Tstage %in% Tstageoi)
	sub_clin = sub_clin[complete.cases(sub_clin$pid),]
	cat("\tFiltered patient number: ", dim(sub_clin)[1], "\n")

}
if (organ == "breast") {
	if (histype !=  "") {
	# For IDC
	sub_clin = sub_clin[sub_clin[,"oncotree_code"] == histype,]
	#sub_clin = sub_clin[sub_clin[,"Oncotree.Code"] %in% histype,]
	sub_clin = sub_clin[complete.cases(sub_clin$pid),]
	cat("\tFiltered patient number: ", dim(sub_clin)[1], "\n")
	}
	if (hrtype != "") {
			if (length(hrtype) != 3) {stop("Not enough info for hormone receptor status!!!")}
			cat("\tER: ", hrtype[1], "\n")
			if (hrtype[1] == "P") {
				sub_clin = sub_clin[sub_clin$ER.Expr == "+",]
				sub_clin = sub_clin[complete.cases(sub_clin$pid),]
			}
			if (hrtype[1] == "N") {
				sub_clin = sub_clin[sub_clin$ER.Expr == "-",]
				sub_clin = sub_clin[complete.cases(sub_clin$pid),]
			}
			cat("\t\tFiltered patient number: ", dim(sub_clin)[1], "\n")

			cat("\tPR: ", hrtype[2], "\n")
			if (hrtype[2] == "P") {
				sub_clin = sub_clin[sub_clin$PR.Expr == "+",]
				sub_clin = sub_clin[complete.cases(sub_clin$pid),]
			}
			if (hrtype[2] == "N") {
				sub_clin = sub_clin[sub_clin$PR.Expr == "-",]
				sub_clin = sub_clin[complete.cases(sub_clin$pid),]
			}
			cat("\t\tFiltered patient number: ", dim(sub_clin)[1], "\n")

			cat("\tHer2: ", hrtype[3], "\n")
			if (hrtype[3] == "P") {
				sub_clin = sub_clin[sub_clin$Her2.Expr == "+",]
				sub_clin = sub_clin[complete.cases(sub_clin$pid),]
			}
			if (hrtype[3] == "N") {
				sub_clin = sub_clin[sub_clin$Her2.Expr == "-",]
				sub_clin = sub_clin[complete.cases(sub_clin$pid),]
			}
			cat("\t\tFiltered patient number: ", dim(sub_clin)[1], "\n")
	}
}

if (dim(sub_clin)[1] <= 3) {stop("TOO SMALL SAMPLE SIZE!!!")}

cat("Loading genome annotation data...\n")
st = Sys.time()
annot = as.data.frame(read_excel(paste(work_dir, annot_file, sep = ""))) 
# annot = read.table(paste(data_dir, annot_file, sep = ""), header = TRUE)
# print(head(annot))
print(Sys.time()-st)

cat("Load gene signature...\n")
sign = read.table(paste(sign_dir, sign_file, sep = ""), header = TRUE, row.names = 1)


#################################################################################
if (selfmap) {
	## PREP for sig.score with MAPPPING
	# data = [r:sample, c:gene/probe]
	# annot = [r:gene/probe, c: genome annotation] 
	# x = [r: genes in signature, c: probe+EntrezGene.ID, coefficient]
	ssdata = t(expr) # sig.score data
	ssdata = as.data.frame(ssdata)
	cat("Expression data dimension, sample: probe ", dim(ssdata), "\n")
	# print(ssdata[1:9,1:6])

	ssannot = annot # sig.score annotation
	rownames(ssannot) = ssannot$Probe_Id
	ssannot$Probe_Id = NULL
	ssannot[,"RealEntrezGene.ID"] = ssannot[,"EntrezGene.ID"]
	ssannot[,"EntrezGene.ID"] = ssannot[,"ILMN_Gene"] # Trick the code to use EntrezGene.ID for mapping
	cat("Annotation data dimension, probe: annot ", dim(ssannot), "\n")
	# print(all(ssannot$probe == rownames(ssannot))) # Check the probe is exactly equal to rownames
	# print(head(ssannot))

	sign_gene = rownames(sign) 
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

	## RUN sig.score
	cat("Run sig.score in all patient samples\n")
	sig_score <- sig.score(x=sssign, data=ssdata, annot=ssannot, do.mapping=TRUE, signed=TRUE, verbose=TRUE)
	# print(head(sig_score$score))

} else {
	## PREP for sig.score WITHOUT MAPPPING
	# data = [r:sample, c:gene]
	# annot = [r:gene, c: genome annotation] 
	# x = [r: genes in signature, c: probe+EntrezGene.ID, coefficient]
	
	
	rownames(expr) = expr[,1]
	ssdata = expr[, 3:dim(expr)[2]]
	ssdata = as.data.frame(t(ssdata))
	cat("Expression data dimension, sample: probe ", dim(ssdata), "\n")
	print(ssdata[1:9,1:6])

	ssannot = expr[,1:2]
	cat("Annotation data dimension, probe: annot ", dim(ssannot), "\n")
	print(head(ssannot))

	sign_gene = toupper(rownames(sign))  
	cat("Signature gene number: ", length(sign_gene),"\n")
	anno_gene = rownames(expr)
	ol_gene = intersect(sign_gene, anno_gene)
	unava_gene = setdiff(sign_gene, ol_gene)
	cat("Unavailable",length(unava_gene) ," genes: ", unava_gene, "\n")
	sub_annot = ssannot[ssannot$Hugo_Symbol %in% sign_gene,]
	cat("Available probe: ", dim(sub_annot)[1], "\n")
	sssign = as.data.frame(matrix(ncol = 3, nrow = dim(sub_annot)[1])) # Initialize sig.score x
	colnames(sssign) = c("probe", "EntrezGene.ID", "coefficient")
	rownames(sssign) = rownames(sub_annot)
	sssign$probe = rownames(sub_annot)
	sssign$EntrezGene.ID = sub_annot$Entrez_Gene_Id

	sssign$coefficient = sign[sssign$probe,"logfc"]
	cat("Run sig.score in all patient samples\n")
	sig_score = sig.score(x=sssign, data=ssdata, annot=ssannot, do.mapping=FALSE, signed=TRUE, verbose=TRUE)
	# print(head(sssign))
}
ssin = list("data" = ssdata, "annot" = ssannot, "x" = sssign)
if (sig_save) {
	st = Sys.time()
	cat("Saving sig.score input and output to ", res_dir, "\n")
	saveRDS(ssin, file = paste(res_dir, "sig_score_inputs.RDS", sep = ""))
	saveRDS(sig_score, file = paste(res_dir, "sig_score_outputs.RDS", sep = ""))
	print(Sys.time()-st)
}


#################################################################################
cat("Generate subtype signature score and survival data\n")
sc_res = as.data.frame(sig_score$score)
colnames(sc_res) = "sig_score"
sc_res$pid = rownames(sc_res)
sub_scres = sc_res[sc_res$pid %in% sub_clin$pid,]
sub_sigscore = sub_scres

if (gp_gene != "") {
	cat("Using ", gp_gene, " as the group criteria...\n")
	sub_expr = singene_expr(gene = gp_gene, expr = expr, annot = annot, subdf = sub_scres, map = selfmap)
	sub_scres = sub_expr
	sub_scres[,"gpvalue"] = sub_scres[,gp_gene]
	hist_xlab = gp_gene
} else {
	sub_scres[,"gpvalue"] = sub_scres[,"sig_score"]
	hist_xlab = "Signature Score"
}

#################################################################################
## Histogram for sig.score
if (FALSE) {
	cat("Generate histogram plot of signature score\n")
	title = paste(db_name, sg_name, sep = "  ")
	sc_hist = ggplot(sub_scres, aes(x=gpvalue)) + 
		geom_histogram(color="darkblue", fill="lightblue", binwidth = 0.02) +
		labs(title=title, x=hist_xlab, y = "Count") + 
		theme_classic()
	sc_hist_bld = ggplot_build(sc_hist)
	sc_hist_data = sc_hist_bld$data[[1]]
	zero_x = sc_hist_data[sc_hist_data$count == 0,"x"]
	right_zero_count = sum(sc_hist_data[sc_hist_data$x>=zero_x[1],"count"])
	# print(right_zero_count)
	sc_hist = sc_hist + geom_vline(xintercept = zero_x[1], size = 1, colour = "grey",linetype = "dashed")
	if (gp_app == "oneqcut") {
		cat("Group the patient by one quantile cutoff: ", qcut,"\n")
		qcov = quantile(sub_scres$gpvalue, c(1-qcut)) # quantile cutoff value
		# Add line for cutoff value
		sc_hist = sc_hist + 
			geom_vline(xintercept = qcov[[1]], size = 1, colour = "purple",linetype = "dotdash")
		sc_hist = sc_hist + annotate("text", label = paste("Single cutoff value:\n", format(qcov[[1]],digit = 3), "(", 1-qcut, ")"),
					     hjust = 0, x = qcov[[1]], y = max(sc_hist_data$count)*0.96, size = 4.5, colour = "black")
		hist_tif = paste(res_dir, db_name, sg_name, "single_quantile_cutoff.tiff", sep = "_")

	}
	if (gp_app == "symqcut") {
		cat("Group the patient by one quantile cutoff with symmetric manner: ", qcut,"\n")

		qcov = quantile(sub_scres$gpvalue, c(qcut, 1-qcut)) # quantile cutoff value

		# Add line for cutoff value
		sc_hist = sc_hist + 
			geom_vline(xintercept = qcov[[1]], size = 1, colour = "purple",linetype = "dotdash")
		sc_hist = sc_hist + annotate("text", label = paste("Left cutoff value:\n", format(qcov[[1]],digit = 3), "(", qcut, ")"),
					     hjust = 0, x = qcov[[1]], y = max(sc_hist_data$count)*0.96, size = 4.5, colour = "black")
		sc_hist = sc_hist + 
			geom_vline(xintercept = qcov[[2]], size = 1, colour = "purple",linetype = "dotdash")
		sc_hist = sc_hist + annotate("text", label = paste("Right cutoff value:\n", format(qcov[[2]],digit = 3), "(", 1-qcut, ")"),
					     hjust = 0, x = qcov[[2]], y = max(sc_hist_data$count)*0.96, size = 4.5, colour = "black")
		hist_tif = paste(res_dir, db_name, sg_name, "symmetric_quantile_cutoff.tiff", sep = "_")
	}

	sc_hist = sc_hist + annotate("text", label = paste("Right side counts:", right_zero_count[1], 
							   "\n Percentile:", format(right_zero_count[1]/dim(sub_scres)[1]*100, digit = 4)), 
			 x = zero_x[1], y = max(sc_hist_data$count), size = 4.5, colour = "black")
	ggsave(sc_hist, file = hist_tif, width = 9, height = 6, units = "in", device = "tiff")
}
#################################################################################
## Assign survival data
cat("Extract survival and treatment information from clinical data to subtype sig.score data frame\n")
sub_scres$ost = 0
sub_scres$ose = 0


sub_clin = sub_clin[sub_clin$pid %in% sub_scres$pid,]
# print(head(sub_clin))
# print(dim(sub_clin))
sub_clin = sub_clin[sub_clin$ost != "--",]
# For overall survival
sub_scres$ost = sub_clin$ost[match(sub_scres$pid, sub_clin$pid)]
sub_scres$ose = sub_clin$ose[match(sub_scres$pid, sub_clin$pid)]
sub_scres$ost = as.numeric(as.character(sub_scres$ost))
sub_scres$ose = as.numeric(as.character(sub_scres$ose))

# For disease-free survival

# For treatment type



sub_scres = sub_scres[!is.na(sub_scres$ost),]
## Add all the other factors


# print(head(sub_scres))
print(dim(sub_scres))
# print(dim(sub_clin))
# q(save = "no")

#################################################################################
# Add correlation gene expressions SEPERATIVELY
if(corr_gene != "") { cat("Extract gene expression from expression data to subtype sig.score data frame for correlation...\n")
	sub_corr = sub_sigscore
	for (ig in 1:length(corr_gene)) {
		sub_corr = singene_expr(gene = corr_gene[ig], expr = expr, annot = annot, subdf = sub_corr, map = selfmap)
	}

	cat("Run correlation analysis between sig.score and other genes\n")
	sub_corr$pid = NULL
	sub_cor = as.matrix(sub_corr)
	subcorres = rcorr(sub_cor)
	# cor_plot = paste(, "/", "tex_singene_ER_cor.jpeg", sep = "") # change the folder name according to subtype and sig.score
	corr_tif = paste(res_dir, db_name, sg_name, "gene_correlation_v1.tiff", sep = "_")
	tiff(corr_tif, res = 180, width = 9, heigh = 6, units = "in") # save corrplot jpeg file name
	corrplot(subcorres$r, type="upper",
	         p.mat = subcorres$P, sig.level = 0.0001) ## Specialized the insignificant value according to the significant level
	dev.off()

    for (ig in 1:length(corr_gene)) { cat("Generate scatter plot of correlation between sig.score and", corr_gene[ig], "\n")
		scat_cor <- ggscatter(sub_corr, x = "sig_score", y = corr_gene[ig], # genes correlate with sig.score
		    color = "black", fill="darkgray",shape = 21, size = 2, #xlim=c(3,5),
	        add = "reg.line",  # Add regressin line
	        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
	        conf.int = TRUE, # Add confidence interval
	        cor.coef = TRUE, # Add correlation coefficient.
	        cor.coeff.args = list(method = "pearson", size = 4.5, label.sep = "\n") #label.x = 3,
		   )
		g <- scat_cor + theme_classic()+theme(axis.text=element_text(size=14, color = "#000000"),
                           axis.title=element_text(size=14))+xlab("Tex sig.score")
		scat_tif = paste(res_dir, db_name, sg_name,corr_gene[ig], "corr_scat.tiff", sep = "_")
		ggsave(g, file = scat_tif)
	}

}

if (cox_reg) {
	cat("\tSingle-variant Cox analysis with ", gptype, "\n")
	res.cox <- coxph(Surv(ost, ose) ~ gpvalue, data=sub_scres)
	prescox = summary(res.cox) # Print results of COX
	sum_txt = paste(res_dir, gptype, "_summary_cox_results.txt", sep = "")
	sink(sum_txt)
	print(summary(res.cox))
	sink()
	cox_rds <- paste(res_dir, gptype, "_survana_cox_res.rds", sep="")
	cat("\t\tP-value from logRank test: ", (prescox$sctest["pvalue"]), "\n") 
	saveRDS(res.cox, file=cox_rds)
}

if (mutation_corr) {
	mut = read.csv(paste(mutation_file, "laml_sampleSummary.csv", sep = ""))
	mut_corr = mut[,c("Tumor_Sample_Barcode", "total")]
	mut_corr$CD8A = sub_scres$CD8A[match(mut_corr$Tumor_Sample_Barcode, sub_scres$pid)]
	mut_corr$Tex = sub_scres$sig_score[match(mut_corr$Tumor_Sample_Barcode, sub_scres$pid)]
	mut_corr = mut_corr[!is.na(mut_corr$CD8A),]
	g <- ggscatter(mut_corr, x = "Tex", y = "total", # genes correlate with sig.score
			    color = "black", fill="darkgray",shape = 21, size = 2, #xlim=c(3,5),
		        add = "reg.line",  # Add regressin line
		        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
		        conf.int = TRUE, # Add confidence interval
		        cor.coef = TRUE, # Add correlation coefficient.
		        cor.coeff.args = list(method = "pearson", size = 4.5, label.sep = "\n") #label.x = 3,
			   )
	g = g +  scale_y_continuous(trans = log2_trans(),
                               breaks = trans_breaks("log2", function(x) 2^x),
                               labels = trans_format("log2", math_format(2^.x)))
	g = g + theme_classic()+theme(axis.text=element_text(size=14, color = "#000000"),
                           axis.title=element_text(size=14))+xlab("Tex sig.score")+ylab("Total mutation burden")
	ggsave(g, file = paste(res_dir, "tex_mutation_corr_scat.tiff", sep = ""), width = 5, height= 5)

	g <- ggscatter(mut_corr, x = "CD8A", y = "total", # genes correlate with sig.score
			    color = "black", fill="darkgray",shape = 21, size = 2, #xlim=c(3,5),
		        add = "reg.line",  # Add regressin line
		        add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
		        conf.int = TRUE, # Add confidence interval
		        cor.coef = TRUE, # Add correlation coefficient.
		        cor.coeff.args = list(method = "pearson", size = 4.5, label.sep = "\n") #label.x = 3,
			   )
	g = g + scale_y_continuous(trans = log2_trans(),
                               breaks = trans_breaks("log2", function(x) 2^x),
                               labels = trans_format("log2", math_format(2^.x)))
	g = g + theme_classic()+theme(axis.text=element_text(size=14, color = "#000000"),
                           axis.title=element_text(size=14))+xlab("CD8A")+ylab("Total mutation burden")
	ggsave(g, file = paste(res_dir, "cd8_mutation_corr_scat.tiff", sep = ""), width = 5, height = 5)


}
stop()
##########################################################################
## Assign groups

if (optimal_cutoff) {
	#optimal cutoff for group gene
	res.cox <- coxph(Surv(ost, ose) ~ gpvalue, data=sub_scres)
	res.cox <- cutp(res.cox)$gpvalue
	data.table::setorder(res.cox, "gpvalue")
	cutp_res <-(res.cox)[,"U"]
	colnames(cutp_res)[1] = "test_score"
	cutp_res$cutoff <- as.numeric(rownames(cutp_res))
	g <- ggplot(cutp_res, aes(cutoff, test_score)) + geom_col() + labs(y = "Log-rank Test Score") + theme_classic()+theme(axis.text=element_text(size=14, color = "#000000"),
                           axis.title=element_text(size=14))
	g <- g + geom_vline(xintercept = cutp_res$cutoff[cutp_res$test_score == max(cutp_res$test_score)], size = 0.5, colour = "red",linetype = "dotdash")
	g = g +  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
	 
	CD8_cutoff = res.cox$gpvalue[res.cox$U == max(res.cox$U)]
	if (length(CD8_cutoff) >1 ) {CD8_cutoff = CD8_cutoff[1]}
	print(paste(gp_gene, " max cutoff: " , CD8_cutoff, sep = ""))

	ggsave(g, file = paste(res_dir, gp_gene, "_log_rank_cutp_plot.tiff", sep = ""),  dpi = 300, width = 9, height = 6, units = "in", device = "tiff")
	
	sub_scres$group = "High"
	sub_scres$group[sub_scres$CD8A < CD8_cutoff] = "Low"

	survana(data = sub_scres, plot = res_dir, gptype = "CD8A")
	
	#optimal cutoff for signature
	sub_scres$gpvalue = sub_scres$sig_score
	res.cox <- coxph(Surv(ost, ose) ~ gpvalue, data=sub_scres)
	res.cox <- cutp(res.cox)$gpvalue
	data.table::setorder(res.cox, "gpvalue")
	cutp_res <-(res.cox)[,"U"]
	colnames(cutp_res)[1] = "test_score"
	cutp_res$cutoff <- as.numeric(rownames(cutp_res))
	c <- ggplot(cutp_res, aes(cutoff, test_score)) + geom_col() + labs(y = "Log-rank Test Score") + theme_classic()+theme(axis.text=element_text(size=14, color = "#000000"),
                           axis.title=element_text(size=14))
	c <- c + geom_vline(xintercept = cutp_res$cutoff[cutp_res$test_score == max(cutp_res$test_score)], size = 0.5, colour = "red",linetype = "dotdash")
	c = c + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
	Tex_cutoff = res.cox$gpvalue[res.cox$U == max(res.cox$U)]

	print(paste("Tex max cutoff: ", Tex_cutoff, sep = ""))

	ggsave(c, file = paste(res_dir, "Tex_log_rank_cutp_plot.tiff", sep = ""),  dpi = 300, width = 9, height = 6, units = "in", device = "tiff")
	
	if (length(Tex_cutoff) >1 ) {Tex_cutoff = Tex_cutoff[1]}

	sub_scres$CD8_group = sub_scres$group

	sub_scres$group = "High"
	sub_scres$group[sub_scres$sig_score < Tex_cutoff] = "Low"

	survana(data = sub_scres, plot = res_dir, gptype = "Tex sig.score")

	#combine correlation plot with log-rank plots
	sub_scres$CD8A=rescale(sub_scres$CD8A, to = c(0, 10))
	sub_scres$sig_score=rescale(sub_scres$sig_score, to = c(0, 10))
	corr_plot<-ggscatter(sub_scres, x = "sig_score", y = "CD8A", # genes correlate with sig.score
                     color = "gray40", shape = 19, size = 1.5, xlim = c(0,10), ylim = c(0,10),
                     add = "reg.line",  # Add regressin line
                     add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                     conf.int = TRUE, # Add confidence interval
                     cor.coef = TRUE, # Add correlation coefficient.
                     cor.coeff.args = list(method = "pearson", label.x = 0, label.sep = "\n")) + border() + theme(plot.margin = margin(0, 0, 0, 0, "pt"))
	corr_plot = corr_plot +  theme_classic()+theme(axis.text=element_text(size=14, color = "#000000"),
                           axis.title=element_text(size=14))+xlab("Tex sig.score")


	tex_plot<- c + clean_theme() +labs(title = NULL)
	gene_plot <- g + clean_theme () + labs(title = NULL) + rotate()
	combine_plot<-ggarrange(tex_plot, NULL, corr_plot, gene_plot, 
	                        ncol = 2, nrow = 2,  align = "hv", 
	                        widths = c(2, 1), heights = c(1, 2),
	                        common.legend = TRUE)
	combine_tiff <- paste(res_dir, "CD8A_Tex_combine_plot.tiff", sep = "_")
	ggsave(combine_plot, file = combine_tiff, dpi = 300, width = 11, height = 10, units = "in", device = "tiff" )

	##survival for four group
	sub_scres$tex_group = sub_scres$group
	sub_scres$group[sub_scres$CD8_group=="High" & sub_scres$tex_group == "High"]="CD8hiTexhi"
	sub_scres$group[sub_scres$CD8_group=="High" & sub_scres$tex_group == "Low"]="CD8hiTexlo"
	sub_scres$group[sub_scres$CD8_group=="Low" & sub_scres$tex_group == "Low"]="CD8loTexlo"
	sub_scres$group[sub_scres$CD8_group=="Low" & sub_scres$tex_group == "High"]="CD8loTexhi"
	sub_scres$group<-factor(sub_scres$group, levels = c("CD8hiTexhi", "CD8hiTexlo", "CD8loTexhi", "CD8loTexlo"))


	surv_csv = paste(res_dir, "optimal_cutoff_survival_analysis_inputs.csv")
	write.csv(sub_scres, file = surv_csv)

	numhh <- sum(sub_scres$CD8_group =="High" & sub_scres$tex_group == "High")
	numhl <- sum(sub_scres$CD8_group =="High" & sub_scres$tex_group == "Low")
	numlh <- sum(sub_scres$CD8_group =="Low" & sub_scres$tex_group == "High")
	numll <- sum(sub_scres$CD8_group =="Low" & sub_scres$tex_group == "Low")


	fit <- survfit(Surv(ost, ose) ~ CD8_group+tex_group, data=sub_scres)
	max_ost = max(sub_scres$ost)
	xlim_max = ceiling(max_ost/1000)*1000
	survcurv <- ggsurvplot(fit, data = sub_scres,
					       xlim = c(0,xlim_max), ylim = c(0.00, 1.00),
					       pval = TRUE, pval.size = 6, pval.coord = c(xlim_max/4*3, 0.95),
					       conf.int = FALSE, conf.int.alpha = 0.2, #legend = "none",
					       xlab = "Time (days)", ylab = "Overall Survival", legend.title = "",
					       legend.labs = c(paste("CD8hiTexhi, n =",numhh), paste("CD8hiTexlo, n =",numhl) ,paste("CD8loTexhi, n =",numlh), paste("CD8loTexlo, n =",numll)),
	#                                      surv.median.line = "hv",
					       ggtheme = theme_classic(),
					       palette = c("#D15466", "#0A9CC7","#D15466","#0A9CC7"), linetype = c(1,1,2,2), 
					       font.x = c(14, "bold"), font.y = c(18, "bold"), 
					       font.tickslab = c(16, "plain", "black"), font.legend = c(18, "bold"),
					       risk.table = FALSE, ncensor.plot = FALSE, legend = c(0.15, 0.14))
	survcurv$plot = survcurv$plot + 
    	theme(legend.text = element_text(size = 14, color = "black", face = "bold"))
	surv_plot <- paste(res_dir,"os_four_survana_group.tiff", sep="")
	ggsave(surv_plot, plot = survcurv$plot, dpi = 300, width = 9, height = 6, units = 'in')

	
} else {
	if (gp_app == "oneqcut") {
		cat("Group the patient by one quantile cutoff: ", qcut,"\n")
		qcov = quantile(sub_scres$gpvalue, c(1-qcut)) # quantile cutoff value
	}		
	if (gp_app == "symqcut") {
		cat("Group the patient by one quantile cutoff with symmetric manner: ", qcut,"\n")
		qcov = quantile(sub_scres$gpvalue, c(qcut, 1-qcut)) # qua
	}
	if (length(qcov) == 1) {
		sub_scres$group = "Medium"
		sub_scres[sub_scres$gpvalue <= qcov[[1]],"group"] = "Low"  
		sub_scres[sub_scres$gpvalue > qcov[[1]],"group"] = "High"
		surv_csv = paste(res_dir, "single_qcut_survival_analysis_inputs.csv")
		write.csv(sub_scres, file = surv_csv)
	}
	if (length(qcov) == 2) {
		sub_scres$group = "Medium"
		sub_scres[sub_scres$gpvalue <= qcov[[1]],"group"] = "Low"
		sub_scres[sub_scres$gpvalue >= qcov[[2]],"group"] = "High"
		surv_csv = paste(res_dir, "symmetric_qcut_survival_analysis_inputs.csv")
		write.csv(sub_scres, file = surv_csv)
		sub_scres = sub_scres[sub_scres$group != "Medium",]
		}
	if (length(qcov) > 2) {stop("Mulitple cutoffs!!!")}

	survana(data = sub_scres, plot = res_dir, gptype = gptype)
}
if (group_in_priMarker) {
	high = subset(sub_scres, group == "High")
	low = subset(sub_scres, group == "Low")

	high$CD8_group = high$group
	high$gpvalue = high$sig_score
	qcov = quantile(high$sig_score, c(1-qcut)) 
	high$group = "Medium"
	high[high$gpvalue <= qcov[[1]],"group"] = "Low"  
	high[high$gpvalue > qcov[[1]],"group"] = "High"

	survana(data = high, plot = res_dir, gptype = "CD8_high")

	low$CD8_group = low$group
	low$gpvalue = low$sig_score
	qcov = quantile(low$sig_score, c(1-qcut)) 
	low$group = "Medium"
	low[low$gpvalue <= qcov[[1]],"group"] = "Low"  
	low[low$gpvalue > qcov[[1]],"group"] = "High"
	survana(data = low, plot = res_dir, gptype = "CD8_low")

	sub_scres = rbind(high, low)

	sub_scres$tex_group<-sub_scres$group
	sub_scres$group[sub_scres$CD8_group=="High" & sub_scres$tex_group == "High"]="CD8hiTexhi"
	sub_scres$group[sub_scres$CD8_group=="High" & sub_scres$tex_group == "Low"]="CD8hiTexlo"
	sub_scres$group[sub_scres$CD8_group=="Low" & sub_scres$tex_group == "Low"]="CD8loTexlo"
	sub_scres$group[sub_scres$CD8_group=="Low" & sub_scres$tex_group == "High"]="CD8loTexhi"
	sub_scres$group<-factor(sub_scres$group, levels = c("CD8hiTexhi", "CD8hiTexlo", "CD8loTexhi", "CD8loTexlo"))

	numhh <- sum(sub_scres$CD8_group =="High" & sub_scres$tex_group == "High")
	numhl <- sum(sub_scres$CD8_group =="High" & sub_scres$tex_group == "Low")
	numlh <- sum(sub_scres$CD8_group =="Low" & sub_scres$tex_group == "High")
	numll <- sum(sub_scres$CD8_group =="Low" & sub_scres$tex_group == "Low")


	fit <- survfit(Surv(ost, ose) ~ CD8_group+tex_group, data=sub_scres)
	max_ost = max(sub_scres$ost)
	xlim_max = ceiling(max_ost/1000)*1000
	survcurv <- ggsurvplot(fit, data = sub_scres,
					       xlim = c(0,xlim_max), ylim = c(0.00, 1.00),
					       pval = TRUE, pval.size = 6, pval.coord = c(xlim_max/4*3, 0.95),
					       conf.int = FALSE, conf.int.alpha = 0.2, #legend = "none",
					       xlab = "Time (days)", ylab = "Overall Survival", legend.title = "",
					       legend.labs = c(paste("CD8hiTexhi, n =",numhh), paste("CD8hiTexlo, n =",numhl) ,paste("CD8loTexhi, n =",numlh), paste("CD8loTexlo, n =",numll)),
	#                                      surv.median.line = "hv",
					       ggtheme = theme_classic(),
					       palette = c("#D15466", "#0A9CC7","#D15466","#0A9CC7"), linetype = c(1,1,2,2), 
					       font.x = c(14, "bold"), font.y = c(18, "bold"), 
					       font.tickslab = c(16, "plain", "black"), font.legend = c(18, "bold"),
					       risk.table = FALSE, ncensor.plot = FALSE, legend = c(0.2, 0.17))
	surv_plot <- paste(res_dir,"os_four_survana_group_within_cd8.tiff", sep="")
	ggsave(surv_plot, plot = survcurv$plot, dpi = 300, width = 9, height = 6, units = 'in')

}




