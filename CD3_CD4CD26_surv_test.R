# This script is for survaval analysis according to any customized signature score of subtype samples
# Weihua Guo, Ph.D.
# 7/17/2019
# Jiayi Tan
# 07/01/2019

#rm(list = ls(all.names = TRUE))

survana <- function(data, type, gptype = "Sig.score", plot = "", csv = "", 
		    cox = "", coxfac = c("age","grade","tsize","node_stat"), multicox = TRUE, gdoi = 0) {
	# NOTE: plot/csv/cox should be a folder directory where to store the results
	# NOTE: gptype should be the definition or basis of the grouping method, default is signature score.

	if(type == "os") {
		cat("Analyzing overall survival...\n")
		surv_data = data[,c("group" ,"gpvalue" ,"ost", "ose", coxfac)]
		lab = "Overall survival"
	}
	if(type == "rfs") {
		cat("Analyzing relapse-free survival...\n")
		surv_data = data[,c("group" ,"gpvalue" ,"rfst" ,"rfse" , coxfac)]
		lab = "Relapse-free survival"
	}
	if(type == "dmfs") {
		cat("Analyzing distal metastasis-free survival...\n")
		surv_data = data[,c("group" ,"gpvalue" ,"dmfs" ,"dmfe" , coxfac)]
		lab = "Distal Metastasis-free survival"
	}

	colnames(surv_data)[3:4] <- c("time", "status")
	surv_data$status = as.numeric(surv_data$status)
	numh = sum(surv_data$group == "High")
	numl = sum(surv_data$group == "Low")
	if (dim(surv_data)[1] != numh + numl) {stop("More than two groups!!!")}
#	print(head(surv_data))
	fit <- survfit(Surv(time, status) ~ group, data=surv_data)
	fit_df <- data.frame(time = fit$time, n.risk = fit$n.risk,
			     n.event = fit$n.event, n.censor = fit$n.censor,
			     surv = fit$surv, upper = fit$upper, lower = fit$lower)
	print(head(fit_df))
	if (cox != "") {
		#if (sum(surv_data$status == 1)>0) {
		if (multicox) {
			cat("\tMulti-variant Cox analysis with ", coxfac, "\n")
			if (gdoi == 0) {
				if(db_name != "tcga_brca") {
					res.cox <- coxph(Surv(time, status) ~ gpvalue+age+grade+tsize+node_stat, data=surv_data)
				} else {
					res.cox <- coxph(Surv(time, status) ~ gpvalue+age+grade+node_stat, data=surv_data)
				}	
			} else {
				if(db_name != "tcga_brca") {
					res.cox <- coxph(Surv(time, status) ~ gpvalue+age+tsize+node_stat, data=surv_data)
					} else {
						res.cox <- coxph(Surv(time, status) ~ gpvalue+age+node_stat, data=surv_data)
					} 
			}
		} else {
				cat("\tSingle-variant Cox analysis with ", gptype, "\n")
				res.cox <- coxph(Surv(time, status) ~ gpvalue, data=surv_data)
	    }
		prescox = summary(res.cox) # Print results of COX
		sum_txt = paste(cox, lab, gptype, "_summary_cox_results.txt", sep = "")
		sink(sum_txt)
		print(summary(res.cox))
		sink()
		cox_rds <- paste(cox, lab, gptype, "_survana_cox_res.rds", sep="")
		cat("\t\tP-value from logRank test: ", (prescox$sctest["pvalue"]), "\n") 
		saveRDS(res.cox, file=cox_rds)
		#} else {cat("WARNING: NO event happened in this group!\n")}
	}

	if (csv != "") {
		surv_rescsv <- paste(csv, lab, gptype, "_survana_res.csv", sep="")
		write.csv(fit_df, file=surv_rescsv)}
	if (plot != "") {
		surv_plot <- paste(plot, lab, gptype, "_survana_curves.tiff", sep="")
		max_st = max(surv_data$time)
		xlim_max = ceiling(max_st/1000)*1000
		survcurv <- ggsurvplot(fit, data = surv_data,
				       #xlim = c(0,11000),  hrtype is ER+
				       xlim = c(0,xlim_max), ylim = c(0.00, 1.00),  
				       pval = TRUE, pval.size = 6, pval.coord = c(xlim_max/4*3, 0.1),
				       conf.int = FALSE, conf.int.alpha = 0.2, #legend = "none",
				       xlab = "Time (days)", ylab = lab, legend.title = "",
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
# Personally prefer to having a independent section for all the parameters or variables
# which may be tuned for different inputs and tasks

# For each variable, please at least leave a simple comments describing what this 
# variable represent for

#################################################################################
# work_dir = "/home/weihua/mnts/group_plee/Weihua/surv_validation/" # working directory/path for survival validation
work_dir = "//Bri-net/citi/Peter Lee Group/Weihua/surv_validation/"
db_name = "metabric"
#db_name = "tcga_brca"
#sg_name = "loi_trm" # Loi's TRM sig
sg_name = "CD26_update_0417" # Colt's Tex sig from breast tissue c2
#sg_name = "mamma" # mamma sig
expr_type = "median" # ilid: raw data from EGA, median: raw median data from cbioportal, medianz: zscore from cbioportal
selfmap = FALSE # NOTE: ilid/tcga require this as TRUE; median as FALSE

# data_dir = "/home/weihua/mnts/group_plee/Weihua/metabric_use/" # directory/path for public data
data_dir = paste(work_dir, db_name, "/", sep = "") # generate the directory with all the public data
sign_dir = paste(work_dir, sg_name, "/", sep = "") # generate the directory with signatures and corresponding results


if (db_name == "metabric") {
	if (expr_type == "ilid") {expr_file = "metabric_expr_ilid.RDS"}
	if (expr_type == "median") {expr_file = "data_expression_median.RDS"}
	if (expr_type == "medianz") {expr_file = "data_mRNA_median_Zscores.RDS"}
	clin_rds = "merge_clin_info_corrected.RDS" # clinical information with merged disease-free survival
	annot_file = "HumanHT_12_v30_R3_cleaned_v2.xlsx" # Microarray/Genome annotation
}
if (db_name == "tcga_brca") {
#expr_file = "tcga_brca_log2trans_fpkm_uq_v2.RDS" # Expression file
#expr_file = "breast_Primary_Tumor_htseq.counts_merged_tcga_v2.RDS"
expr_file = "tcga_brca_log2trans_fpkm_uq_v2.RDS" # Expression file

clin_rds = "07212019_tcga_clinical_info.RDS" # clinical information with merged disease-free survival 
annot_file = "gencode.gene.info.v22.xlsx" # Microarray/Genome annotation
}

#sign_file = "loi_trm_signature.txt" # Signature file Loi's TRM
#sign_file = "tex_signature_colt_c2.txt" # Signature file Colt's Tex
#sign_file = "cluster0_t_pca9_15_ER_raw_v1.txt"
sign_file = "cd26_c1_all_tumor_v1.txt"
#sign_file = "tumor_ln_imputated_tl_v1.txt"
#sign_file = "tumor_ln_normal_pbmc_abc_v1.txt"
#sign_file = "tumor_ln_normal_pbmc_raw_v1.txt"
#sign_file = "tumor_ln_unimputated_tlni_v1.txt"
#sign_file = "tumor_obc_v1.txt"
#sign_file = "mamma_signature_v1.txt" # Signature file mamma


histype = "IDC" # histology type: IDC/DCIS
pamst = ""#c("LumA", "LumB") # PAM50 status: LumA/LumB/Basal/Normal/Her2
gdoi = 0 #c(1) # Grade of interest: 1/2/3
hrtype = c("N", "N", "N") # N: Negative, P: Positive, "-": DON'T CARE
sig_save = FALSE
gp_app = "symqcut"#"symqcut" # oneqcut: one quantile cutoff (upper percential), symqcut: symmetric quantile cutoff
qcut = 0.5 #0.25 # This is TOP quantile for oneqcut approach
gp_gene = "CD3D" # Group gene used for categorizing the cohort(if run cox regression of single gene)
# Default "": use signature score AE", "STAT1") # Genes need to be correlated with signature scores
gptype = ""
corr_gene = "" #c("CD8A", "CD3G", "ITG
trt_type = "" #c("ct", "rt", "ht") # check the correlation between sig.score and treatment
strata = 2
group_in_priMarker = TRUE



#################################################################################
# Work for experiment records
res_folder = "sym50_CD4CD26_withinCD3D_TNBC_cbpt" # NOTE: Please change this folder name to identify your experiments
res_dir = paste(sign_dir, res_folder, "/", sep ="")
dir.create(file.path(sign_dir, res_folder), showWarnings = FALSE)
# COPY the used script to the result folder for recording what experiment was run
### !!!Please change the script_dir to the folder directory where this script is located
script_dir = "~/GitHub/surv_analysis/"
script_name = "CD3_CD4CD26_surv_test.R"
file.copy(paste(script_dir, script_name, sep = ""), res_dir)  

#################################################################################
cat("Loading expression data...\n")
st = Sys.time()
## Please use either the full path of the file or change the work directory here
expr = readRDS(paste(data_dir, expr_file, sep = ""))
#expr = readRDS("metabric_expr_ilid.RDS") # When test the script using metabric
#expr = readRDS("tcga_brca_log2trans_fpkm_uq_v2.RDS") # When test the script using tcga
#expr = readRDS("data_expression_median.RDS") # When test the script using cBioportal
#expr = readRDS("tcga_portal_data_expr_v2.RDS")
print(Sys.time()-st)
# print(meta_expr[1:9,1:6]) # Check the input in terminal
expr<-as.data.frame(expr)

cat("Loading clinical data...\n")
st = Sys.time()
# print(head(clin_info))
## To merge the RFS
if (FALSE) {
	clin_file = "merge_clin_info.xlsx" # clinical information
	clin_info = as.data.frame(read_excel(paste(data_dir, clin_file, sep = ""),1))
	cat("Generate the relapse-free survival...\n")
	cat("NOTE: the script will be quited after this...\n")
	clin_info$TOR = ifelse(clin_info$TLR > clin_info$TDR, clin_info$TDR, clin_info$TLR)
	clin_info$OR = ifelse(clin_info$TLR > clin_info$TDR, clin_info$DR, clin_info$LR)
# 	print(clin_info[sample(dim(clin_info)[1], 9),]) # Random select 9 patients to check
	rds_file = paste(data_dir, clin_rds, sep = "")
	saveRDS(clin_info, file = rds_file)
	cat("Save the clinical information to ", rds_file, "\n")
	q(save = "no")
}
clin_info = readRDS(paste(data_dir, clin_rds, sep = ""))
#clin_info = read_excel("1018_tcga_pam50_clin_rsem.xlsx", sheet = 1)
#clin_info = readRDS("merge_clin_info_v3.RDS") # When test the script
#clin_info = read_excel("tcga_portal_clin_info_v2.xlsx", sheet= 1 )
#clin_info = read_excel("07212019_tcga_clinical_info_early.xlsx", sheet = 1)
#clin_info = read_excel("08272019_tcga_pam50_clin.xlsx", sheet = 1)
#clin_info = read_excel("ega_clin_info_mole_sub.xlsx", sheet = 1)
#clin_info = as.data.frame(read_excel(paste(data_dir, clin_rds, sep = "")))
# saveRDS(clin_info, file = paste(data_dir, "07212019_tcga_clinical_info.RDS", sep = ""))
print(Sys.time()-st)

cat("Start to filter by clinical info...\n")
sub_clin = clin_info
cat("\tOriginal patient number: ", dim(sub_clin)[1], "\n")
if (histype != "") {
	# For IDC
	sub_clin = sub_clin[sub_clin[,"oncotree_code"] == histype,]
	sub_clin = sub_clin[complete.cases(sub_clin$pid),]
	cat("\tFiltered patient number: ", dim(sub_clin)[1], "\n")
}


if (gdoi != 0) {
	# print(head(sub_clin))
	sub_clin = sub_clin[sub_clin[,"grade"] %in% gdoi,]
	sub_clin = sub_clin[complete.cases(sub_clin$pid),]
	cat("\tFiltered patient number: ", dim(sub_clin)[1], "\n")
}
if (pamst != "") {
	cat("Using PAM50 as molecular subtype classifier: ", pamst, "\n")
	#sub_clin = subset(sub_clin, Pam50Subtype == pamst)
	sub_clin = subset(sub_clin, pam50 %in% pamst)
	sub_clin = sub_clin[complete.cases(sub_clin$pid),]
	cat("\tFiltered patient number: ", dim(sub_clin)[1], "\n")
} else {
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
annot = as.data.frame(read_excel(paste(data_dir, annot_file, sep = ""))) 
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

	## RUN sig.score
	cat("Run sig.score in all patient samples\n")
	sig_score <- sig.score(x=sssign, data=ssdata, annot=ssannot, do.mapping=TRUE, signed=TRUE, verbose=TRUE)
	# print(head(sig_score$score))

} else {
	## PREP for sig.score WITHOUT MAPPPING
	# data = [r:sample, c:gene]
	# annot = [r:gene, c: genome annotation] 
	# x = [r: genes in signature, c: probe+EntrezGene.ID, coefficient]
	if (expr_type == "median") {rownames(expr) = expr[,1]}
	ssdata = expr[, 3:dim(expr)[2]]
	ssdata = as.data.frame(t(ssdata))
	cat("Expression data dimension, sample: probe ", dim(ssdata), "\n")
	print(ssdata[1:9,1:6])

	ssannot = expr[,1:2]
	cat("Annotation data dimension, probe: annot ", dim(ssannot), "\n")
	print(head(ssannot))

	sign_gene = rownames(sign)
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
cat("Generate histogram plot of signature score\n")
title = paste(db_name, sg_name, pamst, sep = "  ")
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
	hist_tif = paste(res_dir, db_name, sg_name, pamst, "single_quantile_cutoff.tiff", sep = "_")

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
	hist_tif = paste(res_dir, db_name, sg_name, pamst, "symmetric_quantile_cutoff.tiff", sep = "_")
}

sc_hist = sc_hist + annotate("text", label = paste("Right side counts:", right_zero_count[1], 
						   "\n Percentile:", format(right_zero_count[1]/dim(sub_scres)[1]*100, digit = 4)), 
		 x = zero_x[1], y = max(sc_hist_data$count), size = 4.5, colour = "black")
# hist_tif = paste(sign_dir, db_name, sg_name, pamst, ".tiff", sep = "_")
ggsave(sc_hist, file = hist_tif, width = 9, height = 6, units = "in", device = "tiff")

#################################################################################
## Assign survival data
cat("Extract survival and treatment information from clinical data to subtype sig.score data frame\n")
sub_scres$ost = 0
sub_scres$ose = 0
sub_scres$rfst = 0
sub_scres$rfse = 0
sub_scres$dmfs = 0
sub_scres$dmfe = 0
if(db_name != "tcga_brca") {
	sub_scres$CT = 0
	sub_scres$RT = 0
	sub_scres$HT = 0
	sub_scres$ct = "YES"
	sub_scres$rt = "YES"
	sub_scres$ht = "YES"
}


sub_clin = sub_clin[sub_clin$pid %in% sub_scres$pid,]
# print(head(sub_clin))
# print(dim(sub_clin))

# For overall survival
if (db_name == "metabric") {
	#sub_scres[sub_clin$pid,"ost"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"T"]
	#sub_scres[sub_clin$pid,"ose"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"DeathBreast"]
	#sub_scres[sub_clin$pid,"rfst"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"TOR"]
	#sub_scres[sub_clin$pid,"rfse"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"OR"]

	sub_scres$ost = sub_clin$T[match(sub_scres$pid, sub_clin$pid)]
	sub_scres$ose = sub_clin$DeathBreast[match(sub_scres$pid, sub_clin$pid)]
	sub_scres$rfst = sub_clin$TOR[match(sub_scres$pid, sub_clin$pid)]
	sub_scres$rfse = sub_clin$OR[match(sub_scres$pid, sub_clin$pid)]
	sub_scres$dmfs = sub_clin$DMFS[match(sub_scres$pid, sub_clin$pid)]
	sub_scres$dmfe = sub_clin$DMFE[match(sub_scres$pid, sub_clin$pid)]

}

if (db_name == "tcga_brca") {
	sub_scres$ost=sub_clin$T[match(sub_scres$pid, sub_clin$pid)]
	sub_scres$ose=sub_clin$DeathBreast[match(sub_scres$pid, sub_clin$pid)]

}
#sub_scres[sub_clin$pid,"ose"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"ose"]
sub_scres$ost <- as.numeric(sub_scres$ost)
sub_scres$ose <- as.numeric(sub_scres$ose)
sub_scres$dmfs <- as.numeric(sub_scres$dmfs)
sub_scres$dmfe <- as.numeric(sub_scres$dmfe)
# For disease-free survival


# For treatment type
if(FALSE) {
	if(db_name != "tcga_brca") {
		sub_scres[sub_clin$pid,"CT"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"CT"]
		sub_scres[sub_clin$pid,"RT"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"RT"]
		sub_scres[sub_clin$pid,"HT"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"HT"]
		sub_scres$ct[sub_scres$CT == "NO/NA"] = "NO"
		sub_scres$rt[sub_scres$RT == "NO/NA"] = "NO"
		sub_scres$ht[sub_scres$HT == "NO/NA"] = "NO"
	}
}




if(trt_type != "") { 
	for (it in 1:length(trt_type)) { cat("Generate sig.score box plot with", trt_type[it], "treatment\n")
		sig_trt <- ggboxplot(sub_scres, x = trt_type[it], y = "sig_score",
		               color = trt_type[it], palette = "jco",
		               add = "jitter") + stat_compare_means(method = "t.test")
        box_tif = paste(res_dir, db_name, sg_name, pamst, trt_type[it], "corr.tiff", sep = "_")
        ggsave(sig_trt, file = box_tif)
	}
}

## Add all the other factors
#sub_scres[sub_clin$pid,"age"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"age_at_diagnosis"]
#sub_scres[sub_clin$pid,"grade"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"grade"]
#sub_scres$grade=as.numeric(sub_scres$grade)
#sub_scres[sub_clin$pid,"tsize"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"tumor_size"]
#sub_scres$tsize=as.numeric(sub_scres$tsize)
#sub_scres[sub_clin$pid,"node_stat"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"Lymph.Nodes.Positive"]
#sub_scres$node_stat=as.numeric(sub_scres$node_stat)

sub_scres$age = sub_clin$age_at_diagnosis[match(sub_scres$pid, sub_clin$pid)]
sub_scres$grade = sub_clin$grade[match(sub_scres$pid, sub_clin$pid)]
sub_scres$tsize = sub_clin$tumor_size[match(sub_scres$pid, sub_clin$pid)]
sub_scres$node_stat = sub_clin$Lymph.Nodes.Positive[match(sub_scres$pid, sub_clin$pid)]
sub_scres$menopausal_State = sub_clin$Inferred.Menopausal.State[match(sub_scres$pid, sub_clin$pid)]

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
	corr_tif = paste(res_dir, db_name, sg_name, pamst, "gene_correlation_v1.tiff", sep = "_")
	tiff(corr_tif, res = 180, width = 9, heigh = 6, units = "in") # save corrplot jpeg file name
	corrplot(subcorres$r, type="upper",
	         p.mat = subcorres$P, sig.level = 0.0001) ## Specialized the insignificant value according to the significant level
	dev.off()

    for (ig in 1:length(corr_gene)) { cat("Generate scatter plot of correlation between sig.score and", corr_gene[ig], "\n")
		scat_cor<-ggscatter(sub_corr, x = "sig_score", y = corr_gene[ig], # genes correlate with sig.score
		   color = "black", shape = 21, size = 2, 
		   add = "reg.line",  # Add regressin line
		   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
		   conf.int = TRUE, # Add confidence interval
		   cor.coef = TRUE, # Add correlation coefficient.
		   cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
		   )
		scat_tif = paste(res_dir, db_name, sg_name, pamst, corr_gene[ig], "corr_scat.tiff", sep = "_")
		ggsave(scat_cor, file = scat_tif)
	}

}

#################################################################################
## Assign groups
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

survana(data = sub_scres, type = "os", plot = res_dir, gptype = gptype, 
	cox = res_dir, multicox = TRUE, coxfac = c("age","tsize","grade", "node_stat"), gdoi = gdoi)
if (db_name != "tcga_brca") {
	survana(data = sub_scres, type = "rfs", plot = res_dir, gptype = gptype, 
		cox = res_dir, coxfac = c("age","tsize", "grade", "node_stat"), gdoi = gdoi)
	survana(data = sub_scres, type = "dmfs", plot = res_dir, gptype = gptype, 
		cox = res_dir, coxfac = c("age","tsize", "grade", "node_stat"), gdoi = gdoi)}

if (strata ==2) {

	if (group_in_priMarker) {
		high = subset(sub_scres, group == "High")
		low = subset(sub_scres, group == "Low")

		high$CD3_group = high$group
		high$gpvalue = high$sig_score
		qcov = quantile(high$sig_score, c(1-qcut)) 
		high$group = "Medium"
		high[high$gpvalue <= qcov[[1]],"group"] = "Low"  
		high[high$gpvalue > qcov[[1]],"group"] = "High"

		survana(data = high, type = "os", plot = res_dir, gptype = "CD3high",
			cox = "", coxfac = c("age","tsize", "grade", "node_stat"), gdoi = gdoi)
		survana(data = high, type = "rfs", plot = res_dir, gptype = "CD3high",
			cox = "", coxfac = c("age","tsize", "grade", "node_stat"), gdoi = gdoi)
		survana(data = high, type = "dmfs", plot = res_dir, gptype = "CD3high",
			cox = "", coxfac = c("age","tsize", "grade", "node_stat"), gdoi = gdoi)

		low$CD3_group = low$group
		low$gpvalue = low$sig_score
		qcov = quantile(low$sig_score, c(1-qcut)) 
		low$group = "Medium"
		low[low$gpvalue <= qcov[[1]],"group"] = "Low"  
		low[low$gpvalue > qcov[[1]],"group"] = "High"
		survana(data = low, type = "os", plot = res_dir, gptype = "CD3low",
			cox = "", coxfac = c("age","tsize", "grade", "node_stat"), gdoi = gdoi)
		survana(data = low, type = "rfs", plot = res_dir, gptype = "CD3low",
			cox = "", coxfac = c("age","tsize", "grade", "node_stat"), gdoi = gdoi)
		survana(data = low, type = "dmfs", plot = res_dir, gptype = "CD3low",
			cox = "", coxfac = c("age","tsize", "grade", "node_stat"), gdoi = gdoi)
		
		sub_scres = rbind(high, low)

		sub_scres$CD26_group<-sub_scres$group
		sub_scres$group[sub_scres$CD3_group=="High" & sub_scres$CD26_group == "High"]="CD3hiCD26hi"
		sub_scres$group[sub_scres$CD3_group=="High" & sub_scres$CD26_group == "Low"]="CD3hiCD26lo"
		sub_scres$group[sub_scres$CD3_group=="Low" & sub_scres$CD26_group == "Low"]="CD3loCD26lo"
		sub_scres$group[sub_scres$CD3_group=="Low" & sub_scres$CD26_group == "High"]="CD3loCD26hi"
		sub_scres$group<-factor(sub_scres$group, levels = c("CD3hiCD26hi", "CD3hiCD26lo", "CD3loCD26hi", "CD3loCD26lo"))

		numhh <- sum(sub_scres$CD3_group =="High" & sub_scres$CD26_group == "High")
		numhl <- sum(sub_scres$CD3_group =="High" & sub_scres$CD26_group == "Low")
		numlh <- sum(sub_scres$CD3_group =="Low" & sub_scres$CD26_group == "High")
		numll <- sum(sub_scres$CD3_group =="Low" & sub_scres$CD26_group == "Low")


		fit <- survfit(Surv(ost, ose) ~ CD3_group+CD26_group, data=sub_scres)
		survcurv <- ggsurvplot(fit, data = sub_scres,
						       xlim = c(0,11000), ylim = c(0.00, 1.00),
						       pval = TRUE, pval.size = 6, pval.coord = c(7500, 0.1),
						       conf.int = FALSE, conf.int.alpha = 0.2, #legend = "none",
						       xlab = "Time (days)", ylab = "Overall Survival", legend.title = "",
						       legend.labs = c(paste("CD3hiCD26hi, n =",numhh), paste("CD3hiCD26lo, n =",numhl) ,paste("CD3loCD26hi, n =",numlh), paste("CD3loCD26lo, n =",numll)),
		#                                      surv.median.line = "hv",
						       ggtheme = theme_classic(),
						       palette = c("#D15466", "#0A9CC7","#D15466","#0A9CC7"), linetype = c(1,1,2,2), 
						       font.x = c(14, "bold"), font.y = c(18, "bold"), 
						       font.tickslab = c(16, "plain", "black"), font.legend = c(18, "bold"),
						       risk.table = FALSE, ncensor.plot = FALSE, legend = c(0.2, 0.17))
		surv_plot <- paste(res_dir,"os_four_survana_cd26group_within_CD3.tiff", sep="")
		ggsave(surv_plot, plot = survcurv$plot, dpi = 300, width = 9, height = 6, units = 'in')

		fit <- survfit(Surv(rfst, rfse) ~ CD3_group+CD26_group, data=sub_scres)
		survcurv <- ggsurvplot(fit, data = sub_scres,
						       xlim = c(0,11000), ylim = c(0.00, 1.00),
						       pval = TRUE, pval.size = 6, pval.coord = c(7500, 0.1),
						       conf.int = FALSE, conf.int.alpha = 0.2, #legend = "none",
						       xlab = "Time (days)", ylab = "Relapse-free Survival", legend.title = "",
						       legend.labs = c(paste("CD3hiCD26hi, n =",numhh), paste("CD3hiCD26lo, n =",numhl) ,paste("CD3loCD26hi, n =",numlh), paste("CD3loCD26lo, n =",numll)),
		#                                      surv.median.line = "hv",
						       ggtheme = theme_classic(),
						       palette = c("#D15466", "#0A9CC7","#D15466","#0A9CC7"), linetype = c(1,1,2,2), 
						       font.x = c(14, "bold"), font.y = c(18, "bold"), 
						       font.tickslab = c(16, "plain", "black"), font.legend = c(18, "bold"),
						       risk.table = FALSE, ncensor.plot = FALSE, legend = c(0.2, 0.17))
		surv_plot <- paste(res_dir,"rfs_four_survana_cd26group_within_CD3.tiff", sep="")
		ggsave(surv_plot, plot = survcurv$plot, dpi = 300, width = 9, height = 6, units = 'in')

		fit <- survfit(Surv(dmfs, dmfe) ~ CD3_group+CD26_group, data=sub_scres)
		survcurv <- ggsurvplot(fit, data = sub_scres,
						       xlim = c(0,11000), ylim = c(0.00, 1.00),
						       pval = TRUE, pval.size = 6, pval.coord = c(7500, 0.1),
						       conf.int = FALSE, conf.int.alpha = 0.2, #legend = "none",
						       xlab = "Time (days)", ylab = "Distal Metastasis-free Survival", legend.title = "",
						       legend.labs = c(paste("CD3hiCD26hi, n =",numhh), paste("CD3hiCD26lo, n =",numhl) ,paste("CD3loCD26hi, n =",numlh), paste("CD3loCD26lo, n =",numll)),
		#                                      surv.median.line = "hv",
						       ggtheme = theme_classic(),
						       palette = c("#D15466", "#0A9CC7","#D15466","#0A9CC7"), linetype = c(1,1,2,2), 
						       font.x = c(14, "bold"), font.y = c(18, "bold"), 
						       font.tickslab = c(16, "plain", "black"), font.legend = c(18, "bold"),
						       risk.table = FALSE, ncensor.plot = FALSE, legend = c(0.2, 0.17))
		surv_plot <- paste(res_dir,"dmfs_four_survana_group_within_CD3.tiff", sep="")
		ggsave(surv_plot, plot = survcurv$plot, dpi = 300, width = 9, height = 6, units = 'in')

	} else {
		sub_scres[,"CD3_group"] = sub_scres[,"group"]
		sub_scres[,"gpvalue"] = sub_scres[,"sig_score"]
		qcov = quantile(sub_scres$gpvalue, c(1-qcut))
		#qcov<-tex_qcov
		sub_scres$group = "Medium"
		sub_scres[sub_scres$gpvalue <= qcov[[1]],"group"] = "Low"  
		sub_scres[sub_scres$gpvalue > qcov[[1]],"group"] = "High"
		
		numhh <- sum(sub_scres$CD3_group =="High" & sub_scres$group == "High")
		numhl <- sum(sub_scres$CD3_group =="High" & sub_scres$group == "Low")
		numlh <- sum(sub_scres$CD3_group =="Low" & sub_scres$group == "High")
		numll <- sum(sub_scres$CD3_group =="Low" & sub_scres$group == "Low")

		fit <- survfit(Surv(ost, ose) ~ CD3_group+group, data=sub_scres)
		survcurv <- ggsurvplot(fit, data = sub_scres,
						       xlim = c(0,11000), ylim = c(0.00, 1.00),
						       pval = TRUE, pval.size = 6, pval.coord = c(8500, 0.12),
						       conf.int = FALSE, conf.int.alpha = 0.2, #legend = "none",
						       xlab = "Time (days)", ylab = "Overall Survival", legend.title = "",
						       legend.labs = c(paste("CD3hiCD26hi, n =",numhh), paste("CD3hiCD26lo, n =",numhl) ,paste("CD3loCD26hi, n =",numlh), paste("CD3loCD26lo, n =",numll)),
		#                                       surv.median.line = "hv",
						       ggtheme = theme_classic(),
						       palette = c("#D15466", "#0A9CC7","#D15466","#0A9CC7"), linetype = c(1,1,2,2), 
						       font.x = c(14, "bold"), font.y = c(18, "bold"), 
						       font.tickslab = c(16, "plain", "black"), font.legend = c(18, "bold"),
						       risk.table = FALSE, ncensor.plot = FALSE, legend = c(0.2, 0.17))
		surv_plot <- paste(res_dir,"os_four_survana.tiff", sep="")
		ggsave(surv_plot, plot = survcurv$plot, dpi = 300, width = 9, height = 6, units = 'in')
		

		fit <- survfit(Surv(rfst, rfse) ~ CD3_group+group, data=sub_scres)
		survcurv <- ggsurvplot(fit, data = sub_scres,
						       xlim = c(0,11000), ylim = c(0.00, 1.00),
						       pval = TRUE, pval.size = 6, pval.coord = c(8500, 0.12),
						       conf.int = FALSE, conf.int.alpha = 0.2, #legend = "none",
						       xlab = "Time (days)", ylab = "Relapse-free Survival", legend.title = "",
						       legend.labs = c(paste("CD3hiCD26hi, n =",numhh), paste("CD3hiCD26lo, n =",numhl) ,paste("CD3loCD26hi, n =",numlh), paste("CD3loCD26lo, n =",numll)),
		#                                      surv.median.line = "hv",
						       ggtheme = theme_classic(),
						       palette = c("#D15466", "#0A9CC7","#D15466","#0A9CC7"), linetype = c(1,1,2,2), 
						       font.x = c(14, "bold"), font.y = c(18, "bold"), 
						       font.tickslab = c(16, "plain", "black"), font.legend = c(18, "bold"),
						       risk.table = FALSE, ncensor.plot = FALSE, legend = c(0.2, 0.17))
		surv_plot <- paste(res_dir,"rfs_four_survana.tiff", sep="")
		ggsave(surv_plot, plot = survcurv$plot, dpi = 300, width = 9, height = 6, units = 'in')

	}
}
#################################################################################
