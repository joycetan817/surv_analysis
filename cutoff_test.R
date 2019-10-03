rm(list = ls(all.names = TRUE))

# Extract single gene expression from expression data for grouping samples
singene_expr = function (gene, expr, annot, subdf, caltype = "mean", map = TRUE) {
	if (map) {
		# caltype: calculation type for multiple probes, mean, median, max, min
		gene_info = subset(annot, ILMN_Gene == gene) # ILMN_Gene is the column of gene name which is used to extract gene probe ID
		if (dim(gene_info)[1] == 0) {stop("No MATCHED gene found!!!")}
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


cutoff_screen <- function (data, type, qcut, gp_app, gp_val) {
	cat("Generate survival/cox pvalue from all cutoffs\n")
	pval_cut<-as.data.frame(matrix(ncol=2, nrow=length(qcut)))
	colnames(pval_cut)[1:2]= c("pval", "log_score")
	rownames(pval_cut) = qcut

	for (cut in 1: length(qcut)) {
		temp_scres <- data
		temp_scres$gpvalue <- temp_scres[,gp_val]
		if (gp_app == "oneqcut") {
			qcov = quantile(temp_scres$gpvalue, c(1-qcut[cut])) # one quantile cutoff
			temp_scres$group = "Medium"
			temp_scres[temp_scres$gpvalue <= qcov[[1]],"group"] = "Low"  
			temp_scres[temp_scres$gpvalue > qcov[[1]],"group"] = "High"
		}
		if (gp_app == "symqcut") {
			qcov = quantile(temp_scres$gpvalue, c(qcut[cut], 1-qcut[cut])) # symmetric quantile cutoff
			temp_scres$group = "Medium"
			temp_scres[temp_scres$gpvalue <= qcov[[1]],"group"] = "Low"
			temp_scres[temp_scres$gpvalue >= qcov[[2]],"group"] = "High"
			temp_scres = temp_scres[temp_scres$group != "Medium",]
		}
		if (type == "surv") {
			fit <- survfit(Surv(ost, ose) ~ group, data = temp_scres)
			pval_cut[cut, 1]=-log10(surv_pvalue(fit)[2])
		} 
		if (type == "cox_reg") {
		res.cox <- coxph(Surv(ost, ose) ~ group, data=temp_scres)
		prescox = summary(res.cox)
		pval_cut[cut, 1]=-log10(prescox$sctest["pvalue"])
		pval_cut[cut, 2]=prescox$sctest["test"]
		}
	}

	if (type == "surv") {
		pval_cut[2] = NULL
	}
	pval_csv = paste(res_dir, gp_val, "cut_pval.csv")
	write.csv(pval_cut, file= pval_csv)
################################################################################
	pval_cut$cutoff <- as.numeric(rownames(pval_cut))
	if( type == "cox_reg") {
		c <- ggplot(pval_cut, aes(cutoff, log_score))
		c <- c + geom_col() + labs(title = paste(gp_val, db_name, pamst), x ="cutoff", y = "log-rank test score") + coord_cartesian(xlim=c(0, 1)) + theme_classic()
		c <- c + geom_vline(xintercept = pval_cut$cutoff[pval_cut$log_score == max(pval_cut$log_score)], size = 0.5, colour = "red",linetype = "dotdash")
		ggsave(c, file = paste(res_dir, gp_val, "_logscore_barplot.tiff"),  dpi = 300, width = 9, height = 6, units = "in", device = "tiff")
		
		d <- ggplot(pval_cut, aes(cutoff, pval))
		d <- d + geom_col() + labs(title = paste(gp_val, db_name, pamst), x ="cutoff", y = "adjust P value") + coord_cartesian(xlim=c(0, 1)) + theme_classic()
		ggsave(d, file = paste(res_dir, gp_val, "_pval_barplot.tiff"),  dpi = 300, width = 9, height = 6, units = "in", device = "tiff")
		
	} else {
		c <- ggplot(pval_cut, aes(cutoff, pval)) 
		c <- c + geom_col() + labs(title = paste(gp_val, db_name, pamst), x ="cutoff", y = "adjust P value") + coord_cartesian(xlim=c(0, 1)) + theme_classic()
		c <- c + geom_vline(xintercept = pval_cut$cutoff[pval_cut$pval == max(pval_cut$pval)], size = 0.5, colour = "red",linetype = "dotdash")
		ggsave(c, file = paste(res_dir, gp_val, "_pval_barplot.tiff"),  dpi = 300, width = 9, height = 6, units = "in", device = "tiff")
	}
	return(c) 
}

suppressMessages(library(genefu))
suppressMessages(library(readxl))
suppressMessages(library(ggplot2))
suppressMessages(library(survival))
suppressMessages(library(survminer))


work_dir = "//Bri-net/citi/Peter Lee Group/Weihua/surv_validation/"
#db_name = "metabric"
db_name = "tcga_brca"
sg_name = "tex_brtissue" # Colt's Tex sig from breast tissue c2
pval_name = "cut_pval"
expr_type = "" # ilid: raw data from EGA, median: raw median data from cbioportal, medianz: zscore from cbioportal
selfmap = TRUE # NOTE: ilid/tcga requires this as TRUE; median as FALSE

if (db_name == "metabric") {
	if (expr_type == "ilid") {expr_file = "metabric_expr_ilid.RDS"}
	if (expr_type == "median") {expr_file = "data_expression_median.RDS"}
	if (expr_type == "medianz") {expr_file = "data_mRNA_median_Zscores.RDS"}
	clin_rds = "merge_clin_info_v3.RDS" # clinical information with merged disease-free survival
	annot_file = "HumanHT_12_v30_R3_cleaned_v2.xlsx" # Microarray/Genome annotation
}
if (db_name == "tcga_brca") {
expr_file = "tcga_brca_log2trans_fpkm_uq_v2.RDS" # Expression file
clin_rds = "08272019_tcga_pam50_clin.xlsx" # clinical information with merged disease-free survival 
annot_file = "gencode.gene.info.v22.xlsx" # Microarray/Genome annotation
}

data_dir = paste(work_dir, db_name, "/", sep = "") # generate the directory with all the public data
sign_dir = paste(work_dir, sg_name, "/", sep = "") # generate the directory with signatures and corresponding results
pval_dir = paste(work_dir, pval_name, "/", sep = "")

sign_file = "tex_signature_colt_c2.txt" # Signature file Colt's Tex

histype = "IDC" # histology type: IDC/DCIS
pamst = c("LumA","LumB") # PAM50 status: LumA/LumB/Basal/Normal/Her2
hrtype = "" #c("N", "N", "N") # N: Negative, P: Positive, "-": DON'T CARE
gp_gene = "CD8A"
gptype = "CD8A"
qcut <- seq(from = 0.05, to = 0.95, by = 0.001) # screen cutoff range
gp_app = "oneqcut"#"symqcut" # oneqcut: one quantile cutoff (upper percential), symqcut: symmetric quantile cutoff
pval = "cox_reg" #  p value get from "surv" or "cox_reg"
gp_val = c("sig_score", "CD8A")

res_folder = "CD8A_cox_LumA_B_IDC_tcga_corr" # NOTE: Please change this folder name to identify your experiments
res_dir = paste(pval_dir, res_folder, "/", sep ="")
dir.create(file.path(pval_dir, res_folder), showWarnings = FALSE)

cat("Loading expression data...\n")
st = Sys.time()
#expr = readRDS(paste(data_dir, expr_file, sep = ""))
#expr = readRDS("metabric_expr_ilid.RDS")
#expr = readRDS("data_expression_median.RDS") # When test the script using cBioportal
expr = readRDS("tcga_brca_log2trans_fpkm_uq_v2.RDS") # When test the script using tcga
print(Sys.time()-st)
expr<-as.data.frame(expr)

cat("Loading clinical data...\n")
#clin_info = readRDS(paste(data_dir, clin_rds, sep = ""))
#clin_info = readRDS("merge_clin_info_v3.RDS")
clin_info = read_excel("08272019_tcga_pam50_clin.xlsx", sheet = 1)



cat("Start to filter by clinical info...\n")
sub_clin = clin_info

cat("\tOriginal patient number: ", dim(sub_clin)[1], "\n")
if (histype != "") {
	sub_clin = sub_clin[sub_clin[,"oncotree_code"] == histype,]
	sub_clin = sub_clin[complete.cases(sub_clin$pid),]
	cat("\tFiltered patient number: ", dim(sub_clin)[1], "\n")
}

if (pamst != "") {
	cat("Using PAM50 as molecular subtype classifier: ", pamst, "\n")
	#sub_clin = sub_clin[sub_clin[,"Pam50Subtype"] %in% pamst,]
	sub_clin = subset(sub_clin, Pam50Subtype %in% pamst)
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
annot = as.data.frame(read_excel(paste(data_dir, annot_file, sep = ""))) 

cat("Load gene signature...\n")
sign = read.table(paste(sign_dir, sign_file, sep = ""), header = TRUE, row.names = 1)


#################################################################################
if (selfmap) {
	ssdata = t(expr) # sig.score data
	ssdata = as.data.frame(ssdata)
	cat("Expression data dimension, sample: probe ", dim(ssdata), "\n")

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


#################################################################################
cat("Generate subtype signature score and survival data\n")
sc_res = as.data.frame(sig_score$score)
colnames(sc_res) = "sig_score"
sc_res$pid = rownames(sc_res)
sub_scres = sc_res[sc_res$pid %in% sub_clin$pid,]
sub_sigscore = sub_scres

if (gp_gene != "") {
	cat("Using ", gp_gene, " as the group criteria...\n")
	sub_expr = singene_expr(gene = gp_gene, expr = expr, annot = annot, subdf = sub_scres, caltype = "mean", map = selfmap)
	sub_scres = sub_expr
	sub_scres[,"gpvalue"] = sub_scres[,gp_gene]
} else {
	sub_scres[,"gpvalue"] = sub_scres[,"sig_score"]
}

#################################################################################
## Assign survival data
cat("Extract survival and treatment information from clinical data to subtype sig.score data frame\n")
sub_scres$ost = sub_clin$T[match(sub_scres$pid, sub_clin$pid)]
sub_scres$ose = as.numeric(sub_clin$DeathBreast[match(sub_scres$pid, sub_clin$pid)])
## Add all the other factors
#sub_scres$age = as.numeric(sub_clin$age_at_diagnosis[match(sub_scres$pid, sub_clin$pid)])
#sub_scres$grade = as.numeric(sub_clin$grade[match(sub_scres$pid, sub_clin$pid)])
#sub_scres$tsize = as.numeric(sub_clin$tumor_size[match(sub_scres$pid, sub_clin$pid)])
#sub_scres$node_stat = as.numeric(sub_clin$Lymph.Nodes.Positive[match(sub_scres$pid, sub_clin$pid)])


#################################################################################

tex_plot <- cutoff_screen(data = sub_scres, type = pval, qcut = qcut, gp_app = gp_app, gp_val = gp_val[1])
gene_plot <- cutoff_screen(data = sub_scres, type = pval, qcut = qcut, gp_app = gp_app, gp_val = gp_val[2])

corr_plot <- ggscatter(sub_scres, x = "sig_score", y = "CD8A", # genes correlate with sig.score
           color = "gray40", shape = 19, size = 1.5,
           add = "reg.line",  # Add regressin line
           add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
           conf.int = TRUE, # Add confidence interval
           cor.coef = TRUE, # Add correlation coefficient.
           cor.coeff.args = list(method = "spearman", label.x = 3, label.sep = "\n")) + border() + theme(plot.margin = margin(0, 0, 0, 0, "pt"))

tex_plot<- tex_plot + clean_theme() +labs(title = NULL)
gene_plot <- gene_plot + clean_theme () + labs(title = NULL) + rotate()
combine_plot<-ggarrange(tex_plot, NULL, corr_plot, gene_plot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)
combine_tiff <- paste(res_dir, gp_val[1], gp_val[2], "combine_plot.tiff", sep = "_")
ggsave(combine_plot, file = combine_tiff, dpi = 300, width = 11, height = 10, units = "in", device = "tiff" )