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

survana <- function(data, type, gptype = "Sig.score", plot = "", csv = "", 
		    cox = "", coxfac = c("age","grade","tsize","node_stat"), multicox = TRUE) {
	# NOTE: plot/csv/cox should be a folder directory where to store the results
	# NOTE: gptype should be the definition or basis of the grouping method, default is signature score.

	if(type == "os") {
		cat("Analyzing overall survival...\n")
		surv_data = data[,c("group","ost", "ose", coxfac)]
		lab = "Overall survival"
	}
	if(type == "rfs") {
		cat("Analyzing relapse-free survival...\n")
		surv_data = data[,c("group","rfst","rfse", coxfac)]
		lab = "Relapse-free survival"
	}
	colnames(surv_data)[2:3] <- c("time", "status")
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
		if (multicox) {
			cat("\tMulti-variant Cox analysis with ", coxfac, "\n")
			res.cox <- coxph(Surv(time, status) ~ group+age+grade+tsize+node_stat, data=surv_data)
		} else {
			cat("\tSingle-variant Cox analysis with ", gptype, "\n")
			res.cox <- coxph(Surv(time, status) ~ group, data=surv_data)
		}
		prescox = summary(res.cox) # Print results of COX
		sum_txt = paste(cox, lab, gptype, "_summary_cox_results.txt", sep = "")
		sink(sum_txt)
		print(summary(res.cox))
		sink()
		cox_rds <- paste(cox, lab, gptype, "_survana_cox_res.rds", sep="")
		cat("\t\tP-value from logRank test: ", (prescox$sctest["pvalue"]), "\n")
		saveRDS(res.cox, file=cox_rds)}
	if (csv != "") {
		surv_rescsv <- paste(csv, lab, gptype, "_survana_res.csv", sep="")
		write.csv(fit_df, file=surv_rescsv)}
	if (plot != "") {
		surv_plot <- paste(plot, lab, gptype, "_survana_curves.tiff", sep="")
		survcurv <- ggsurvplot(fit, data = surv_data,
				       xlim = c(0,10000), ylim = c(0.00, 1.00),
				       pval = TRUE, pval.size = 6, pval.coord = c(0, 0.2),
				       conf.int = TRUE, conf.int.alpha = 0.2,
				       xlab = "Time (days)", ylab = lab, legend.title = gptype,
				       legend.labs = c("Low", "High"),
#                                      surv.median.line = "hv",
				       ggtheme = theme_classic(),
				       palette = c("#E7B800", "#2E9FDF"), 
				       font.x = c(14, "bold"), font.y = c(18, "bold"), 
				       font.tickslab = c(16, "plain"), font.legend = c(18, "bold"),
				       risk.table = FALSE, ncensor.plot = FALSE)
		note_on_plot <- paste("nHigh = ", numh, "\t", "nLow = ", numl, "\t")
		survcurv$plot <- survcurv$plot + annotate("text", x=2000, y=0.1, label=note_on_plot, size = 6)
		ggsave(surv_plot, plot = survcurv$plot, dpi = 120, width = 9, height = 6, units = 'in')
	}
}

# Please load all the packages at the very begining of each script
cat("Loading genefu library...\n")
suppressMessages(library(genefu))
suppressMessages(library(readxl))
suppressMessages(library(ggplot2))
suppressMessages(library(survival))
suppressMessages(library(survminer))

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
clin_rds = "merge_clin_info_v3.RDS" # clinical information with merged disease-free survival
annot_file = "HumanHT_12_v30_R3_cleaned_v2.xlsx" # Microarray/Genome annotation
# sign_file = "trm_tex_brtissue_only_all_markers.xlsx" # Signature file
sign_file = "loi_trm_signature.txt" # Signature file
# 

pamst = "LumB"
gp_app = "oneqcut"
# gp_app = "symqcut"

qcut = 0.25 # This is TOP quantile for oneqcut approach

# Work for experiment records
res_folder = "top25bot75_loi_trm_LumB_all_metabric" # NOTE: Please change this folder name to identify your experiments
res_dir = paste(sign_dir, res_folder, "/", sep ="")
dir.create(file.path(sign_dir, res_folder), showWarnings = FALSE)
# COPY the used script to the result folder for recording what experiment was run
script_dir = "/home/weihua/git_repo/surv_analysis/"
script_name = "surv_public_test.R"
file.copy(paste(script_dir, script_name, sep = ""), res_dir)


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
print(Sys.time()-st)

cat("Start to filter by clinical info...\n")
sub_clin = clin_info
cat("\tOriginal patient number: ", dim(sub_clin)[1], "\n")
# For IDC
# sub_clin = sub_clin[sub_clin[,"oncotree_code"] %in% c("IDC"),]
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
## RUN sig.score
cat("Run sig.score in all patient samples\n")
sig_score <- sig.score(x=sssign, data=ssdata, annot=ssannot, do.mapping=TRUE, signed=TRUE, verbose=TRUE)
# print(head(sig_score$score))
if (FALSE) {
	st = Sys.time()
	cat("Saving sig.score input and output to ", ssin, "\n")
	saveRDS(ssin, file = paste(res_dir, "sig_score_inputs.RDS", sep = ""))
	saveRDS(sig_score, file = paste(res_dir, "sig_score_outputs.RDS", sep = ""))
	print(Sys.time()-st)
}


cat("Generate subtype signature score and survival data\n")
sc_res = as.data.frame(sig_score$score)
colnames(sc_res) = "sig_score"
sc_res$pid = rownames(sc_res)
sub_scres = sc_res[sc_res$pid %in% sub_clin$pid,]


## Histogram for sig.score
cat("Generate histogram plot of signature score\n")
title = paste(db_name, sg_name, pamst, sep = "  ")
sc_hist = ggplot(sub_scres, aes(x=sig_score)) + 
	geom_histogram(color="darkblue", fill="lightblue", binwidth = 0.036) +
	labs(title=title, x="sig.score", y = "Count") + 
	theme_classic()
sc_hist_bld = ggplot_build(sc_hist)
sc_hist_data = sc_hist_bld$data[[1]]
zero_x = sc_hist_data[sc_hist_data$count == 0,"x"]
right_zero_count = sum(sc_hist_data[sc_hist_data$x>=zero_x[1],"count"])
print(right_zero_count)
sc_hist = sc_hist + geom_vline(xintercept = zero_x[1], size = 1, colour = "grey",linetype = "dashed")
if (gp_app == "oneqcut") {
	cat("Group the patient by one quantile cutoff: ", qcut,"\n")
	qcov = quantile(sub_scres$sig_score, c(1-qcut)) # quantile cutoff value
	# Add line for cutoff value
	sc_hist = sc_hist + 
		geom_vline(xintercept = qcov[[1]], size = 1, colour = "purple",linetype = "dotdash")
	sc_hist = sc_hist + annotate("text", label = paste("Single cutoff value:\n", format(qcov[[1]],digit = 3), "(", 1-qcut, ")"),
				     hjust = 0, x = qcov[[1]], y = max(sc_hist_data$count)*0.96, size = 4.5, colour = "black")
	hist_tif = paste(res_dir, db_name, sg_name, pamst, "single_quantile_cutoff.tiff", sep = "_")

}
if (gp_app == "symqcut") {
	cat("Group the patient by one quantile cutoff with symmetric manner: ", qcut,"\n")

	qcov = quantile(sub_scres$sig_score, c(qcut, 1-qcut)) # quantile cutoff value

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
ggsave(sc_hist, file = hist_tif, width = 9, height = 6, units = "in")

## Assign survival data
cat("Extract survival information from clinical data to subtype sig.score data frame\n")
sub_scres$ost = 0
sub_scres$ose = 0
sub_scres$rfst = 0
sub_scres$rfse = 0
# For overall survival
sub_scres[sub_clin$pid,"ost"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"T"]
sub_scres[sub_clin$pid,"ose"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"DeathBreast"]
# For disease-free survival
sub_scres[sub_clin$pid,"rfst"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"TOR"]
sub_scres[sub_clin$pid,"rfse"] = sub_clin[sub_clin$pid %in% sub_scres$pid,"OR"]

## Add all the other factors
sub_scres[sub_clin$pid,"age"] = as.numeric(sub_clin[sub_clin$pid %in% sub_scres$pid,"age_at_diagnosis"])
sub_scres[sub_clin$pid,"grade"] = as.numeric(sub_clin[sub_clin$pid %in% sub_scres$pid,"grade"])
sub_scres[sub_clin$pid,"tsize"] = as.numeric(sub_clin[sub_clin$pid %in% sub_scres$pid,"tumor_size"])
sub_scres[sub_clin$pid,"node_stat"] = as.numeric(sub_clin[sub_clin$pid %in% sub_scres$pid,"Lymph.Nodes.Positive"])


## Assign groups
if (length(qcov) == 1) {
	sub_scres$group = "Medium"
	sub_scres[sub_scres$sig_score <= qcov[[1]],"group"] = "Low"
	sub_scres[sub_scres$sig_score > qcov[[1]],"group"] = "High"
	surv_csv = paste(res_dir, "single_qcut_survival_analysis_inputs.csv")
	write.csv(sub_scres, file = surv_csv)
}
if (length(qcov) == 2) {
	sub_scres$group = "Medium"
	sub_scres[sub_scres$sig_score <= qcov[[1]],"group"] = "Low"
	sub_scres[sub_scres$sig_score >= qcov[[2]],"group"] = "High"
	surv_csv = paste(res_dir, "symmetric_qcut_survival_analysis_inputs.csv")
	write.csv(sub_scres, file = surv_csv)
	sub_scres = sub_scres[sub_scres$group != "Medium",]
}
if (length(qcov) > 2) {stop("Mulitple cutoffs!!!")}
print(head(sub_scres))


survana(data = sub_scres, type = "os", plot = res_dir, gptype = "TRM sig.score", cox = res_dir)
survana(data = sub_scres, type = "rfs", plot = res_dir, gptype = "TRM sig.score", cox = res_dir)

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

# cat("Generate exhuasted T cell sig. score in subtype\n")
# sub_score<-subset(ET_score, ET_score$patient_id %in% subtype_clin$pid)

# cat("Extract survival information from clinical data to subtype sig.score data frame")
# sub_score<-sub_score[order(-sub_score$sig.score),]
# sub_score$OS<-subtype_clin$T[match(subtype_clin$pid, sub_score$patient_id)]
# sub_score$os_status<-subtype_clin$DeathBreast[match(subtype_clin$pid, sub_score$patient_id)]
# subtype_clin$DFS<-ifelse(subtype_clin$TLR>subtype_clin$TDR,subtype_clin$TDR,subtype_clin$TLR)
# subtype_clin$dfs_status<-ifelse(c(subtype_clin$LR==0 & subtype_clin$DR==0),0,1) 
# sub_score$DFS<-subtype_clin$DFS[match(sub_score$patient_id,subtype_clin$pid)]
# sub_score$dfs_status<-subtype_clin$dfs_status[match(sub_score$patient_id,subtype_clin$pid)]

# stop("JUST STOP!!!!!!!")

# cat("Select patient samples by 25/25 cutoff of exhuasted T sig.score")

# library(dplyr)
# top_25<-subset(sub_score, sig.score > quantile(sig.score, prob = 1 - 25/100))
# bot_25<-subset(sub_score, sig.score < quantile(sig.score, prob = 25/100))
# top_25$score="t25%"
# bot_25$score="b25%"
# sub_score<-rbind(top_25,bot_25)



# filename<- paste("c:/Users/jitan/Documents/ET_", subtype, ".xlsx", sep = "")
# library(xlsx)
# write.xlsx(sub_score, filename)




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
