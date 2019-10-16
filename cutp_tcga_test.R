



suppressMessages(library(readxl))
suppressMessages(library(ggplot2))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(survMisc))


work_dir = "//Bri-net/citi/Peter Lee Group/Weihua/surv_validation/"
db_name = "tcga_brca"
sg_name = "ching_validation_res" # Colt's Tex sig from breast tissue c2
#pval_name = "cut_pval"
expr_type = "" # ilid: raw data from EGA, median: raw median data from cbioportal, medianz: zscore from cbioportal
cancer_type = "crc" # "crc"
selfmap = TRUE # NOTE: ilid/tcga requires this as TRUE; median as FALSE



gp_gene = "CD274"
gptype = "Tex sig.score"#"CD8A"
#pval = "surv" #  p value get from "surv" or "cox_reg"
#gp_val = "sig_score"#c("sig_score", "CD8A")
corr_plot = TRUE

if (cancer_type == "skcm") {
	expr_file = "tcga_skcm_expr_fpkm.rds"
	clin_file = "tcga_clin_skcm_clean_v4.xlsx"
}
if (cancer_type == "crc") {
expr_file = "tcga_crc_expr_fpkm.rds" # Expression file
clin_file = "tcga_clin_crc_clean_v4.xlsx" # clinical information with merged disease-free survival
} 
annot_file = "gencode.gene.info.v22.xlsx" # Microarray/Genome annotation


data_dir = paste(work_dir, db_name, "/", sep = "") # generate the directory with all the public data
sign_dir = paste(work_dir, sg_name, "/", sep = "") # generate the directory with signatures and corresponding results


cat("Loading expression data...\n")
st = Sys.time()
expr = readRDS(paste(data_dir, expr_file, sep = ""))
print(Sys.time()-st)

cat("Loading clinical data...\n")
clin_info = read_excel(paste(data_dir, clin_file, sep = ""))
#clin_info = readRDS("merge_clin_info_v3.RDS")
#clin_info = read_excel("08272019_tcga_pam50_clin.xlsx", sheet = 1)

cat("Loading genome annotation data...\n")
annot = as.data.frame(read_excel(paste(data_dir, annot_file, sep = ""))) 

gene_info = subset(annot, ILMN_Gene == gp_gene)
if (dim(gene_info)[1] == 0) {stop("No MATCHED gene found!!!")}
if (sum(gene_info$Probe_Id %in% rownames(expr)) == 0) {stop("NO MATCHED probe found!!!")}
gene_probe = expr[gene_info$Probe_Id,]
gene_probe["pid",] = colnames(gene_probe)
gene_expr <- as.data.frame(t(gene_probe))
colnames(gene_expr)[1] = gp_gene
gene_expr[,"gpvalue"] = gene_expr[,gp_gene]
gene_expr$gpvalue <- as.character(gene_expr$gpvalue)
gene_expr$gpvalue <- as.numeric(gene_expr$gpvalue)

gene_expr$pid <- substr(gene_expr$pid, 1, 12)

cat("Extract survival and treatment information from clinical data to subtype sig.score data frame\n")
gene_expr$ost = as.numeric(clin_info$T[match(gene_expr$pid, clin_info$pid)])
gene_expr$ose = as.numeric(clin_info$ose[match(gene_expr$pid, clin_info$pid)])

res.cox <- coxph(Surv(ost, ose) ~ gpvalue, data=gene_expr)
res.cox <- cutp(res.cox)$gpvalue
data.table::setorder(res.cox, "gpvalue")
cutp_res <-(res.cox)[,"U"]
colnames(cutp_res)[1] = "test_score"
cutp_res$cutoff <- as.numeric(rownames(cutp_res))
g <- ggplot(cutp_res, aes(cutoff, test_score)) + geom_col() + labs(y = "log-rank test score") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
g <- g + geom_vline(xintercept = cutp_res$cutoff[cutp_res$test_score == max(cutp_res$test_score)], size = 0.5, colour = "red",linetype = "dotdash")
max_cutoff <- res.cox$gpvalue[res.cox$U == max(res.cox$U)]
print(max_cutoff)
stop()
ggsave(g, file = "skcm_CD8A_cutp.tiff",  dpi = 300, width = 15, height = 6, units = "in", device = "tiff")