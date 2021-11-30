### Screen for PD-L1- TAM signature in METABRIC
### Weihua Guo, Ph.D.
### 11/02/21

rm(list = ls())
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(survMisc))


data_dir <- "C:/Users/wguo/OneDrive - City of Hope National Medical Center/tmp_works/tam_surv/METABRIC/brca_metabric"
sig_dir <- "C:/Users/wguo/OneDrive - City of Hope National Medical Center/tmp_works/tam_surv"
hr_subtype <- "ER" # ER, regardless of HER2; HER2+, ER-/PR-/HER2+; TNBC, ER-/PR-/HER2-
oncotree <- "IDC" # ALL: use all the samples
strict <- TRUE

res_folder <- paste(hr_subtype, oncotree, "pdl1_neg_tam_rf_metabric", sep = "_") # Change here to change the output folder names
dir.create(file.path(sig_dir, res_folder), showWarnings = FALSE)
ppf <- paste(sig_dir, "/", res_folder, "/", hr_subtype, "_", oncotree, "_", sep = "")

cat("Reading METABRIC data...\n")
expr <- read.table(paste(data_dir, "/data_expression_median.txt", sep = ""), header = T, sep = "\t",
                   check.names = F, stringsAsFactors = F)
dup_genes <- as.data.frame(table(expr[,1]))
dup_genes <- dup_genes[dup_genes$Freq > 1,]
colnames(expr)[1] <- "gene"

pcln <- read.table(paste(data_dir, "/data_clinical_patient.txt", sep = ""), header = T, sep = "\t",
                   check.names = F, stringsAsFactors = F)
scln <- read.table(paste(data_dir, "/data_clinical_sample.txt", sep = ""), header = T, sep = "\t",
                   check.names = F, stringsAsFactors = F)

mcln <- merge(pcln, scln, by = "PATIENT_ID", all.y = T)
mcln$flag <- TRUE

cat("Filtering the clinical data...\n")
if (hr_subtype == "ER") {
  hr_mask <- mcln$ER_STATUS == "Positive"
} else if (hr_subtype == "HER2") {
  hr_mask <- (mcln$ER_STATUS == "Negative") & (mcln$HER2_STATUS == "Positive") & (mcln$PR_STATUS == "Negative")
} else if (hr_subtype == "TNBC") {
  hr_mask <- (mcln$ER_STATUS == "Negative") & (mcln$HER2_STATUS == "Negative") & (mcln$PR_STATUS == "Negative")
} else {
  hr_mask <- mcln$flag
}

if (oncotree == "IDC") {
  octr_mask <- mcln$ONCOTREE_CODE == "IDC"
} else {
  octr_mask <- mcln$flag
}

use_cln <- mcln[hr_mask & octr_mask,]
use_cln <- use_cln[!is.na(use_cln$RFS_STATUS),]
use_cln <- use_cln[!is.na(use_cln$RFS_MONTHS),]

cat("Reading signature...\n")
sign <- read.table(paste(sig_dir, "/tam_pdl1_neg_mast_loose_marker_v1.txt", sep = ""), header = T, sep = "\t",
                   check.names = F, stringsAsFactors = F)
if (strict) {
  sign <- sign[abs(sign$logfc) >= 1,]
  ppf <- paste(ppf, "strict_", sep = "")
}
cat("Merging all the data...\n")
# sig_genes <- sign_expr$gene
sig_genes <- c("SPP1", "FABP5", "COL1A2") # CHANGE HERE
use_sign <- sign[sign$gene %in% intersect(sig_genes, expr$gene),]
sign_expr <- expr[expr$gene %in% intersect(sig_genes, expr$gene),]

##### MCPcounter ####
use_expr <- sign_expr[,intersect(use_cln$PATIENT_ID, colnames(sign_expr))]
use_cln <- use_cln[use_cln$PATIENT_ID %in% colnames(use_expr),]
rownames(use_cln) <- use_cln$PATIENT_ID
use_df <- use_cln[colnames(use_expr),]
use_df$MCP_sign <- colMeans(use_expr)
use_df$RFS_STATUS <- ifelse(str_detect(use_df$RFS_STATUS, "0"), 0, 1)

fit <- coxph(Surv(RFS_MONTHS, RFS_STATUS) ~ AGE_AT_DIAGNOSIS + GRADE + TUMOR_SIZE + TUMOR_STAGE + MCP_sign, data = use_df)
forest_gg <- ggforest(fit)
ggsave(paste(ppf, 'rfs_MCP_multivariant.png', sep = ''), forest_gg, dpi = 300, width = 6, height = 4)

fit <- coxph(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_sign, data = use_df)
cutp_res <- cutp(fit)

out_df <- use_df
out_df <- merge(use_df, as.data.frame(cutp_res[[1]]), by = "MCP_sign")
write.csv(out_df, paste(ppf, "rfs_MCP_use_df_with_cutp.csv", sep = ""))

use_df$MCP_group50 <- ifelse(use_df$MCP_sign >= quantile(use_df$MCP_sign, 0.5), "High", "Low")
fit <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_group50, data = use_df)
surv_gg <- ggsurvplot(fit, data = use_df, pval = TRUE,
                      title = "50/50",
                      legend.title = 'PD-L1- TAM',
                      legend.labs = c('High', 'Low'))
png(paste(ppf, 'rfs_MCP_50.png', sep = ''), res = 300, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

use_df$MCP_group2575 <- ifelse(use_df$MCP_sign >= quantile(use_df$MCP_sign, 0.25), "High", "Low")
fit <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_group2575, data = use_df)
surv_gg <- ggsurvplot(fit, data = use_df, pval = TRUE,
                      title = "25/75",
                      legend.title = 'PD-L1- TAM',
                      legend.labs = c('High', 'Low'))
png(paste(ppf, 'rfs_MCP_2575.png', sep = ''), res = 300, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

use_df$MCP_group7525 <- ifelse(use_df$MCP_sign >= quantile(use_df$MCP_sign, 0.75), "High", "Low")
fit <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_group7525, data = use_df)
surv_gg <- ggsurvplot(fit, data = use_df, pval = TRUE,
                      title = "75/25",
                      legend.title = 'PD-L1- TAM',
                      legend.labs = c('High', 'Low'))
png(paste(ppf, 'rfs_MCP_7525.png', sep = ''), res = 300, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

use_df$MCP_group2525 <- ifelse(use_df$MCP_sign >= quantile(use_df$MCP_sign, 0.75), "High", 
                               ifelse(use_df$MCP_sign <= quantile(use_df$MCP_sign, 0.25), "Low", "Int"))
tmp_df <- use_df[use_df$MCP_group2525 != "Int",]
fit <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_group2525, data = tmp_df)
surv_gg <- ggsurvplot(fit, data = tmp_df, pval = TRUE,
                      title = "25/25",
                      legend.title = 'PD-L1- TAM',
                      legend.labs = c('High', 'Low'))
png(paste(ppf, 'rfs_MCP_2525.png', sep = ''), res = 300, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

use_df$MCP_group3030 <- ifelse(use_df$MCP_sign >= quantile(use_df$MCP_sign, 0.70), "High", 
                               ifelse(use_df$MCP_sign <= quantile(use_df$MCP_sign, 0.30), "Low", "Int"))
tmp_df <- use_df[use_df$MCP_group3030 != "Int",]
fit <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_group3030, data = tmp_df)
surv_gg <- ggsurvplot(fit, data = tmp_df, pval = TRUE,
                      title = "30/30",
                      legend.title = 'PD-L1- TAM',
                      legend.labs = c('High', 'Low'))
png(paste(ppf, 'rfs_MCP_3030.png', sep = ''), res = 300, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

## TODO: average/CD68


##### MCPcounter/CD68 ####
use_expr <- sign_expr[,intersect(use_cln$PATIENT_ID, colnames(sign_expr))]
use_cln <- use_cln[use_cln$PATIENT_ID %in% colnames(use_expr),]
rownames(use_cln) <- use_cln$PATIENT_ID
use_df <- use_cln[colnames(use_expr),]
use_df$MCP_sign <- colMeans(use_expr)
use_df$MCP_sign <- use_df$MCP_sign/t(expr[expr$gene == "CD68", colnames(use_expr)])
use_df$RFS_STATUS <- ifelse(str_detect(use_df$RFS_STATUS, "0"), 0, 1)

fit <- coxph(Surv(RFS_MONTHS, RFS_STATUS) ~ AGE_AT_DIAGNOSIS + GRADE + TUMOR_SIZE + TUMOR_STAGE + MCP_sign, data = use_df)
forest_gg <- ggforest(fit)
ggsave(paste(ppf, 'rfs_MCP_CD68_multivariant.png', sep = ''), forest_gg, dpi = 300, width = 6, height = 4)

use_df$MCP_group50 <- ifelse(use_df$MCP_sign >= quantile(use_df$MCP_sign, 0.5), "High", "Low")
fit <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_group50, data = use_df)
surv_gg <- ggsurvplot(fit, data = use_df, pval = TRUE,
                      title = "50/50",
                      legend.title = 'PD-L1- TAM',
                      legend.labs = c('High', 'Low'))
png(paste(ppf, 'rfs_MCP_CD68_50.png', sep = ''), res = 300, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

use_df$MCP_group2575 <- ifelse(use_df$MCP_sign >= quantile(use_df$MCP_sign, 0.25), "High", "Low")
fit <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_group2575, data = use_df)
surv_gg <- ggsurvplot(fit, data = use_df, pval = TRUE,
                      title = "25/75",
                      legend.title = 'PD-L1- TAM',
                      legend.labs = c('High', 'Low'))
png(paste(ppf, 'rfs_MCP_CD68_2575.png', sep = ''), res = 300, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

use_df$MCP_group7525 <- ifelse(use_df$MCP_sign >= quantile(use_df$MCP_sign, 0.75), "High", "Low")
fit <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_group7525, data = use_df)
surv_gg <- ggsurvplot(fit, data = use_df, pval = TRUE,
                      title = "75/25",
                      legend.title = 'PD-L1- TAM',
                      legend.labs = c('High', 'Low'))
png(paste(ppf, 'rfs_MCP_CD68_7525.png', sep = ''), res = 300, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

use_df$MCP_group2525 <- ifelse(use_df$MCP_sign >= quantile(use_df$MCP_sign, 0.75), "High", 
                               ifelse(use_df$MCP_sign <= quantile(use_df$MCP_sign, 0.25), "Low", "Int"))
tmp_df <- use_df[use_df$MCP_group2525 != "Int",]
fit <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_group2525, data = tmp_df)
surv_gg <- ggsurvplot(fit, data = tmp_df, pval = TRUE,
                      title = "25/25",
                      legend.title = 'PD-L1- TAM',
                      legend.labs = c('High', 'Low'))
png(paste(ppf, 'rfs_MCP_CD68_2525.png', sep = ''), res = 300, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

use_df$MCP_group3030 <- ifelse(use_df$MCP_sign >= quantile(use_df$MCP_sign, 0.70), "High", 
                               ifelse(use_df$MCP_sign <= quantile(use_df$MCP_sign, 0.30), "Low", "Int"))
tmp_df <- use_df[use_df$MCP_group3030 != "Int",]
fit <- survfit(Surv(RFS_MONTHS, RFS_STATUS) ~ MCP_group3030, data = tmp_df)
surv_gg <- ggsurvplot(fit, data = tmp_df, pval = TRUE,
                      title = "30/30",
                      legend.title = 'PD-L1- TAM',
                      legend.labs = c('High', 'Low'))
png(paste(ppf, 'rfs_MCP_CD68_3030.png', sep = ''), res = 300, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

##### Univariate screening ####
use_expr <- rbind(sign_expr, expr[expr$gene %in% c("CD68", "CD274"),])
rownames(use_expr) <- use_expr$gene
use_expr <- use_expr[,intersect(use_cln$PATIENT_ID, colnames(use_expr))]
use_cln <- use_cln[use_cln$PATIENT_ID %in% colnames(use_expr),]

sig_genes <- rownames(use_expr)
use_df <- as.data.frame(t(use_expr))
use_df$PATIENT_ID <- rownames(use_df)
use_df <- merge(use_df, use_cln, by = "PATIENT_ID")
use_df$RFS_STATUS <- ifelse(str_detect(use_df$RFS_STATUS, "0"), 0, 1)

cat('\t\tForest plot -- Univariant...MANUALL\n')
#https://www.r-bloggers.com/2016/12/cox-proportional-hazards-model/
covariates <- sig_genes
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(RFS_MONTHS, RFS_STATUS)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = use_df)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$logtest["pvalue"], digits=2)
                         logrank.test<-signif(x$logtest["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         res<-c(beta, HR, HR.confint.lower,HR.confint.upper, logrank.test, p.value)
                         names(res)<-c("beta", "HR", "HRL", "HRU", "logrank", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res$sig <- sig_genes
res$p_format <- str_c('p-val = ', res$p.value)
write.csv(res, paste(ppf, 'rfs_univariant.csv', sep = ''))

form_gg <- ggplot(res, aes(x = HR, y = sig)) +
  geom_point() +
  geom_errorbar(aes(xmin = HRL, xmax = HRU), width = 0.2) +
  geom_vline(xintercept = 1.0, linetype = 'dashed', color = 'red') +
  geom_text(aes(label = p_format), hjust = -1, vjust = -1.5) +
  labs(x = 'Hazard ratio (univariant tests, overall survival)', y = 'Signatures') +
  theme_classic()
ggsave(paste(ppf, 'rfs_univariant.png', sep = ''), form_gg, dpi = 300, width = 6, height = 9)

##### Per CD68 ####
ratio_df <- use_df
ratio_df[,sig_genes] <- ratio_df[,sig_genes]/ratio_df[,"CD68"]
covariates <- sig_genes[sig_genes != "CD68"]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(RFS_MONTHS, RFS_STATUS)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = ratio_df)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$logtest["pvalue"], digits=2)
                         logrank.test<-signif(x$logtest["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         res<-c(beta, HR, HR.confint.lower,HR.confint.upper, logrank.test, p.value)
                         names(res)<-c("beta", "HR", "HRL", "HRU", "logrank", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res$sig <- covariates
res$p_format <- str_c('p-val = ', res$p.value)
write.csv(res, paste(ppf, 'rfs_univariant_per_cd68.csv', sep = ''))

form_gg <- ggplot(res, aes(x = HR, y = sig)) +
  geom_point() +
  geom_errorbar(aes(xmin = HRL, xmax = HRU), width = 0.2) +
  geom_vline(xintercept = 1.0, linetype = 'dashed', color = 'red') +
  geom_text(aes(label = p_format), hjust = -1, vjust = -1.5) +
  labs(x = 'Hazard ratio (univariant tests, overall survival)', y = 'Signatures') +
  theme_classic()
ggsave(paste(ppf, 'rfs_univariant_per_cd68.png', sep = ''), form_gg, dpi = 300, width = 6, height = 9)

