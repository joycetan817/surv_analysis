# For TAM CD68
# Weihua Guo, Ph.D
# 02/04/2021

rm(list = ls())
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(cocor))
suppressMessages(library(scales))
suppressMessages(library(survminer))
suppressMessages(library(survival))
suppressMessages(library(survMisc))
suppressMessages(library(rstatix))

data_dir <- 'C:/Users/wguo/OneDrive - City of Hope National Medical Center/tmp_works/tam_surv/'
pdl1_status <- 'neg'
co <- "tp25"
cov <- 0.25
plot_pf <- paste(data_dir, pdl1_status, '/PDL1_', pdl1_status, '_TAM_', co, sep = '')

use_df <- read.csv(paste(data_dir, pdl1_status, '_tam.csv', sep = ''), row.names = 1)
clinical_df <- read.csv(paste(data_dir, 'brca_metabric_clinical_data.csv', sep = ''), 
                        header = T, row.names = 1, check.names = T)
rownames(clinical_df) <- str_replace(rownames(clinical_df), '-', '.')
clinical_df$Pam50...Claudin.low.subtype
use_df$pam50_claudin <- clinical_df[rownames(use_df), "Pam50...Claudin.low.subtype"]

nm_df <- use_df
nm_df$group_tam <- ifelse(nm_df$gpvalue > quantile(nm_df$gpvalue, 1-cov), 'High', 
                           ifelse(nm_df$gpvalue < quantile(nm_df$gpvalue, cov), 'Low', 'Medium'))
nm_df$group_cd68 <- ifelse(nm_df$CD68 > quantile(nm_df$CD68, 1-cov), 'High', 
                           ifelse(nm_df$CD68 < quantile(nm_df$CD68, cov), 'Low', 'Medium'))
nm_df$group_cd8A <- ifelse(nm_df$CD8A > quantile(nm_df$CD8A, 1-cov), 'High', 
                           ifelse(nm_df$CD8A < quantile(nm_df$CD8A, cov), 'Low', 'Medium'))
nm_df <- nm_df[nm_df$group_tam != "Medium",]



cat('\t\tTAM\n')
fit <- survfit(Surv(ost, ose) ~ group_tam, data = nm_df)
surv_gg <- ggsurvplot(fit, data = nm_df, pval = TRUE,
                      title = paste("PD-L1", pdl1_status, "-", co),
                      legend.title = 'TAM',
                      legend.labs = c('High', 'Low'),
                      risk.table = TRUE)
png(paste(plot_pf, '_nm_os_km.png', sep = '_'), res = 600, width = 6, height = 6, units = 'in')
print(surv_gg)
gar <- dev.off()

fit <- survfit(Surv(rfst, rfse) ~ group_tam, data = nm_df)
surv_gg <- ggsurvplot(fit, data = nm_df, pval = TRUE,
                      title = paste("PD-L1", pdl1_status, "-", co),
                      legend.title = 'TAM',
                      legend.labs = c('High', 'Low'),
                      risk.table = TRUE)
png(paste(plot_pf, '_nm_rfs_km.png', sep = '_'), res = 600, width = 6, height = 6, units = 'in')
print(surv_gg)
gar <- dev.off()

cat('Hazard ratio -- Multivariant\n')
fit <- coxph(Surv(ost, ose) ~ age+grade+tsize+pam50_claudin+gpvalue, data = nm_df)
sum_fit <- summary(fit)
sum_fit_df <- cbind(sum_fit$coefficients, sum_fit$conf.int)
write.csv(sum_fit_df, paste(plot_pf, '_multiv_os.csv', sep = '_'))
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, '_multiv_os.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)

cat('Hazard ratio -- Multivariant\n')
fit <- coxph(Surv(rfst, rfse) ~ age+grade+tsize+pam50_claudin+gpvalue, data = nm_df)
sum_fit <- summary(fit)
sum_fit_df <- cbind(sum_fit$coefficients, sum_fit$conf.int)
write.csv(sum_fit_df, paste(plot_pf, '_multiv_rfs.csv', sep = '_'))
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, '_multiv_rfs.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)

cat('\t\tTAM-CD68\n')
sub_df <- nm_df[nm_df$group_cd68 != 'Medium',]
sub_df$CD68_TAM <- str_c('CD68', sub_df$group_cd68, '-TAM', sub_df$group_tam)
fit <- survfit(Surv(ost, ose) ~ CD68_TAM, data = sub_df)
surv_gg <- ggsurvplot(fit, data = sub_df, pval = TRUE,
                      title = paste("PD-L1", pdl1_status, "-", co),
                      legend.title = 'CD68-TAM',
                      risk.table = TRUE)
png(paste(plot_pf, '_cd68_tam_os_km.png', sep = '_'), res = 600, width = 7.5, height = 6, units = 'in')
print(surv_gg)
gar <- dev.off()
pres <- pairwise_survdiff(Surv(ost, ose) ~ CD68_TAM, data=sub_df)
write.csv(pres$p.value, paste(plot_pf, '_pair_cd68_tam_os_bhadj.csv', sep = ''))

fit <- survfit(Surv(rfst, rfse) ~ CD68_TAM, data = sub_df)
surv_gg <- ggsurvplot(fit, data = sub_df, pval = TRUE,
                      title = paste("PD-L1", pdl1_status, "-", co),
                      legend.title = 'CD68-TAM',
                      risk.table = TRUE)
png(paste(plot_pf, '_cd68_tam_rfs_km.png', sep = '_'), res = 600, width = 7.5, height = 6, units = 'in')
print(surv_gg)
gar <- dev.off()
pres <- pairwise_survdiff(Surv(rfst, rfse) ~ CD68_TAM, data=sub_df)
write.csv(pres$p.value, paste(plot_pf, '_pair_cd68_tam_rfs_bhadj.csv', sep = ''))

cat('\t\tTAM-CD8A\n')
sub_df <- nm_df[nm_df$group_cd8A != 'Medium',]
sub_df$CD8A_TAM <- str_c('CD8A', sub_df$group_cd8A, '-TAM', sub_df$group_tam)
fit <- survfit(Surv(ost, ose) ~ CD8A_TAM, data = sub_df)
surv_gg <- ggsurvplot(fit, data = sub_df, pval = TRUE,
                      title = paste("PD-L1", pdl1_status, "-", co),
                      legend.title = 'CD8A-TAM',
                      risk.table = TRUE)
png(paste(plot_pf, '_cd8a_tam_os_km.png', sep = '_'), res = 600, width = 7.5, height = 6, units = 'in')
print(surv_gg)
gar <- dev.off()
pres <- pairwise_survdiff(Surv(ost, ose) ~ CD8A_TAM, data=sub_df)
write.csv(pres$p.value, paste(plot_pf, '_pair_cd8a_tam_os_bhadj.csv', sep = ''))

fit <- survfit(Surv(rfst, rfse) ~ CD8A_TAM, data = sub_df)
surv_gg <- ggsurvplot(fit, data = sub_df, pval = TRUE,
                      title = paste("PD-L1", pdl1_status, "-", co),
                      legend.title = 'CD8A-TAM',
                      risk.table = TRUE)
png(paste(plot_pf, '_cd8a_tam_rfs_km.png', sep = '_'), res = 600, width = 7.5, height = 6, units = 'in')
print(surv_gg)
gar <- dev.off()
pres <- pairwise_survdiff(Surv(rfst, rfse) ~ CD8A_TAM, data=sub_df)
write.csv(pres$p.value, paste(plot_pf, '_pair_cd8a_tam_rfs_bhadj.csv', sep = ''))

cat('\t\tTAM-PAM50\n')
pam_df <- nm_df[nm_df$pam50_claudin %in% c("LumA", "LumB"),]
pam_df$pam50_TAM <- str_c(pam_df$pam50, '-', pam_df$group_tam)
fit <- survfit(Surv(ost, ose) ~ pam50_TAM, data = pam_df)
surv_gg <- ggsurvplot(fit, data = pam_df, pval = TRUE,
                      title = paste("PD-L1", pdl1_status, "-", co),
                      legend.title = 'TAM',
                      risk.table = TRUE)
png(paste(plot_pf, '_pam50_os_km.png', sep = '_'), res = 600, width = 7.5, height = 6, units = 'in')
print(surv_gg)
gar <- dev.off()
pres <- pairwise_survdiff(Surv(ost, ose) ~ pam50_TAM, data=pam_df)
write.csv(pres$p.value, paste(plot_pf, '_pair_sub_pam50_os_bhadj.csv', sep = ''))

fit <- survfit(Surv(rfst, rfse) ~ pam50_TAM, data = pam_df)
surv_gg <- ggsurvplot(fit, data = pam_df, pval = TRUE,
                      title = paste("PD-L1", pdl1_status, "-", co),
                      legend.title = 'TAM',
                      risk.table = TRUE)
png(paste(plot_pf, '_pam50_rfs_km.png', sep = '_'), res = 600, width = 7.5, height = 6, units = 'in')
print(surv_gg)
gar <- dev.off()
pres <- pairwise_survdiff(Surv(rfst, rfse) ~ pam50_TAM, data=pam_df)
write.csv(pres$p.value, paste(plot_pf, '_pair_sub_pam50_rfs_bhadj.csv', sep = ''))
