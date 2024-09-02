rm(list = ls())

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(rstatix))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(circlize))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(forestmodel))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggwordcloud))

data_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/brinet_Joyce_backup_240814/TCGA_pancancer"
score_file <- "all_gene_signature_pancancer_chrom_instab_update1017.csv"
timer_file <- "infiltration_estimation_for_tcga.csv"
ksdr_dir <- "Kassandra_TCGA"

expr_id <- "hypoxia_tex_corr_240902"

dir.create(file.path(data_dir, expr_id), showWarnings = FALSE)
expDir <- paste(data_dir, expr_id, sep = "/")
ppf <- paste(data_dir, "/", expr_id, "/", expr_id, "_", sep = "")
set.seed(669)


ksdr_merge_flag <- FALSE
b_tex_corr_flag <- FALSE
b_tex_corr_vis_flag <- TRUE
ratio_surv_flag <- FALSE

if (ksdr_merge_flag) {
	all_files <- list.files(paste(data_dir, ksdr_dir, sep = "/"))
	i <- 0
	for (ik in all_files) {
		tmp_df <- read.table(paste(data_dir, ksdr_dir, ik, sep = "/"), sep = "\t", header = T, check.names = F, row.names = 1)
		if (i == 0) {
			ksdr_df <- tmp_df
		} else {
			ksdr_df <- cbind(ksdr_df, tmp_df)
		}
		i <- i+1
	}
	write.csv(ksdr_df,paste(data_dir, "all_kassadra_dataframe.csv", sep = "/"))
	print(dim(ksdr_df))
	id_df <- as.data.frame(str_split_fixed(colnames(ksdr_df), "-", n = 4))
	id_df$col <- str_c(id_df$V1, "-", id_df$V2, "-", id_df$V3, "-", id_df$V4)
	id_df$new_col <- str_c(id_df$V1, "-", id_df$V2, "-", id_df$V3)
	id_df <- id_df[id_df$V4 == "01",]
	pt_ksdr_df <- ksdr_df[,id_df$col]
	colnames(pt_ksdr_df) <- id_df$new_col
	write.csv(pt_ksdr_df,paste(data_dir, "pt_kassadra_dataframe.csv", sep = "/"))
}

if (ratio_surv_flag) {
	score_df <- read.csv(paste(data_dir, score_file, sep = "/"), sep = ",")
	score_df$ID <- str_replace_all(score_df$pid, "\\.", "-")
	score_df$tex_cd8a_ratio <- score_df$Tex/score_df$CD8A
	score_df$tex_cd8a_ratio_group <- ifelse(score_df$tex_cd8a_ratio > median(score_df$tex_cd8a_ratio), "High", "Low")
	print(head(score_df))
	tmp_lfit <- survfit(Surv(ost, ose) ~ tex_cd8a_ratio_group, data = score_df, id = ID)
	tmp_lgg <- ggsurvplot(tmp_lfit, pval = T, ggtheme = theme_bw(), legend = "right") + 
		labs(title = "All cancer types", x = "Time (day)", color = "Tex/CD8A ratio (median)") 		
	png(paste(ppf, 'all_cancer_tex_cd8a_ratio_median_kmplot.png', sep = ""), res = 300, width = 9, height = 4, units = 'in')
	print(tmp_lgg)
	gar <- dev.off()
	tmp_cox <- coxph(Surv(ost, ose) ~ tex_cd8a_ratio_group, data=score_df, id=ID, ties="breslow")
	write.csv(summary(tmp_cox)$coefficients, paste(ppf, 'all_cancer_tex_cd8a_ratio_median_cox.csv', sep = ""))

	png(paste(ppf, 'all_cancer_tex_cd8a_ratio_median_forestmodel_hr.png', sep = ""), res = 300, width = 9, height = 4, units = 'in')
	print(forest_model(coxph(Surv(ost, ose) ~ tex_cd8a_ratio_group, data=score_df))+labs(title = "All cancer types"))
	gar <- dev.off()
	png(paste(ppf, 'all_cancer_tex_cd8a_ratio_median_forestmodel_ratio.png', sep = ""), res = 300, width = 9, height = 4, units = 'in')
	print(forest_model(coxph(Surv(ost, ose) ~ tex_cd8a_ratio, data=score_df))+labs(title = "All cancer types"))
	gar <- dev.off()


	for (ic in unique(score_df$cancertype)) {
		tmp_df <- score_df[score_df$cancertype == ic,]
		tmp_df$tex_cd8a_ratio <- tmp_df$Tex/tmp_df$CD8A
		tmp_df$tex_cd8a_ratio_group <- ifelse(tmp_df$tex_cd8a_ratio > median(tmp_df$tex_cd8a_ratio), "High", "Low")
		print(head(tmp_df))
		tmp_lfit <- survfit(Surv(ost, ose) ~ tex_cd8a_ratio_group, data = tmp_df, id = ID)
		tmp_lgg <- ggsurvplot(tmp_lfit, pval = T, ggtheme = theme_bw(), legend = "right") + 
			labs(title = ic, x = "Time (day)", color = "Tex/CD8A ratio (median)") 		
		png(paste(ppf, ic, '_tex_cd8a_ratio_median_kmplot.png', sep = ""), res = 300, width = 9, height = 4, units = 'in')
		print(tmp_lgg)
		gar <- dev.off()
		tmp_cox <- coxph(Surv(ost, ose) ~ tex_cd8a_ratio_group, data=tmp_df, id=ID, ties="breslow")
		write.csv(summary(tmp_cox)$coefficients, paste(ppf, ic, '_tex_cd8a_ratio_median_cox.csv', sep = ""))

		png(paste(ppf, ic, '_tex_cd8a_ratio_median_forestmodel_hr.png', sep = ""), res = 300, width = 9, height = 4, units = 'in')
		print(forest_model(coxph(Surv(ost, ose) ~ tex_cd8a_ratio_group, data=tmp_df))+labs(title = ic))
		gar <- dev.off()
		png(paste(ppf, ic, '_tex_cd8a_ratio_median_forestmodel_ratio.png', sep = ""), res = 300, width = 9, height = 4, units = 'in')
		print(forest_model(coxph(Surv(ost, ose) ~ tex_cd8a_ratio, data=tmp_df))+labs(title = ic))
		gar <- dev.off()
	}


}

if (b_tex_corr_flag) {
	score_df <- read.csv(paste(data_dir, score_file, sep = "/"), sep = ",")
	score_df$ID <- str_replace_all(score_df$pid, "\\.", "-")
	print(dim(score_df))
	timer_df <- read.csv(paste(data_dir, timer_file, sep = "/"), sep = ",", check.names = F, row.names = 1)
	timer_df$`B cell_CIBERSORT` <- rowSums(timer_df[,c("B cell memory_CIBERSORT", "B cell naive_CIBERSORT", "B cell plasma_CIBERSORT")])
	timer_df$`B cell_CIBERSORT-ABS` <- rowSums(timer_df[,c("B cell memory_CIBERSORT-ABS", "B cell naive_CIBERSORT-ABS", "B cell plasma_CIBERSORT-ABS")])
	timer_cols <- colnames(timer_df)
	timer_df$SID <- str_sub(rownames(timer_df), -2,-1)
	timer_df <- timer_df[timer_df$SID == "01",]
	timer_df$ID <- str_sub(rownames(timer_df), 1, -4)
	print(dim(timer_df))
	ksdr_df <- as.data.frame(t(read.csv(paste(data_dir, "pt_kassadra_dataframe.csv", sep = "/"), sep = ",", check.names = F, row.names = 1)))
	ksdr_df <- ksdr_df[,!str_detect(colnames(ksdr_df), "_std")]
	colnames(ksdr_df) <- str_c(str_replace_all(colnames(ksdr_df), "_", " "), "_Kassandra")
	ksdr_cols <- colnames(ksdr_df)
	ksdr_df$ID <- rownames(ksdr_df)
	print(dim(ksdr_df))
	print(colnames(ksdr_df))

	merge_df <- merge(score_df, timer_df, by = "ID", all.x = T)
	merge_df <- merge(merge_df, ksdr_df, by = "ID", all.x = T)
	write.csv(merge_df, paste(data_dir, "all_scores.csv", sep = "/"))
	print(dim(merge_df))
	print(colnames(merge_df))
	print(table(merge_df$cancertype))
	if (FALSE) { #CORR
	cust_cols <- c("Tex", "PDCD1", "CD274", "TLS", "CD8A", "TMB", "hypoxia")
	cor_res <- merge_df %>%
		group_by(cancertype) %>%
		cor_test(vars = cust_cols, vars2 = c(timer_cols, ksdr_cols, "TLS", "hypoxia"), method = "spearman") %>%
		add_significance('p')
	write.csv(cor_res, paste(data_dir, "spearman_correlation_pre_cancer_all_scores_result.csv", sep = "/"))
	cor_res <- merge_df %>%
		group_by(cancertype) %>%
		cor_test(vars = cust_cols, vars2 = c(timer_cols, ksdr_cols, "TLS", "hypoxia"), method = "pearson")  %>%
		add_significance('p')
	write.csv(cor_res, paste(data_dir, "pearson_correlation_per_cancer_all_scores_result.csv", sep = "/"))

	cor_res <- merge_df %>%
		cor_test(vars = cust_cols, vars2 = c(timer_cols, ksdr_cols, "TLS", "hypoxia"), method = "spearman") %>%
		add_significance('p')
	write.csv(cor_res, paste(data_dir, "spearman_correlation_all_scores_result.csv", sep = "/"))
	cor_res <- merge_df %>%
		cor_test(vars = cust_cols, vars2 = c(timer_cols, ksdr_cols, "TLS", "hypoxia"), method = "pearson") %>%
		add_significance('p')
	write.csv(cor_res, paste(data_dir, "pearson_correlation_all_scores_result.csv", sep = "/"))
	} # CORR

	gath_df <- gather(merge_df, "celltype", "value", c(timer_cols, ksdr_cols, "TLS", "hypoxia"))
	gath_df$keep_id <- str_c(gath_df$cancertype, "+", gath_df$celltype)
	stat_df <- gath_df %>%
		group_by(cancertype, celltype, group) %>%
		summarize(n = n(), 
			  avg = mean(value),
			  median = median(value),
			  min = min(value),
			  max = max(value),
			  sd = sd(value))
	stat_df <- stat_df %>%
		group_by(cancertype, celltype) %>%
		mutate(n_zero_sd = sum(sd == 0))
	stat_df$keep_id <- str_c(stat_df$cancertype, "+", stat_df$celltype)
	keep_ids <- stat_df$keep_id[stat_df$n_zero_sd == 0]
	print(head(stat_df))
	write.csv(stat_df, paste(data_dir, "tex_group_cell_scores_stat_df.csv", sep = "/"))

	tex_pw_res <- gath_df %>% 
		filter(keep_id %in% keep_ids) %>%
		group_by(cancertype, celltype) %>%
		pairwise_wilcox_test(value ~ group, detailed = T) %>%
		add_significance('p') %>%
		adjust_pvalue('p') %>%
		add_significance('p.adj')
	write.csv(tex_pw_res, paste(data_dir, "tex_group_one2one_wilcox_result.csv", sep = "/"))
}

if (b_tex_corr_vis_flag) {
	cor_res <- read.csv(paste(data_dir, "pearson_correlation_all_scores_result.csv", sep = "/"), row.names = 1)
	cor_res$tool <- str_split_fixed(cor_res$var2, "_", n = 2)[,2]
	cor_res$tool[cor_res$tool == ""] <- "Cust"
	cor_res$cell <- str_split_fixed(cor_res$var2, "_", n = 2)[,1]
	cor_res$clean_cell <- cor_res$cell
	cor_res$clean_cell <- str_replace_all(cor_res$clean_cell, "cells", "cell")
	cor_res$clean_cell <- str_replace_all(cor_res$clean_cell, "s ", " ")
	cor_res$clean_cell <- str_replace_all(cor_res$clean_cell, "(.*)s$", "\\1")
	cor_res$clean_cell[str_detect(cor_res$cell, "ibrobl")] <- "Cancer associated fibroblast"
	cor_res$clean_cell[str_detect(cor_res$cell, "endritic")] <- "Dendritic cell"
	cor_res$clean_cell[str_detect(cor_res$cell, "NK cell")] <- "NK cell"
	cor_res$clean_cell[str_detect(cor_res$cell, "Mast")]  <- "Mast cell"
	cor_res$clean_cell[str_detect(cor_res$cell, "Treg")]  <- "Treg"
	cor_res$clean_cell[str_detect(cor_res$cell, "CD4 T")]  <- "T cell CD4+"
	cor_res$clean_cell[str_detect(cor_res$cell, "dothel")] <- "Endothelial cell"
	cor_res$clean_cell[cor_res$cell == "CD8 T cells"]  <- "T cell CD8+"

	print(head(cor_res))
	print(unique(cor_res$clean_cell))


	tex_res <- cor_res[cor_res$var1 == "hypoxia",]
	sum_res <- tex_res %>%
		group_by(clean_cell) %>%
		summarize(avg_cor = mean(cor), 
			  n_tool = n(), 
			  median_cor = median(cor)) %>%
		filter(n_tool >= 3)
	print(sum_res)
	sum_res$abs_avg_cor <- abs(sum_res$avg_cor)
	sum_res$cor_direction <- ifelse(sum_res$avg_cor<0, "Negative", "Positive")

	gg <- ggplot(sum_res, aes(label = clean_cell, size = abs_avg_cor, color = avg_cor)) +
		geom_text_wordcloud(shape = "pentagon") +
		scale_size_area(max_size = 20) +
		scale_color_gradient2(midpoint = 0, high = "firebrick", low = "dodgerblue") +
		theme_minimal() +
		theme(legend.position = "right")
	ggsave(paste(ppf, "pearson_hypoxia_clean_cell_cloud.png", sep = ""), dpi = 300, width = 12, height = 9)

	for (it in unique(cor_res$tool)) {
		cat(it, "\n")
		tmp_res <- cor_res[cor_res$tool == it,]
		tmp_res <- tmp_res[tmp_res$var1 == "hypoxia",]
#		print(head(tmp_res))
#		print(dim(tmp_res))

		dot_gg <- ggplot(tmp_res, aes(x = reorder(cell, cor), y = cor)) +
		geom_point(aes(color = cor, size = p.signif)) +
		scale_size_manual(values = c("ns" = 0.1, "*" = 2, "**" = 4, "***" = 6, "****" = 8)) +
		scale_color_gradient2(midpoint = 0, high = "firebrick", low = "dodgerblue") +
		labs(x = paste("Cell types in", it), color = "PCC with\nhypoxia signature\nscores", size = "Statistic\nsignificance", y = "Tex") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1),
		      strip.text.x = element_text(angle = 45), legend.position = "right"
		)
		ggsave(paste(ppf, it, "_pearson_hypoxia_vs_b_dot_all_score.png", sep = ""), dpi = 300, width = 1.8+0.5*nrow(tmp_res), height = 6)

	}

	if (FALSE) {#B2
	tex_b_res <- cor_res[cor_res$var1 == "Tex",]
	tex_b_res <- tex_b_res[str_detect(tex_b_res$var2, "B cell"),]

	print(head(tex_b_res))

	dot_gg <- ggplot(tex_b_res, aes(x = cell, y = cor)) +
		geom_point(aes(color = cor, size = p.signif)) +
		scale_size_manual(values = c("ns" = 0.1, "*" = 2, "**" = 4, "***" = 6, "****" = 8)) +
		scale_color_gradient2(midpoint = 0, high = "firebrick", low = "dodgerblue") +
		facet_grid(.~tool, space = "free", scales = "free") +
		labs(x = "Types of B cells", color = "Pearson correlation coefficients with\nTex signature scores", size = "Statistic significance", y = "Tex") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1),
		      strip.text.x = element_text(angle = 45), legend.position = "top"
		)
	ggsave(paste(ppf, "pearson_hypoxia_vs_b_dot_all_score.png", sep = ""), dpi = 300, width = 12, height = 6)
	}#B2

	cor_res <- read.csv(paste(data_dir, "pearson_correlation_per_cancer_all_scores_result.csv", sep = "/"), row.names = 1)
	cor_res$tool <- str_split_fixed(cor_res$var2, "_", n = 2)[,2]
	cor_res$tool[cor_res$tool == ""] <- "Cust"
	cor_res$cell <- str_split_fixed(cor_res$var2, "_", n = 2)[,1]
	print(head(cor_res))
	cor_res$clean_cell <- cor_res$cell
	cor_res$clean_cell <- str_replace_all(cor_res$clean_cell, "cells", "cell")
	cor_res$clean_cell <- str_replace_all(cor_res$clean_cell, "s ", " ")
	cor_res$clean_cell <- str_replace_all(cor_res$clean_cell, "(.*)s$", "\\1")
	cor_res$clean_cell[str_detect(cor_res$cell, "ibrobl")] <- "Cancer associated fibroblast"
	cor_res$clean_cell[str_detect(cor_res$cell, "endritic")] <- "Dendritic cell"
	cor_res$clean_cell[str_detect(cor_res$cell, "NK cell")] <- "NK cell"
	cor_res$clean_cell[str_detect(cor_res$cell, "Mast")]  <- "Mast cell"
	cor_res$clean_cell[str_detect(cor_res$cell, "Treg")]  <- "Treg"
	cor_res$clean_cell[str_detect(cor_res$cell, "CD4 T")]  <- "T cell CD4+"
	cor_res$clean_cell[str_detect(cor_res$cell, "dothel")] <- "Endothelial cell"
	cor_res$clean_cell[cor_res$cell == "CD8 T cells"]  <- "T cell CD8+"

	for (ic in unique(cor_res$cancertype)) {
		cat(ic, "\n")
		tex_res <- cor_res[cor_res$var1 == "hypoxia",]
		tex_res <- tex_res[tex_res$cancertype == ic,]
		sum_res <- tex_res %>%
			group_by(clean_cell) %>%
			summarize(avg_cor = mean(cor), 
				  n_tool = n(), 
				  median_cor = median(cor)) %>%
			filter(n_tool >= 3)
		print(sum_res)
		sum_res$abs_avg_cor <- abs(sum_res$avg_cor)
		sum_res$cor_direction <- ifelse(sum_res$avg_cor<0, "Negative", "Positive")

		gg <- ggplot(sum_res, aes(label = clean_cell, size = abs_avg_cor, color = avg_cor)) +
			geom_text_wordcloud(shape = "pentagon") +
			scale_size_area(max_size = 20) +
			scale_color_gradient2(midpoint = 0, high = "firebrick", low = "dodgerblue") +
			labs(title = ic) +
			theme_minimal() +
			theme(legend.position = "right")
		ggsave(paste(ppf, ic, "_pearson_hypoxia_clean_cell_cloud.png", sep = ""), dpi = 300, width = 12, height = 9)
	}

	for (it in unique(cor_res$tool)) {
		cat(it, "\n")
		tmp_res <- cor_res[cor_res$tool == it,]
		tmp_res <- tmp_res[tmp_res$var1 == "hypoxia",]
#		print(head(tmp_res))
#		print(dim(tmp_res))

		dot_gg <- ggplot(tmp_res, aes(x = reorder(cell, cor), y = cancertype)) +
		geom_point(aes(color = cor, size = p.signif)) +
		scale_size_manual(values = c("ns" = 0.1, "*" = 2, "**" = 4, "***" = 6, "****" = 8)) +
		scale_color_gradient2(midpoint = 0, high = "firebrick", low = "dodgerblue") +
		labs(x = paste("Cell types in", it), color = "PCC with\nhypoxia signature\nscores", size = "Statistic\nsignificance", y = "Cancer types") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1),
		      strip.text.x = element_text(angle = 45), legend.position = "right"
		)
		ggsave(paste(ppf, it, "_pearson_hypoxia_vs_b_dot_per_cancertype_score.png", sep = ""), dpi = 300, width = 1.8+0.5*length(unique(tmp_res$cell)), height = 6)


		cor_df <- spread(tmp_res[,c("cancertype", "var2", "cor")], "var2", "cor")
		col_max <- max(c(max(tmp_res$cor, na.rm = T), abs(min(tmp_res$cor, na.rm = T))))
		rownames(cor_df) <- cor_df$cancertype
		cor_df$cancertype <- NULL
		col_df <- as.data.frame(str_split_fixed(colnames(cor_df), "_", n = 2))
		colnames(col_df) <- c("cell", "tool")
		color_fun = colorRamp2(c(-col_max, 0, col_max), c("dodgerblue", "white", "firebrick"))
		colnames(cor_df) <- col_df$cell

		hm <- Heatmap(cor_df, name = "PCC to hypoxia", col = color_fun, rect_gp = gpar(col = "white", lwd = 2))
		png(paste(ppf, it, "_pearson_hypoxia_vs_b_heatmap_per_cancertype_score.png", sep = ""), res = 300, width = 1.8+0.5*length(unique(tmp_res$cell)), height = 7, units = 'in')
		print(draw(hm))
		gar <- dev.off()
		print(dim(cor_df))
	}

	q(save = "no")

	if (FALSE) {#B1
	tex_b_res <- cor_res[cor_res$var1 == "Tex",]
	tex_b_res <- tex_b_res[str_detect(tex_b_res$var2, "B cell"),]
	print(head(tex_b_res))

	dot_gg <- ggplot(tex_b_res, aes(x = cell, y = cancertype)) +
		geom_point(aes(color = cor, size = p.signif)) +
		scale_size_manual(values = c("ns" = 0.1, "*" = 2, "**" = 4, "***" = 6, "****" = 8)) +
		scale_color_gradient2(midpoint = 0, high = "firebrick", low = "dodgerblue") +
		facet_grid(.~tool, space = "free", scales = "free") +
		labs(x = "Types of B cells", color = "Pearson correlation coefficients with\nTex signature scores", size = "Statistic significance", y = "Cancer types") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1),
		      strip.text.x = element_text(angle = 45), legend.position = "top"
		)
	ggsave(paste(ppf, "pearson_hypoxia_vs_b_dot_per_cancertype_score.png", sep = ""), dpi = 300, width = 12, height = 7)

	cor_df <- spread(tex_b_res[,c("cancertype", "var2", "cor")], "var2", "cor")
	col_max <- max(c(max(tex_b_res$cor), abs(min(tex_b_res$cor))))
	rownames(cor_df) <- cor_df$cancertype
	cor_df$cancertype <- NULL
	col_df <- as.data.frame(str_split_fixed(colnames(cor_df), "_", n = 2))
	colnames(col_df) <- c("B_cell_type", "tool")
	color_fun = colorRamp2(c(-col_max, 0, col_max), c("dodgerblue", "white", "firebrick"))

	hm <- Heatmap(cor_df, name = "PCC to Tex", col = color_fun, column_split = col_df$tool)
	png(paste(ppf, "pearson_hypoxia_vs_b_heatmap_per_cancertype_score.png", sep = ""), res = 300, width = 12, height = 7, units = 'in')
	print(draw(hm))
	gar <- dev.off()
	print(dim(cor_df))
	}#B1
}
