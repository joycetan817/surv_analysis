# For Colt's NanoString Data
# Prefer to MS OS
# Weihua Guo. Ph.D.
# 07/22/2019

rm(list = ls(all.names = TRUE))

library(readxl)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(corrplot)

data_dir = "//bri-net/citi/Peter Lee Group/Weihua/PD1_CD39_nanostring/pure_data"
plot_dir = "//bri-net/citi/Peter Lee Group/Weihua/PD1_CD39_nanostring/results"
plot_dir = "~/temp_work_athome"
expr_file = "ns_expr_log2.xlsx" # Expression
annot_file = "fold_ihc_clinical_annot.xlsx" # annotation
cscore_file = "ns_cell_score_raw_log2.xlsx" # Raw cell score
# cscore_file = "ns_cell_score_relative_log2.xlsx" # Relative cell score
merge_out = "ns_expr_annot_cellscore_merge"

scale_flag = FALSE
cell_score_hmplot = FALSE
expr_hmplot = FALSE
corr_flag = TRUE

raw.expr = as.data.frame(read_excel(paste(data_dir, "/", expr_file, sep = ""), ))
rownames(raw.expr) = raw.expr$pid
expr = raw.expr
expr$pid = NULL

raw.annot = as.data.frame(read_excel(paste(data_dir, "/", annot_file, sep = ""), ))
annot = raw.annot
rownames(annot) = annot$pid
annot$pid = NULL

raw.cscore = as.data.frame(read_excel(paste(data_dir, "/", cscore_file, sep = ""), ))
cscore = raw.cscore
rownames(cscore) = cscore$pid
cscore$pid = NULL

raw.merge = merge(raw.annot, raw.expr, by = "pid")
raw.merge = merge(raw.merge, raw.cscore, by = "pid")

# saveRDS(raw.merge, file = paste(data_dir, "/", merge_out, ".RDS", sep = ""))
# write.csv(raw.merge, file = paste(data_dir, "/", merge_out, ".csv", sep = ""))

## Start to draw the heatmap
if (cell_score_hmplot) {
# For cell score
if (scale_flag) {
	pl_cscore = t(scale(cscore))
	valname = "Z-score"
	cscore_hm_tiff = "cell_score_raw_log2_zscore_heatmap_v1.tiff"
	cscore_clin_hm_tiff = "cell_score_raw_log2_zscore_clin_heatmap_v1.tiff"
} else {
	pl_cscore = t(cscore)
	valname = "Cell score (log2)"
	cscore_hm_tiff = "cell_score_raw_log2_heatmap_v1.tiff"
	cscore_clin_hm_tiff = "cell_score_raw_log2_clin_heatmap_v1.tiff"
}

str_ci_color = c("#FFFF66","#66FFFF") # Stroma/CI color code
phen_ann = HeatmapAnnotation(
	pdl1_ihc = anno_barplot(annot[,c("CD274_stroma", "CD274_epi")],
		gp = gpar(fill = str_ci_color), height = unit(1.5, "cm")),
	cd8_ihc = anno_barplot(annot[,c("CD8_stroma", "CD8_epi")],
		gp = gpar(fill = str_ci_color), height = unit(1.5, "cm")),
	cd20_ihc = anno_barplot(annot[,c("CD20_stroma", "CD20_epi")],
		gp = gpar(fill = str_ci_color), height = unit(1.5, "cm")),
	Tex = anno_barplot(annot$Tex_FACS,
		gp = gpar(fill = "purple"), height = unit(1.2, "cm")),
	gap = unit(0.18, "cm")
)

ldg = Legend(labels = c("Stroma", "Cancer Island"), 
	title = "IHC", legend_gp = gpar(fill = str_ci_color))

names(phen_ann) = c("PD-L1 (IHC)", "CD8 (IHC)", "CD20 (IHC)", "PD-L1+CD39+\nTex (FACS%)")

cshm = Heatmap(pl_cscore, name = valname,
	show_column_name = FALSE,
	column_order = order(annot$Tex_FACS),
	top_annotation = phen_ann,
	heatmap_legend_param = list(direction = "horizontal"))

cshm_save = paste(plot_dir, "/", cscore_hm_tiff, sep = "")
tiff(cshm_save, res = 180, width = 9, heigh = 9, units = 'in')
draw(cshm, annotation_legend_list = ldg, merge_legend = TRUE,
	heatmap_legend_side = "bottom", 
	annotation_legend_side = "bottom")
gar = dev.off()

clin_ann = annot[,c("KI67","grade","stage")]
clin_ann$KI67 = factor(clin_ann$KI67)
clin_ann$stage = factor(clin_ann$stage)
clin_ann$grade = factor(clin_ann$grade)
stage_col = colorRampPalette(brewer.pal(6,"YlOrRd"))(length(levels(clin_ann$stage)))
names(stage_col) = levels(clin_ann$stage)
grade_col = colorRampPalette(brewer.pal(9,"Blues"))(length(levels(clin_ann$stage)))
names(grade_col) = levels(clin_ann$stage)
ki67_col = colorRampPalette(brewer.pal(9,"PuOr"))(length(levels(clin_ann$KI67)))
names(ki67_col) = levels(clin_ann$KI67)

annot_lgd_set = list(labels_gp = gpar(fontsize = 15), 
	title_gp = gpar(fontsize = 18), grid_height = unit(0.81,"cm"))

ha = HeatmapAnnotation(df = clin_ann,
	col = list(KI67 = ki67_col, grade = grade_col, stage = stage_col))
names(ha) = c("Ki67", "Grade", "Stage")
cshm_clin = Heatmap(pl_cscore, name = valname,
	show_column_name = FALSE,
	column_order = order(annot$Tex_FACS),
	top_annotation = ha,
	heatmap_legend_param = list(direction = "horizontal"))

cshm_save = paste(plot_dir, "/", cscore_clin_hm_tiff, sep = "")
tiff(cshm_save, res = 180, width = 9, heigh = 9, units = 'in')
draw(cshm_clin, merge_legend = TRUE,
	heatmap_legend_side = "bottom",
	annotation_legend_side = "bottom")
gar = dev.off()
}

if (expr_hmplot) {
# For expression
if (scale_flag) {
	pl_expr = t(scale(expr))
	valname = "Z-score"
	expr_hm_tiff = "expr_norm_log2_zscore_heatmap_v1.tiff"
	expr_clin_hm_tiff = "expr_norm_log2_zscore_clin_heatmap_v1.tiff"
} else {
	pl_expr = t(expr)
	valname = "Expression (log2)"
	expr_hm_tiff = "expr_norm_log2_heatmap_v1.tiff"
	expr_clin_hm_tiff = "expr_norm_log2_clin_heatmap_v1.tiff"
}

str_ci_color = c("#FFFF66","#66FFFF") # Stroma/CI color code
phen_ann = HeatmapAnnotation(
	pdl1_ihc = anno_barplot(annot[,c("CD274_stroma", "CD274_epi")],
		gp = gpar(fill = str_ci_color), height = unit(1.5, "cm")),
	cd8_ihc = anno_barplot(annot[,c("CD8_stroma", "CD8_epi")],
		gp = gpar(fill = str_ci_color), height = unit(1.5, "cm")),
	cd20_ihc = anno_barplot(annot[,c("CD20_stroma", "CD20_epi")],
		gp = gpar(fill = str_ci_color), height = unit(1.5, "cm")),
	Tex = anno_barplot(annot$Tex_FACS,
		gp = gpar(fill = "purple"), height = unit(1.2, "cm")),
	gap = unit(0.18, "cm")
)

ldg = Legend(labels = c("Stroma", "Cancer Island"), 
	title = "IHC", legend_gp = gpar(fill = str_ci_color))

names(phen_ann) = c("PD-L1 (IHC)", "CD8 (IHC)", "CD20 (IHC)", "PD-L1+CD39+\nTex (FACS%)")

ephm = Heatmap(pl_expr, name = valname,
	show_column_name = FALSE,
	show_row_name = FALSE,
	column_order = order(annot$Tex_FACS),
	top_annotation = phen_ann,
	heatmap_legend_param = list(direction = "horizontal"),
	width = unit(9, "cm"), height = unit(9, "cm"))

ephm_save = paste(plot_dir, "/", expr_hm_tiff, sep = "")
tiff(ephm_save, res = 180, width = 6, heigh = 7, units = 'in')
draw(ephm, annotation_legend_list = ldg, merge_legend = TRUE,
	heatmap_legend_side = "bottom", 
	annotation_legend_side = "bottom")
gar = dev.off()

clin_ann = annot[,c("KI67","grade","stage")]
clin_ann$KI67 = factor(clin_ann$KI67)
clin_ann$stage = factor(clin_ann$stage)
clin_ann$grade = factor(clin_ann$grade)
stage_col = colorRampPalette(brewer.pal(6,"YlOrRd"))(length(levels(clin_ann$stage)))
names(stage_col) = levels(clin_ann$stage)
grade_col = colorRampPalette(brewer.pal(9,"Blues"))(length(levels(clin_ann$stage)))
names(grade_col) = levels(clin_ann$stage)
ki67_col = colorRampPalette(brewer.pal(9,"PuOr"))(length(levels(clin_ann$KI67)))
names(ki67_col) = levels(clin_ann$KI67)

annot_lgd_set = list(labels_gp = gpar(fontsize = 15), 
	title_gp = gpar(fontsize = 18), grid_height = unit(0.81,"cm"))

ha = HeatmapAnnotation(df = clin_ann,
	col = list(KI67 = ki67_col, grade = grade_col, stage = stage_col))
names(ha) = c("Ki67", "Grade", "Stage")
ephm_clin = Heatmap(pl_expr, name = valname,
	show_column_name = FALSE,
	show_row_name = FALSE,
	column_order = order(annot$Tex_FACS),
	top_annotation = ha,
	heatmap_legend_param = list(direction = "horizontal"),
	width = unit(9, "cm"), height = unit(9, "cm"))

ephm_save = paste(plot_dir, "/", expr_clin_hm_tiff, sep = "")
tiff(ephm_save, res = 180, width = 6, heigh = 6, units = 'in')
draw(ephm_clin, merge_legend = TRUE,
	heatmap_legend_side = "bottom",
	annotation_legend_side = "bottom")
gar = dev.off()
}

if (corr_flag) {
	# Corr2CellScore
	cs_annot = merge(raw.cscore, raw.annot, by = "pid")
	rownames(cs_annot) = cs_annot$pid
	cs_annot$pid = NULL
	cs_pheno = cs_annot[,1:(ncol(cs_annot)-7)]
	print(class(cs_pheno))
	cs_pheno = as.matrix(cs_pheno)
	cs_pheno_corres = rcorr(cs_pheno)
	csph_r = cs_pheno_corres[['r']]
	csph_p = cs_pheno_corres[['P']]
	csr_csv = "correlation_cell_score_2_ihc_facs_r.csv"
	csp_csv = "correlation_cell_score_2_ihc_facs_pval.csv"
	write.csv(csph_r, file = paste(plot_dir, "/", csr_csv, sep = ""))
	write.csv(csph_p, file = paste(plot_dir, "/", csp_csv, sep = ""))


	colscale <- colorRampPalette(c("blue", "white", "red")) 
	csph_cpsave = "corrplot_cell_score_ihc_facs_sign.tiff"
	tiff(paste(plot_dir, "/", csph_cpsave, sep = ""), 
		res = 180, width = 9, heigh = 9, units = 'in')
	corrplot(csph_r, p.mat = csph_p, sig.level = 0.01, insig = "blank",
		col = colscale(100), tl.col = "black")
	gar = dev.off()
}