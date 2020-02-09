# This script is for scRNA-seq analysis using Seurat package
# Jiayi Tan
# 01/22/2020

rm(list = ls(all.names = TRUE))

batch_corr = function (sc_obj, group) {
	sample.list <- SplitObject(sc_obj, split.by = group)
	for (i in 1:length(sample.list)) {
    sample.list[[i]] <- FindVariableFeatures(sample.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
	}

	reference.list <- sample.list
	integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
	sc_cor <- IntegrateData(anchorset = integrate.anchors, dims = 1:30)
	DefaultAssay(sc_cor) <- "integrated"
	return (sc_cor)
}

sr_clst = function(sc_obj, ndim, resdir) {
	sc_obj <- ScaleData(sc_obj)
	sc_obj <- RunPCA(sc_obj)
	pca_vizdim = VizDimLoadings(sc_obj, dim = 1:6, reduction = "pca")
	ggsave(plot = pca_vizdim, filename = paste(resdir, "pca_vizdim_plot.tiff", sep = ""), width = 10, height = 15)
	
	pca_dim = DimPlot(sc_obj, reduction = "pca")
	ggsave(plot = pca_dim, filename = paste(resdir, "pca_dim_plot.tiff", sep = ""), width = 9, height = 6)

	tiff(paste(resdir, "pca_heatmap.tiff", sep = ""), res = 180, width = 12, height = 16, units = "in")
	DimHeatmap(sc_obj, dims = 1:15, cells = 500, balanced = TRUE)
	gar = dev.off()

	sc_obj = JackStraw(sc_obj, num.replicate = 100)
	sc_obj = ScoreJackStraw(sc_obj, dims = 1:20)
	js_plot = JackStrawPlot(sc_obj, dims = 1:20)
	ggsave(plot = js_plot, filename = paste(resdir, "JackStraw_plot.tiff", sep = ""), width = 9, heigh = 6)

	elb_plot = ElbowPlot(sc_obj, ndims = 30)
	ggsave(plot = elb_plot, filename = paste(resdir, "elbow_plot.tiff", sep = ""), width = 9, heigh = 6)

	saveRDS(sc_obj, paste(resdir, "seurat_obj_cor_bf_pc.RDS", sep = ""))

	sc_obj = FindNeighbors(sc_obj, dims = 1:ndim)
	sc_obj = FindClusters(sc_obj, resolution = 0.5)

	sc_obj = RunUMAP(sc_obj, dims = 1:ndim)
	cls_umap = DimPlot(sc_obj, reduction = "umap")
	ggsave(plot = cls_umap, filename = paste(resdir, ndim, "_cluster_umap_plot.tiff", sep = ""), width = 9, heigh = 6)
	btc = DimPlot(sc_obj, reduction = "umap", group.by = "orig.ident")
	ggsave(plot = btc, filename = paste(resdir, "pca", ndim, "_batch_check_aft_cor.tiff", sep = ""), width = 9, heigh = 6)


	expr = c("PTPRC", "CD3E", "CD4", "CD8A", "TCF7", "IL7R", "ITGAE", "PDCD1", "ENTPD1", "DPP4", "FOXP3")
	expr_ol = FeaturePlot(sc_obj, features = expr, reduction = "umap",
					    min.cutoff = "q05", max.cutoff = "q95", ncol = 4, cols = c("orange1", "grey", "purple"))
	ggsave(plot = expr_ol, filename = paste(resdir, "expression_overlay.tiff", sep = ""),  width = 16, heigh = 10)

	cluster.markers <- FindAllMarkers(sc_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	write.csv(cluster.markers, paste(resdir, "pca", ndim, "_cluster_marker.csv", sep = ""))
	return (sc_obj)
}



library(dplyr)
library(Seurat)
library(ggplot2)

sc_obj.data <- Read10X(data.dir = "C:/Users/jitan/Documents/New folder/")
st = Sys.time()
sc_obj <- CreateSeuratObject(counts = rna_sc_obj, project = "sc_objER")
sc_obj
print(Sys.time()-st)

work_dir = "//bri-net/citi/Peter Lee Group/Joyce/single_cell/seurat_analysis/"
data_name = "Dana_tumor" #"Dana_TNBC" #Public scRNA raw data
cluster = "immune_res" # all immune cells 
reclst = "tcell_res" #recluster in t cells
data_dir = paste(work_dir, data_name, "/", sep = "") 
res_dir = paste(work_dir, data_name, cluster, "/", sep = "/")
reclst_res_dir = paste(work_dir, data_name, reclst, "/", sep = "/")



sc_obj <- readRDS(paste(work_dir, data_name, "BC3_5_TNBC_all_filter_seurat_project.RDS", sep= "/"))
sc_obj <- read.csv(paste(work_dir, data_name, "raw_corrected.csv", sep= "/"))
sc_obj <- CreateSeuratObject(counts = sc_obj, project = "sc_BC")
saveRDS(sc_obj, paste(res_dir, "all_filter_seurat_project.RDS", sep = ""))

batch_check = TRUE #if batch check is needed
batch_cor = TRUE #if batch correction is needed
ndim = 15 ## choose PC
recluster = TRUE

###########QC
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT.")
featplot <- VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot=featplot, file = paste(res_dir, "feature_plot_bf_filter.tiff", sep = ""),  width = 9, height = 6)
plot1 <- FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
featscat <- CombinePlots(plots = list(plot1, plot2))
ggsave(plot = featscat, file = paste(res_dir, "feature_scatter_bf_filter.tiff", sep = ""),  width = 9, height = 6)


sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 15) ###check feature plot
featplot_filt <- VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(plot = featplot_filt, file = paste(res_dir, "feature_plot_aft_filter.tiff", sep = ""),  width = 9, height = 6)



####Normalizing the data
sc_obj <- NormalizeData(sc_obj)

if (batch_check) {
	btc <- sc_obj
	btc <- FindVariableFeatures(btc, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(btc)
	btc <- ScaleData(btc, features = all.genes)
	btc <- RunPCA(btc, features = VariableFeatures(object = btc))
	btc <- RunUMAP(btc, reduction = "pca", dims = 1:30)
	bcplot <- DimPlot(btc, reduction = "umap", group.by = "orig.ident")
	ggsave(plot=bcplot, file = paste(res_dir, "batch_check_bf_cor.tiff", sep = ""),  width = 9, height = 6)
}


if (batch_cor) {
	sc_cor = batch_corr(sc_obj, group = "orig.ident")
	sc_cor = sc_clst(sc_cor, ndim = 15, resdir = res_dir)
} else {
	sc_obj = sc_clst(sc_obj, ndim = 15, resdir = res_dir)
}


if (recluster) {
	tc <- subset(sc_cor, ident = c(0,4,5,6,8)) ###recluster t cell based on CD3 and CD45
	sub_check <- DimPlot(tc)
	ggsave(plot = sub_check, filename = paste(reclst_res_dir, "subset_check.tiff", sep = ""),  width = 9, heigh = 6)
	saveRDS(tc, paste(reclst_res_dir, "seurat_obj_tcell_cluster.RDS", sep = ""))

if (batch_cor) {
	tc_cor = batch_corr(tc, group="orig.ident")
	tc_cor = sr_clst(tc_cor, ndim = 15, resdir = reclst_res_dir)
}


		tc_cor <- ScaleData(tc_cor)
		tc_cor <- RunPCA(tc_cor)
		pca_vizdim = VizDimLoadings(tc_cor, dim = 1:6, reduction = "pca")
		ggsave(plot = pca_vizdim, filename = paste(reclst_res_dir, "pca_vizdim_plot_cor.tiff", sep = ""), width = 10, height = 15)
		
		pca_dim = DimPlot(tc_cor, reduction = "pca")
		ggsave(plot = pca_dim, filename = paste(reclst_res_dir, "pca_dim_plot_cor.tiff", sep = ""), width = 9, height = 6)

		tiff(paste(reclst_res_dir, "pca_heatmap_cor.tiff", sep = ""), res = 180, width = 12, height = 16, units = "in")
		DimHeatmap(tc_cor, dims = 1:15, cells = 500, balanced = TRUE)
		gar = dev.off()

		tc_cor = JackStraw(tc_cor, num.replicate = 100)
		tc_cor = ScoreJackStraw(tc_cor, dims = 1:20)
		js_plot = JackStrawPlot(tc_cor, dims = 1:20)
		ggsave(plot = js_plot, filename = paste(reclst_res_dir, "JackStraw_plot_cor.tiff", sep = ""), width = 9, heigh = 6)

		elb_plot = ElbowPlot(tc_cor, ndims = 20)
		ggsave(plot = elb_plot, filename = paste(reclst_res_dir, "elbow_plot_cor.tiff", sep = ""), width = 9, heigh = 6)

		saveRDS(tc_cor, paste(reclst_res_dir, "seurat_obj_cor_bf_pc.RDS", sep = ""))

		tc_cor = FindNeighbors(tc_cor, dims = 1:ndim)
		tc_cor = FindClusters(tc_cor, resolution = 0.5)

		tc_cor = RunUMAP(tc_cor, dims = 1:ndim)
		cls_umap = DimPlot(tc_cor, reduction = "umap")
		ggsave(plot = cls_umap, filename = paste(reclst_res_dir, "pca", ndim, "_cluster_umap_plot_cor.tiff", sep = ""), width = 9, heigh = 6)
		btc = DimPlot(tc_cor, reduction = "umap", group.by = "orig.ident")
		ggsave(plot = btc, filename = paste(reclst_res_dir, ndim, "batch_check_aft_cor.tiff", sep = ""), width = 9, heigh = 6)


		expr = c("PTPRC", "CD3E", "CD4", "CD8A", "TCF7", "IL7R", "ITGAE", "PDCD1", "ENTPD1", "DPP4", "FOXP3")
		expr_ol = FeaturePlot(tc_cor, features = expr, reduction = "umap",
						    min.cutoff = "q05", max.cutoff = "q95", ncol = 4, cols = c("orange1", "grey", "purple"))
		ggsave(plot = expr_ol, filename = paste(reclst_res_dir, "expression_overlay_cor.tiff", sep = ""),  width = 16, heigh = 10)

		stem_expr = c("ZNF683", "CTSW", "LEF1", "NUCB2", "SCML4", "TXK", "CD55", "S100B", "MAL", "PLAC8")
		stem_ol = FeaturePlot(tc_cor, features = stem_expr, reduction = "umap",
						    min.cutoff = "q05", max.cutoff = "q95", ncol = 4, cols = c("orange1", "grey", "purple"))
		ggsave(plot = stem_ol, filename = paste(reclst_res_dir, "stem_expression_overlay_cor.tiff", sep = ""),  width = 16, heigh = 10)


		vln_feat = VlnPlot(tc_cor, features = expr, slot = "counts", log = TRUE, ncol = 3, pt.size = 0.5)
		ggsave(plot = vln_feat, filename = paste(reclst_res_dir, "violin_gene_feature_cor.tiff", sep = ""), width = 12, heigh = 16)
		cluster.markers <- FindAllMarkers(tc_cor, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
		write.csv(cluster.markers, paste(reclst_res_dir, "pca", ndim, "_cluster_marker_cor.csv", sep = ""))

		} else {

		tc <- FindVariableFeatures(tc, selection.method = "vst", nfeatures = 2000)	
		top10 <- head(VariableFeatures(tc), 10)
		# plot variable features with and without labels
		plot1 <- VariableFeaturePlot(tc)
		plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
		ggsave(plot = plot2, filename = paste(reclst_res_dir, "vari_feature.tiff", sep = ""),  width = 7, heigh = 6)

		all.genes <- rownames(tc)
		tc <- ScaleData(tc, features = all.genes)
		tc <- RunPCA(tc, features = VariableFeatures(object = tc))
		print(tc[["pca"]], dims = 1:5, nfeatures = 5)

		pca_vizdim = VizDimLoadings(tc, dim = 1:6, reduction = "pca")
		ggsave(plot = pca_vizdim, filename = paste(reclst_res_dir, "pca_vizdim_plot.tiff", sep = ""), width = 10, height = 15)
		
		pca_dim = DimPlot(tc, reduction = "pca")
		ggsave(plot = pca_dim, filename = paste(reclst_res_dir, "pca_dim_plot.tiff", sep = ""), width = 9, height = 6)

		tiff(paste(reclst_res_dir, "pca_heatmap.tiff", sep = ""), res = 180, width = 12, height = 16, units = "in")
		DimHeatmap(tc, dims = 1:15, cells = 500, balanced = TRUE)
		gar = dev.off()

		tc = JackStraw(tc, num.replicate = 100)
		tc = ScoreJackStraw(tc, dims = 1:20)
		js_plot = JackStrawPlot(tc, dims = 1:20)
		ggsave(plot = js_plot, filename = paste(reclst_res_dir, "JackStraw_plot.tiff", sep = ""), width = 9, heigh = 6)

		elb_plot = ElbowPlot(tc, ndims = 20)
		ggsave(plot = elb_plot, filename = paste(reclst_res_dir, "elbow_plot.tiff", sep = ""), width = 9, heigh = 6)

		saveRDS(tc, paste(reclst_res_dir, "seurat_obj_tc_bf_pc.RDS", sep = ""))

		tc = FindNeighbors(tc, dims = 1:ndim)
		tc = FindClusters(tc, resolution = 0.5)

		tc = RunUMAP(tc, dims = 1:ndim)
		cls_umap = DimPlot(tc, reduction = "umap")
		ggsave(plot = cls_umap, filename = paste(reclst_res_dir, ndim, "_cluster_umap_plot.tiff", sep = ""), width = 9, heigh = 6)
		btc = DimPlot(tc, reduction = "umap", group.by = "orig.ident")
		ggsave(plot = btc, filename = paste(reclst_res_dir, ndim, "batch_check_aft_cor.tiff", sep = ""), width = 9, heigh = 6)


		expr = c("PTPRC", "CD3E", "CD4", "CD8A", "TCF7", "IL7R", "ITGAE", "PDCD1", "ENTPD1", "DPP4")
		expr_ol = FeaturePlot(tc, features = expr, reduction = "umap",
						    min.cutoff = "q05", max.cutoff = "q95", ncol = 4, cols = c("orange1", "grey", "purple"))
		ggsave(plot = expr_ol, filename = paste(reclst_res_dir, "expression_overlay.tiff", sep = ""),  width = 16, heigh = 10)

		stem_expr = c("ZNF683", "CTSW", "LEF1", "NUCB2", "SCML4", "TXK", "CD55", "S100B", "MAL", "PLAC8")
		stem_ol = FeaturePlot(tc, features = stem_expr, reduction = "umap",
						    min.cutoff = "q05", max.cutoff = "q95", ncol = 4, cols = c("orange1", "grey", "purple"))
		ggsave(plot = stem_ol, filename = paste(reclst_res_dir, "stem_expression_overlay.tiff", sep = ""),  width = 16, heigh = 10)


		vln_feat = VlnPlot(tc, features = expr, slot = "counts", log = TRUE, ncol = 3, pt.size = 0.5)
		ggsave(plot = vln_feat, filename = paste(reclst_res_dir, "violin_gene_feature.tiff", sep = ""), width = 12, heigh = 16)
		cluster.markers <- FindAllMarkers(tc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
		write.csv(cluster.markers, paste(reclst_res_dir, "pca", ndim, "_cluster_marker.csv", sep = ""))



		}
	




}










####Identification of highly variable features (feature selection)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000) # nfeatures = (1000,2000)
top10 <- head(VariableFeatures(sc_obj), 10)
plot1 <- VariableFeaturePlot(sc_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave(plot=plot2, file= "varfeat_plot.tiff",width = 9, height = 6)

####Scaling the data
all.genes <- rownames(sc_obj)
sc_obj <- ScaleData(sc_obj, features = all.genes)

####Perform linear dimensional reduction

sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj))
print(sc_obj[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(sc_obj, dims = 1:2, reduction = "pca")
DimPlot(sc_obj, reduction = "pca")
DimHeatmap(sc_obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sc_obj, dims = 1:15, cells = 500, balanced = TRUE)


####Determine the ‘dimensionality’ of the dataset
sc_obj <- JackStraw(sc_obj, num.replicate = 100)
sc_obj <- ScoreJackStraw(sc_obj, dims = 1:20)
JackStrawPlot(sc_obj, dims = 1:20)

ElbowPlot(sc_obj)
sc_obj <- FindNeighbors(sc_obj, dims = 1:10)

sc_obj <- FindClusters(sc_obj, resolution = 0.5)

head(Idents(sc_obj), 5)

####Run non-linear dimensional reduction (UMAP/tSNE)
sc_obj <- RunUMAP(sc_obj, dims = 1:10)
DimPlot(sc_obj, reduction = "umap")
saveRDS(sc_obj, file = "../output/pbmc_tutorial.rds")

cluster1.markers <- FindMarkers(sc_obj, ident.1 = 1, min.pct = 0.25)
rnasc_obj.markers <- FindAllMarkers(sc_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

rnasc_obj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.csv(rnasc_obj.markers, ".csv")

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet", "unknown1","unknown2","unknown3","unknown4", "unknown5","unknown6")
names(new.cluster.ids) <- levels(sc_obj)
sc_obj <- RenameIdents(sc_obj, new.cluster.ids)
DimPlot(sc_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}
reference.list <- integrate.list[c("tumor", "normal", "LN1", "LN2")]
integrate.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)



LN2<-Read10X(data.dir = "Z:/Group/scmodb/vdj_citeseq_expr_112019/30928_20_LN2_raw_feature_bc_matrix/")
 oldRowName = rownames(LN2$`Antibody Capture`)
 newRowName = str_split_fixed(oldRowName, "_", n = 4)[,3]
 rownames(LN2$`Antibody Capture`) = newRowName


LN2obj = CreateSeuratObject(counts = LN2$`Gene Expression`, project = saName[4])

 
LN2obj[["Protein"]] = CreateAssayObject(counts = LN2$`Antibody Capture`)

LN2obj[["percent.mt"]] = PercentageFeatureSet(LN2obj, pattern = "^MT-")
LN2obj = subset(LN2obj, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & percent.mt < 5)
 
LN2obj = NormalizeData(LN2obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

LN2obj = NormalizeData(LN2obj, assay = "Protein", normalization.method = "CLR")
 
LN2obj = FindVariableFeatures(LN2obj, selection.method = "vst", nfeatures = 2000)

allGenes = rownames(LN1obj[["RNA"]])
LN1obj = ScaleData(LN1obj, features = allGenes, assay = "RNA")

LN1obj = ScaleData(LN1obj, assay = "Protein")
LN1obj = RunPCA(LN1obj, features = VariableFeatures(object =LN1obj)) 
ElbowPlot(LN1obj)
LN1obj = FindNeighbors(LN1obj, dims = 1:15)
LN1obj = FindClusters(LN1obj, resolution = 0.5)


LN1obj = RunUMAP( LN1obj, dims = 1:15)

DimPlot(LN1obj, reduction = "umap")
protFeas = paste0("protein_", rownames(LN1obj[["Protein"]]))
FeaturePlot(LN1obj, features = protFeas, reduction = "umap",
            min.cutoff = "q05", max.cutoff = "q95", ncol = 4, cols = c("orange1", "grey", "purple"))

assign(saName[3], LN1obj)

integrateList = list("tumor" = bc377_xx_TUMOR,
			     "normal" = bc389_xx_TUMOR, 
			     "LN1" = bc392_xx_TUMOR,
			     "LN2" = bc393_xx_TUMOR
			     )