


#data_dir = "//isi-dcnl/user_data/plee/Group/public_data/tcga/"
res_dir = "Y:/Joyce/TCGA/mutation/"
cancertype = "stomach"
res_folder = paste(res_dir, cancertype, "_varscan/", sep = "")
TCGAtype = "STAD"

dir.create(file.path(res_folder), showWarnings = FALSE)
#


maf.file = paste(res_folder, "TCGA.", TCGAtype, ".varscan.maf.gz", sep = "")

mafObj = read.maf(maf = maf.file)

tiff(paste(res_folder, "mutation_summary.tiff", sep = ""), res = 180, width = 9, height = 6, units = "in")
plotmafSummary(maf = mafObj, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

#oncoplot for top ten mutated genes.
tiff(paste(res_folder, "mutation_top10_genes.tiff", sep = ""), res = 180, width = 9, height = 6, units = "in")
oncoplot(maf = mafObj, top = 10)
dev.off()


mafObj.titv = titv(maf = mafObj, plot = FALSE, useSyn = TRUE)
#plot titv summary
tiff(paste(res_folder, "mutation_snp_titv.tiff", sep = ""), res = 180, width = 6, height = 6, units = "in")
plotTiTv(res = mafObj.titv)
dev.off()



if(FALSE) {
	#Shows sample summry.
	sample = getSampleSummary(mafObj)
	sample$Tumor_Sample_Barcode = substr(sample$Tumor_Sample_Barcode, 1,12)
	write.csv(sample, paste(res_folder, "mafObj_sampleSummary.csv", sep = ""))

	#Shows sample summry.
	#Shows gene summary.
	#getGeneSummary(mafObj)
	#shows clinical data associated with samples
	#getFields(mafObj)
	#Writes maf summary to an output file with basename mafObj.
	write.mafSummary(maf = mafObj, basename = res_folder)
}


