
library(maftools)

#data_dir = "//isi-dcnl/user_data/plee/Group/public_data/tcga/"
res_dir = "Y:/Joyce/TCGA_pancancer/mutation/"
cancertype = "breast"
res_folder = paste(res_dir, cancertype, "_varscan/", sep = "")
TCGAtype = "BRCA"

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

lollipopPlot(maf = mafObj, gene = 'TP53', AACol = 'HGVSp_Short', showMutationRate = TRUE)

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

tiff(paste(res_folder, "msomatic_interaction_corr_plot.tiff", sep = ""), res = 180, width = 6, height = 6, units = "in")
somaticInteractions(maf = mafObj, top = 25, pvalue = c(0.05, 0.1))
dev.off()

soma_inter = somaticInteractions(maf = mafObj, top = 25, pvalue = c(0.05, 0.1))
write.csv(soma_inter, paste(res_folder, "somatic_interaction_summary.csv", sep = 


sample_total = getSampleSummary(mafObj)
all_mut_gene=mafObj@data
tex_sig_mut=subset(tex_mut_all, cancer_type == "crc")
sample_reid=sample_total
sample_reid$Tumor_Sample_Barcode = substr(sample_reid$Tumor_Sample_Barcode, 1,12)
sample_reid$Tumor_Sample_Barcode=gsub("-",".",sample_reid$Tumor_Sample_Barcode)

tex_sig_mut=subset(tex_sig_mut, X %in% sample_reid$Tumor_Sample_Barcode)

sub_scres=subset(tex_sig_mut, group !='Medium')
stop()
for (ig in 1:length(mut_gene)) {

    sub_mut=subset(all_mut_gene, Hugo_Symbol== mut_gene[ig])
    sample_total$mut_gene ="FALSE"
    sample_total$mut_gene[sample_total$Tumor_Sample_Barcode %in% sub_mut$Tumor_Sample_Barcode] = "TRUE"

    T=sum(sample_total$mut_gene == "TRUE")
    F=sum(sample_total$mut_gene == "FALSE")

    cat( T, cancertype, "patients had ", mut_gene[ig], " mutation", "\n")

    cat( F, cancertype, "patients didn't have ", mut_gene[ig], " mutation", "\n")


     
    sample_reid=sample_total
    sample_reid$Tumor_Sample_Barcode = substr(sample_reid$Tumor_Sample_Barcode, 1,12)


    sample_reid$Tumor_Sample_Barcode=gsub("-",".",sample_reid$Tumor_Sample_Barcode)

    tex_sig_mut$mut_gene=sample_reid$mut_gene[match(tex_sig_mut$X, sample_reid$Tumor_Sample_Barcode)]




    tex_sig_mut$mut_gene=factor(tex_sig_mut$mut_gene, levels = c("TRUE", "FALSE"))
    g<-ggboxplot(tex_sig_mut, x = "mut_gene", y = "sig_score",
                 color = "mut_gene", 
                 add = "jitter") 

    c =g+theme_classic()+theme(legend.position="none", axis.title.x=element_text(size=14, color = "#000000"),
                               axis.ticks.x = element_blank(),axis.text=element_text(size=10, color = "#000000"),
                               axis.title=element_text(size=14))+ ylab("Tex sig.score") +xlab(paste(mut_gene[ig], "Mutation", sep = " "))
    d<-c+stat_compare_means(method = "t.test", label.x = 1.5,  size = 5, label = "p.format")
    ggsave(d, file = paste(cancertype, "_Tex_sig.score_", mut_gene[ig],"_mutation_tcga_ttes_ER.tiff",sep = "") , width = 4, height=6,dpi= 300, units = "in", device = "tiff")


    col=dim(sample_total)[2]
    colnames(sample_total)[col] = mut_gene[ig]

    col=dim(tex_sig_mut)[2]
    colnames(tex_sig_mut)[col] = mut_gene[ig]

    }


subtype_clin$DFS<-ifelse(subtype_clin$TLR>subtype_clin$TDR,subtype_clin$TDR,subtype_clin$TLR) # dfs includes local relapse and distant replase

if (tex_sig_mut$X %in% test2$x) {tex_sig_mut$pik3ca_mut_fre=test2$fre}
 else if (tex_sig_mut$X %in% sub_mut_id $Tumor_Sample_Barcode) {tex_sig_mut$pik3ca_mut_fre =1}
 else {tex_sig_mut$pik3ca_mut_fre=0}


mut_gene = ""
sub_mut=subset(all_mut_gene, Hugo_Symbol== mut_gene)

sub_mut$Tumor_Sample_Barcode=substr(sub_mut$Tumor_Sample_Barcode, 1, 12)
sub_mut$Tumor_Sample_Barcode=gsub("-", ".", sub_mut$Tumor_Sample_Barcode)
dup_mut=sub_mut[duplicated(sub_mut$Tumor_Sample_Barcode)]
dup_id=subset(sub_mut, Tumor_Sample_Barcode %in% dup_mut$Tumor_Sample_Barcode)
dup_freq=count(dup_id$Tumor_Sample_Barcode)

sub_tex_mut=subset(tex_sig_mut, X %in% sub_mut$Tumor_Sample_Barcode)

sub_tex_mut$TP53_freq=1
sub_tex_mut$TP53_freq[sub_tex_mut$X %in% dup_freq$x]=dup_freq$freq




mut_sum=as.data.frame(matrix(ncol = 5, nrow = length(mut_gene)))
colnames(mut_sum) = c("num_TRUE", "num_FALSE", "mean_TRUE", "mean_FALSE", "p.value")
rownames(mut_sum)=mut_gene


for (ig in 1:length(mut_gene)) {

    sub_mut=subset(all_mut_gene, Hugo_Symbol== mut_gene[ig])
    sample_total$mut_gene ="FALSE"
    sample_total$mut_gene[sample_total$Tumor_Sample_Barcode %in% sub_mut$Tumor_Sample_Barcode] = "TRUE"

    ##cat( T, cancertype, "patients had ", mut_gene[ig], " mutation", "\n")

    ##cat( F, cancertype, "patients didn't have ", mut_gene[ig], " mutation", "\n")


     
    sample_reid=sample_total
    sample_reid$Tumor_Sample_Barcode = substr(sample_reid$Tumor_Sample_Barcode, 1,12)


    sample_reid$Tumor_Sample_Barcode=gsub("-",".",sample_reid$Tumor_Sample_Barcode)

    tex_sig_mut$mut_gene=sample_reid$mut_gene[match(tex_sig_mut$X, sample_reid$Tumor_Sample_Barcode)]




    tex_sig_mut$mut_gene=factor(tex_sig_mut$mut_gene, levels = c("TRUE", "FALSE"))

    T=sum(tex_sig_mut$mut_gene == "TRUE")
    F=sum(tex_sig_mut$mut_gene == "FALSE")
   

    mut_sum[ig,1]=T
    mut_sum[ig,2]=F

    if (T > 1) {
    	stat_test = t.test(sig_score~mut_gene, data = tex_sig_mut)

	    mean_TRUE=stat_test$estimate[1]
	    mean_FALSE=stat_test$estimate[2]
	    p.value=stat_test$p.value

	    mut_sum[ig,3]=mean_TRUE
	    mut_sum[ig,4]=mean_FALSE
	    mut_sum[ig,5]=p.value
	}


	    col=dim(sample_total)[2]
    colnames(sample_total)[col] = mut_gene[ig]

    col=dim(tex_sig_mut)[2]
    colnames(tex_sig_mut)[col] = mut_gene[ig]

    }

all_mut=unique(sig_sum$X)

sum_mut=as.data.frame(matrix(nrow=length(all_mut), ncol=1))
rownames(sum_mut)=all_mut
colnames(sum_mut)=c('number','cancer_type')
for (ig in 1:length(all_mut)) {
	igene=all_mut[ig]
	sub_mut=subset(sig_sum, X == igene)
	cancertype = paste(sub_mut$cancer_type, collapse = ",")
	sum_mut[ig,1]=dim(sub_mut)[1]
	sum_mut[ig,2]=cancertype
}
