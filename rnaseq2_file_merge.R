library("rjson")
skcm_annot <- fromJSON(file = "metadata.cart.2019-10-16.json")


annot<-as.data.frame(matrix(ncol=2, nrow=length(skcm_annot)))	
colnames(annot)[1]="file_name"
colnames(annot)[2]="pid"

for (i in 1:length(skcm_annot)) { 
	annot[i, 1] = skcm_annot[[i]][["file_name"]]	
}

for (i in 1:length(skcm_annot)) { 
	annot[i, 2] = skcm_annot[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
}

annot$sample_id<-substr(annot$pid, 1, 12)

tcga_dir = "Z:/plee/Group/public_data/tcga/mrna_expr/skcm/"
data_type = "NORM_RSEM"#"FPKM_UQ"#"cts" "FPKM_UQ"
data_dir = paste(tcga_dir, data_type, "/", sep ="")
setwd(data_dir)
files = list.files(path = data_dir)
expr <- read.table(files[1], header = TRUE, row.names=1)
#expr[,2:3]=NULL When merge raw counts
colnames(expr)[1]=annot$sample_id[annot$file_name==files[1]]


st = Sys.time()
for (ii in 2:length(files)) {
	expr_temp <- read.table(files[ii], header = TRUE, row.names=1)
	#expr_temp[,2:3]=NULL
	colnames(expr_temp)[1]=annot$sample_id[annot$file_name==files[ii]]
	expr<-cbind(expr, expr_temp)
}
print(Sys.time()-st)



