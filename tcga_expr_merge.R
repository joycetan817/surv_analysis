
tcga_skcm_dir = "Z:/plee/Group/public_data/tcga/mrna_expr/crc/"
data_type = "FPKM_UQ"#"cts" "FPKM_UQ"
data_dir = paste(tcga_skcm_dir, data_type, "/", sep ="")
setwd(data_dir)
files = list.files(path = data_dir)
expr <- read.table(files[1], row.names=1)

file <-as.data.frame(list.files(path = data_dir))
colnames(file)[1]="file_name"
sample_annot<- read.table("Z:/plee/Group/public_data/tcga/mrna_expr/crc/info/gdc_sample_sheet.2019-10-03.tsv", sep="\t", header = TRUE) 
sample_annot$File.Name<-sub(".gz", "", sample_annot$File.Name)
file$pid<-sample_annot$Sample.ID[match(file$file_name, sample_annot$File.Name)]
file$pid<-as.character(file$pid)
colnames(expr)[1]=file$pid[file$file_name==files[1]]

st = Sys.time()
for ( ii in 2:length(files)) {
	expr_temp <- read.table(files[ii], row.names=1)
	colnames(expr_temp)[1]=file$pid[file$file_name==files[ii]]
	expr<-cbind(expr, expr_temp)
}
print(Sys.time()-st)
