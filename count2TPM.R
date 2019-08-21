setwd("C:\\Users\\liter_frye\\Desktop")
library(data.table)
#---------本程序可完成人类基因组count转换TPM-----------
#homo为自动文件
homo=read.table("gene_length.txt",header=T)
#f为文件第一列为gene_id,第二列为gene_name(第二列是什么不重要，但是需要一个空行，也可自行添加)
f=fread("GSM2758471_PJ016.filtered.matrix.txt",header=F)
f=data.frame(f)
f[,2]=0
colnames(f)[1:2]=c('gene_id','length')
count=merge(f,homo,by='gene_id')
count[,2]=count[,ncol(count)]
count=count[,-ncol(count)]
tpm=data.frame(gene_id=count[,1])
for(i in 3:ncol(count)){
  col1=data.frame(count[,i]/count$length)
  col1[which(is.na(col1)),1]=0
  col1[sapply(col1, is.infinite),]=0
  tpm[,i-1]=col1*1000000/sum(col1)
}
write.table(tpm,"tpm_patient16.txt",col.names =F,row.names = F)
