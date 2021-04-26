install.packages('DeconRNASeq')
library(DeconRNASeq)

data(multi_tissue)
datasets <- x.data[,2:11]
signatures <- x.signature.filtered.optimal[,2:6]
proportions <- fraction

result <- DeconRNASeq(datasets,signatures,fig=F)
final <- result$out.all

# own analysis
my.datasets <- read.table('/Users/ligk2e/Desktop/tmp/kidney/DeconRNASeq/bulk.txt',header=T,sep='\t',
                          row.names=1)
my.signatures <- read.table('/Users/ligk2e/Desktop/tmp/kidney/DeconRNASeq/sc_signature.txt',header=T,
                            sep='\t',row.names=1)
my.result <- DeconRNASeq(my.datasets,my.signatures,fig=F,use.scale = T)
my.final <- my.result$out.all
rownames(my.final) <- colnames(my.datasets)
write.table(my.final,'/Users/ligk2e/Desktop/tmp/kidney/DeconRNASeq/final.txt',sep='\t',col.names = NA)