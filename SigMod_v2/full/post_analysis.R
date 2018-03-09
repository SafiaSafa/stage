
gene_liste = '~/stage/SigMod_v2/selected_genes.tab'
gene_liste = read.table(gene_liste, stringsAsFactors = FALSE, header =
                            TRUE)
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
#########Ensembl name convert to Ensembl IDs For DAVID analysis
table2 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = "external_gene_name",
                values =  gene_liste$gene, 
                mart = grch37)
write.table(table2$ensembl_gene_id,"~/stage/SigMod_v2/gene_liste.txt",sep='\t',row.names=FALSE,col.names=FALSE, quote = F)

#########Ensembl IDs convert to HUGO symbole For Psygenet analyses
table3 <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters = "ensembl_gene_id", values = table2$ensembl_gene_id, mart = grch37)
write.table(table3$hgnc_symbol,"~/stage/SigMod_v2/gene_hugo_symbole.txt",sep='\t',row.names=FALSE,col.names=FALSE, quote = F)

#########Psygenet data Bipolar
data=read.table("~/stage/SigMod_v2/tab3b.tsv",sep='\t', skip=1)
first <- readLines("~/stage/SigMod_v2/tab3b.tsv", n=1)
first <-unlist(strsplit(first, "\t"))
colnames(data) <- first
common_data<-intersect(as.character(data$`Gene Symbol`),table3$hgnc_symbol)

vec<-rep(0,length(common_data))
for (i in 1:length(common_data)){
vec[i]=which(data$`Gene Symbol`==common_data[i],arr.ind = T)
}
write.table(data[vec,c(1,3)],"~/stage/SigMod_v2/common_genes_psygenet.txt",sep='\t',row.names=FALSE,col.names=TRUE, quote = F)

