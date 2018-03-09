
gene_liste = '~/stage/GRIV/selected_genes.tab'
gene_liste = read.table(gene_liste, stringsAsFactors = FALSE, header =
                            TRUE)

#########Psygenet data Bipolar
data=read.table("~/stage/GRIV/HIV-1_Interactions.csv",sep=',', header = TRUE)



common_data<-intersect(as.character(data$Human_GeneSymbol),gene_liste$gene)

length(common_data)
#[1] 59