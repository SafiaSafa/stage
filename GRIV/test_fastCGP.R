

#awk '{print$2" "$1" "$3" "$9}' low_ld.assoc.logistic > ~/stage/GRIV/snp_pval.txt
#var='SNP chr pos p'
#sed -i "1s/.*/$var/" ~/stage/GRIV/snp_pval.txt

source("/home/vcabeli/stage/GRIV/fastCGP/main.R")
results = fastCGP( "~/stage/GRIV/gene_snp.txt", "~/stage/GRIV/snp_pval.txt")
write.table(results,"~/stage/GRIV/computed_gene_p_values.tab",sep='\t',row.names=FALSE, quote = F)

