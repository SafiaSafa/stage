

#awk '{print$2" "$1" "$3" "$9}' low_ld.assoc.logistic > ~/stage/SigMod_v2/snp_pval.txt
#var='SNP chr pos p'
#sed -i "1s/.*/$var/" ~/stage/SigMod_v2/snp_pval.txt

source("/home/vcabeli/stage/SigMod_v2/full/fastCGP/main.R")
results = fastCGP( "~/stage/SigMod_v2/full/gene_snp.txt", "~/stage/SigMod_v2/full/snp_pval.txt")
write.table(results,"~/stage/SigMod_v2/full/computed_gene_p_values.tab",sep='\t',row.names=FALSE, quote = F)

