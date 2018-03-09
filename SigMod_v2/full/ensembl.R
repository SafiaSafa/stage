#Library
#source("https://bioconductor.org/biocLite.R")
#biocLite("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
#source("https://bioconductor.org/biocLite.R")
#biocLite("genomation")
library(genomation)
library(rtracklayer)

# get snp location
#plink --bfile ~/Documents/data/BP_final/BP.B37-final --recode --out BP.B37-final

bed <- as.data.frame(read.table("~/stage/SigMod_v2/full/BP.B37-final.map",header = FALSE, sep="\t",stringsAsFactors=FALSE))
colnames(bed)=c( "chr","RefSNP_id", "start","end")
bed$start=bed$end
gr <- makeGRangesFromDataFrame(bed,
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqnames.field="chr",
                         start.field="start",
                         end.field="end")
#check position SNPs on genome assemble GRCh37
#snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
#my_snps <- snpsById(snps,"rs3934834")
#my_snps

## get all genes annotations
edb <- EnsDb.Hsapiens.v75
gns<-genes(edb)

#fix notation chromosome to: chr
seqlevelsStyle(gr) <- "UCSC" 
seqlevelsStyle(gns) <- "UCSC" 


# Remove pbl mitochondrial DNA
gns <- keepSeqlevels(gns, seqlevels(gns)[1:24])

# Unify genome assembly hg19 <-> grch37
genome(gr) <- "GRCh37"

#increase the size of the gense by 5000(pb)
gns2 <- gns
start(gns2) <- start(gns2) - 10000
end(gns2) <- end(gns2) + 10000

# get the overlapping data, 
hits <- findOverlaps(gr, gns2)

#get the result table of the overlapping data
grolaps <- gr[queryHits(hits),]
gnolaps <- gns[subjectHits(hits),]
mcols(gnolaps)$symbol <- mcols(gnolaps)$gene_name
out <- cbind(as(gnolaps$symbol, "DataFrame"),as(grolaps$RefSNP_id, "DataFrame"))
colnames(out)<-c('gene','SNP')
write.table(out,"~/stage/SigMod_v2/full/gene_snp.txt",row.names=FALSE, quote = F)
