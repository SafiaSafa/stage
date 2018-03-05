#require( igraph )
#if(packageVersion("igraph") < "1.0.0") {
#stop("Need to install igraph version 1.0.0 or above")
#}
library(igraph)
R_dir = '~/stage/SigMod_v2/R/' ## the path to the SigMod_v2\R folder
## load all R functions in this folder
file.sources = list.files(R_dir, pattern='*.R$',
full.names=TRUE, ignore.case=TRUE)

tmp=sapply(file.sources, source, .GlobalEnv)

############## Read gene p-values data ##############
## For example, your p-values data has a tab delimited text format
gene_ps_file = '~/stage/SigMod_v2/computed_gene_p_values.tab'
gene_ps = read.table(gene_ps_file, stringsAsFactors = FALSE, header =
TRUE)
## construct the string-based scored-network:
#source("https://bioconductor.org/biocLite.R")
#biocLite("STRINGdb")
library(STRINGdb)
scored_net = construct_string_net( gene_ps = gene_ps )

## identify module from the scored_net
## scored_net is the return of construct_scored_net function or construct_string_net function
## results will be saved under the working folder
res_info = SigMod_bisection( net=scored_net )


