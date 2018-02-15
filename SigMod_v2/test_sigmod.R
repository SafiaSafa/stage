R_dir = 'R/' ## the path to the SigMod_v2\R folder
## load all R functions in this folder
file.sources = list.files(R_dir, pattern='*.R$',
                          full.names=TRUE, ignore.case=TRUE)
tmp=sapply(file.sources, source, .GlobalEnv)