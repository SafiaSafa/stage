#! /bin/bash


export PATH=$PATH:/media/vcabeliNFS2/safia/bin
liste="0 1 2 5"
for ele in $liste
do
	qsub -varg=${ele} -lnodes=1:ppn=8 ./rank_combined.sh -e ${i}/qsub_erreur${i}.e -o ${i}/qsub_erreur${i}.o
done
