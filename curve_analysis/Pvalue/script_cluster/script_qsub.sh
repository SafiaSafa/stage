#! /bin/bash


export PATH=$PATH:/media/vcabeliNFS2/safia/bin
liste="0 1 2 5"
for ele in $liste
do
	qsub -varg=${ele} -lnodes=1:ppn=8 ./test_python.sh -e ${i}/qsub_erreur${i}.e -o ${i}/qsub_erreur${i}.o
done
