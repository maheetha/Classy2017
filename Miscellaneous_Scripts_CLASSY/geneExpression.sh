#!/usr/bin/bash

FILES=*.finaledited
counter=0
arr=($1/*.new)

for file in $FILES
	do
	samples=${arr[$counter]}
	bsub -oo $file.out -e $file.err Rscript geneExpression.R $file $samples $file.rsem $file
	counter=$((counter+1))
	done
