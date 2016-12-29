#!/usr/bin/bash

<<COMMENTS
The purpose of this code is to intersect the copy number file and the bed file to come up with a single file that has each gene mapped to a bunch of regions,  sample ids,  and copy number polidy numbers. 

COMMENTS

FILES=$1/*.seg.txt

for file in $FILES
	do
	sed '1d' $file > $file.noheader
	cut -f 2-6 $file.noheader > $file.first
	cut -f 1 $file.noheader > $file.second
	paste -d'\t' $file.first $file.second > $file.rearranged
	rm $file.noheader $file.first $file.second
	echo "bedtools intersect  -a $2 -b $file.rearranged -wa -wb > $file.intersected" | qsub 
	done
