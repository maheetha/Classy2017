#!/usr/bin/bash

<<COMMENTS
The purpose of this code is to intersect the copy number file and the bed file to come up with a single file 
that has each gene mapped regions to sample ids, genes, and copy number polidy numbers. 

INPUT: Copy Number File from TCGA Portal usually ending with *.seg.txt 

SAMPLE FORMAT: 

Sample	Chromosome	Start	End	Num_Probes	Segment_Mean
TCGA-05-XXXX-01A-01D-XXXX-XX	1	3218610	8925112	3421	-0.545


OUTPUT: An intersected bed file with 9 regions

SAMPLE FORMAT: (no header in actual ending file)
CHR	START	END	GENE	CHR	START	END	Num_Probes	Segment_Mean	SAMPLE
1	2985731	3355185	PRDM16	1	3218610	3225940	15	-0.182	TCGA-55-8090-01A-11D-2237-01

COMMENTS

while getopts hi:o:b: opt; do
  case $opt in
  i)
      INPUT_FILE=$OPTARG
      ;;
  o)
      OUTPUT_FILE=$OPTARG
      ;;
  b)
      BED_FILE=$OPTARG
      ;;
  h)
     echo "This program has many mandatory and optional arguments:
      
	Mandatory Arguments:
	-i : input file
	-b : bed file, with the first three columns chr, start, and end.
	
	
	Optional Arguments:
	-h : This help message
	-o : if you have a selected output file name
	
	"
	exit 1 
   esac
done

shift $((OPTIND - 1))

if [[ -z ${INPUT_FILE+x} ]] || [[ -z ${BED_FILE+x} ]]
then
	echo "Input or Bed File Missing"
fi

if [ -z ${OUTPUT_FILE+x} ]
then
	OUTPUT_FILE=$INPUT_FILE.intersected
fi

######################################################################
### Rearrangement of files and intersection with bed file ######
### Note: Due to size of file, bedtools command might need #####
### submission to a cluster or higher memory node ##############
######################################################################

sed '1d' $INPUT_FILE > $INPUT_FILE.noheader
cut -f 2-6 $INPUT_FILE.noheader > $INPUT_FILE.first
cut -f 1 $INPUT_FILE.noheader > $INPUT_FILE.second
paste -d'\t' $INPUT_FILE.first $INPUT_FILE.second > $INPUT_FILE.rearranged
rm $INPUT_FILE.noheader $INPUT_FILE.first $INPUT_FILE.second
echo "bedtools intersect  -a $BED_FILE -b $INPUT_FILE.rearranged -wa -wb > $OUTPUT_FILE" | qsub  
done
