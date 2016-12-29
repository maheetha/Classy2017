#!/gsc/bin/bash


<< COMMENTS

The purpose of this script is to take a directory of MAF (Mutation Annotation Files) and create three files from that: a file with truncation mutation information ( file suffix '.trunc'), a file that contains all of the positions at which there is a missense mutation (file suffix '.recur'), and file that contains all of the missense domains (file suffix '.missensedomains')

$1 = should be the directory that contain original mafs
$2 = should be the directory where you want them to go


COMMENTS

FILES=$1/*.maf


for file in $FILES
	do
	echo $file
	bsub -o $file.out -e $file.err ./maf_script.pl --maf_file=$file --gene_index=0 --mut_index=8 --sample_index=15 --codon_index=41 --domain_index=49		
	cut -f16 $file > $file.samples
        cut -c1-12 $file.samples > $file.samplescut
        cut -f1,9 $file > $file.other
        paste -d'\t' $file.samplescut $file.other > $file.new
        rm $file.samplescut $file.other
	done


