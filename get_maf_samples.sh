#!/gsc/bin/bash


<< COMMENTS

The purpose of this script is to take maf files, and get a list of genes and their samples that are mutated in this particular maf. 

COMMENTS

while getopts hi:s:g:m: opt; do
  case $opt in
  i)
      INPUT_FILE=$OPTARG
      ;;
  s)
      SAMPLE_INDEX=$OPTARG
      ;;
  g)
      GENE_INDEX=$OPTARG
      ;;
  m)
      MUTATIONTYPE_INDEX=$OPTARG
      ;;
  h)
     echo "This program has many mandatory and optional arguments:
      
	Mandatory Arguments:
	-i : input file should be a MAF file
	-g : gene index (usually 0)
	-s : sample index 
	-m : mutation type index
	
	Optional Arguments:
	-h : This help message

	"
	exit 1 
   esac
done

shift $((OPTIND - 1))

if [[ -z ${INPUT_FILE+x} ]] || [[ -z ${GENE_INDEX+x} ]] || [[ -z ${SAMPLE_INDEX+x} ]] || [[ -z ${MUTATIONTYPE_INDEX+x} ]]
then
	echo "One of the inputs missing"
fi


cut -f $SAMPLE_INDEX $INPUT_FILE > $INPUT_FILE.samples
cut -c 1-12 $INPUT_FILE.samples > $INPUT_FILE.samplescut
cut -f $GENE_INDEX,$MUTATIONTYPE_INDEX $INPUT_FILE > $INPUT_FILE.other
paste -d'\t' $INPUT_FILE.samplescut $INPUT_FILE.other > $INPUT_FILE.mutatedsamples
rm $INPUT_FILE.samplescut $INPUT_FILE.other $INPUT_FILE.samples
	


