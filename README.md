# Classy: A Bioinformatics Approach to Discovering Novel Tumor Suppressors from Oncogenes 

The purpose of CLASSY is to extract and synthesize important biological imformation for data provided by THE CANCER GENOME ATLAS PROJECT. It extracts biologically relevant features from mutation data provided through Mutation Anotation Files (MAF), copy number data provided through Copy Number Annotation files, and gene expression data provided through Expression files to classify genes as tumor suppressors or oncogenes. The information below includes all of the code used in this process, and provides a step by step explanation as to how to run the code along with the input and output for each. 

The test_data directory has three test files:

BRCA.TP53.maf -- MAF data for TP53 from BRCA Cancer Data

BRCA.TP53.cna -- CNA data for TP53 from BRCA Cancer Data

BRCA.TP53.exp -- EXP data for TP53 from BRCA Cancer Data

The directory also contains supplementary files: 

bed_file -- a sample bed file

MuSic_Matrix -- a matrix provided by MuSic that can tell us which cancer types to aggregate data from. This is completely optional, and was the method we used to filter our data for irrelevant noise. 

# MAF DATA ANALYSIS

get_maf_samples.sh
The purpose of this script is to take maf files, and get a list of genes and their samples that are mutated in this particular maf. Takes in an input maf file, a sample index (usually field 12 in a MAF), a gene index (usually the first field in a MAF file), and a mutation type index (usually field 9 in a MAF file). Note: Please derive the indicies from your MAF file for use in this script.

   USAGE:
   
       bash get_maf_samples.sh -i BRCA.TP53.maf -s 12 -g 1 -m 9
       
   INPUT FILE: BRCA.TP53.maf
   
   OUTPUT FILE: BRCA.TP53.maf.mutatedsamples
   
parse_maf_file.pl
The purpose of this script is to extract mutation data from MAF files. Particularly we are interested in the number of truncating (Loss of Function) mutations, number of recurrent missense mutations at the codon specific level, and the number of recurrent domain mutations. NOTE: Please derive the indices from your MAF file for use in this script.

   USAGE:
   
       perl parse_maf_file.pl --maf_file=BRCA.TP53.maf --gene_index=0 --mut_index=8 --sample_index=15 --codon_index=41 --domain_index=45

   INPUT FILE: BRCA.TP53.maf
   
   OUTPUT FILE: BRCA.TP53.maf.output


# GENE EXPRESSION ANALYSIS

parse_expression_file.R

USAGE:

      Rscript parse_expression_file.R BRCA.TP53.exp BRCA.TP53.mutatedsamples BRCA



**********FOR COPY NUMBER FILES***********

FIRST COMMANDLINE EXECUTION
bash createCNAfile.sh /directory/of/cna/files bed_file

SECOND COMMANDLINE EXECUTION
bash cnascript.sh

THIRD COMMANDLINE EXEUCTION
bash filtercna.sh


********AFTER INDIVIDUAL FILE ANALYSIS********
FIRST COMMANDLINE EXECUTION
bash mergeMafbyCancerType.sh listOfCancers directory_to_mafoutputs director_to_cnaoutputs directory_to_geneexpressionoutputs

SECOND COMMANDLINE EXECUTION
Rscript pancanMerge.R


********CLASSIFICATION*********
Now we have all of our training and test data we can classify

Rscript classy.R trainingdata testdata output_file_name
