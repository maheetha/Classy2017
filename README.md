# Classy: A Bioinformatics Approach to Discovering Novel Tumor Suppressors from Oncogenes 

The purpose of CLASSY is to extract and synthesize important biological imformation for data provided by THE CANCER GENOME ATLAS PROJECT. It extracts biologically relevant features from mutation data provided through Mutation Anotation Files (MAF), copy number data provided through Copy Number Annotation files, and gene expression data provided through Expression files to classify genes as tumor suppressors or oncogenes. The information below includes all of the code used in this process, and provides a step by step explanation as to how to run the code along with the input and output for each. 

The test_data directory has three test files:

BRCA.TP53.maf -- MAF data for TP53 from BRCA Cancer Data

BRCA.TP53.cna -- CNA data for TP53 from BRCA Cancer Data

BRCA.TP53.exp -- EXP data for TP53 from BRCA Cancer Data

TP53.bedfile -- a sample bed file

Note: the BRCA.TP53.cna contains all samples from the BRCA cancer type, however, the bedfile contains only TP53. When you run intersect_bed_cna.sh, please make sure that the intersected file only has TP53. 

The directory also contains optional supplementary files: 

MuSic_Matrix_Test -- a matrix provided by MuSic that can tell us which cancer types to aggregate data from. This is completely optional, and was the method we used to filter our data for irrelevant noise. 

trainingdata -- this is the training data that one could use for classification 

traininggenes -- this is a compilation of genes we selected for our training set 

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

The purpose of this script is to extract features from expression files. Particularly we are interested calculating the average expression value for this gene in mutated samples and non-mutated samples. Within the script, we filter for any mutation that we haven't used in our calculation of the maf features. In the end, we hope to use a ratio between expression value in mutated versus nonmutated samples as one feature.

USAGE:

      Rscript parse_expression_file.R BRCA.TP53.exp BRCA.TP53.mutatedsamples BRCA.TP53

INPUT FILE: BRCA.TP53.exp, BRCA.TP53.mutatedsamples
OUTPUT FILE: BRCA.TP53.rsem


# COPY NUMBER ANALYSIS 

intersect_bed_cna.sh

The purpose of this script is to intersect a bed file and the relevant copy number file in order to match every copy number region with a gene.

USAGE:

      bash intersect_bed_cna.sh -i BRCA.TP53.cna -b bed_file
      
INPUT FILE: BRCA.TP53.cna
OUTPUT FILE: BRCA.TP53.cna.intersected

parse_cna_file.pl

The purpose of this script is to extract features from copy number files. Particularly we are interested looking at a copy number distribution. We create 8 different "buckets" of ploidy: 0-0.5, 0.5-1, 1-1.5, and so on and so forth, until we reach 4.0. Anything above 4.0 is neglected for now, but we hope to (in later versions), to be able to understand how to better handle (or normalize) that data. As features, we hope to calculate the proportion of samples that fall in each of the different buckets. 

USAGE:

     perl parse_cna_file.pl --cna_file=BRCA.TP53.cna.intersected --maf_sample=BRCA.TP53.mutatedsamples
     
INPUT FILE: BRCA.TP53.cna.intersected, BRCA.TP53.mutatesamples

OUTPUT FILE: BRCA.TP53.cna.intersected.output

# PANCAN ANALYSIS

In the directory, we also provide scripts that can merge all data from multiple cancer types. They all follow the same format, but we've separated them into three different scripts, one for MAF data, CNA data, and EXP data separately.  

As an optional endeavor, you can selectively merge from only certain cancertypes. For this we require you to submit a Matrix as a third argument. An example of one such matrix is included in the test file directory. The list should be by gene (in the first column). Each gene should have a row of 1s and 0s, 1s being for the cancer type from which you'd like to aggregate data from. Please follow the example submitted in the directory. The number of cancer types in the example is flexible, as long as they are in alphabetical order. 

USAGE:
 
     Rscript MERGE_MAF_PANCAN.R [InsertNameForDataFrame] [InsertNameForMergedCounts] [InsertOptionalMergeMatrix]
     Rscript MERGE_CNA_PANCAN.R [InsertNameForDataFrame] [InsertNameForMergedCounts] [InsertOptionalMergeMatrix]
     Rscript MERGE_EXP_PANCAN.R [InsertNameForDataFrame] [InsertNameForMergedCounts] [InsertOptionalMergeMatrix]

# CLASSIFICATION
Now that we have all of our features, we can go on with classification. The script CLASSY.R takes in maf counts, cna counts, and exp counts, and classifies all genes that have all features (no missing data) as a tumor supressor or oncogene. 

CLASSY.R allows you to provide your own training data or derives training data. Note that deriving training data from individual cancer types can be difficult due to limited data. We strongly suggest that if you're going to derive your training data, you do that from some form of PANCAN or multiple cancer type data, because our list of training genes might not be mutated in some cancers.

If you would like to derive your training data, the first argument should be a file that includes a list of genes and their classification. Classification should be binary (-1, and 1 ; or 0 and 1). You must also input a sixth argument "-d", to indicate that you would like to derive training data. 

If you would like to use the training data we've provided, the first argument should be the training data file. You should not have a sixth argument. 

The fifth argument is the name you'd like to give your classifcation.

USAGE:

     Rscript CLASSY.R trainingdata BRCA.TP53.maf.output BRCA.TP53.cna.intersected.output BRCA.TP53.rsem BRCA.classy
     Rscript CLASSY.R traininggenes PANCAN.maf PANCAN.cna.intersected.output PANCAN.rsem PANCAN.classy -d

INPUT: traininggenes (or trainingdata), a MAF counts file, a CNA counts file, and a RSEM counts file. 
OUTPUT: A merged data frame with all genes, their features, and classfications

# Miscellaneous Files

The miscellaneous directory contains files that were used to loop through multiple cancer type initially. We wish to present this as full transparency about all the scripts we use, but we do not intend this to be downloaded and used. 

#Questions: 

Please direct all questions to maheetha.b@gmail.com
