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

   USAGE:
   
       bash get_maf_samples.sh -i BRCA.TP53.maf -s 12 -g 1 -m 9

SECOND COMMANDLINE EXECUTION:
bash mergerecur.sh -- if the files coming from createMAFinput.sh are named differently, you might want to edit the script

THIRD COMMANDLINE EXECUTION:
Rscript mergemafdata.R list_of_all_cancers_that_you_have_maf_data_for 

FOURTH COMMANDLINE EXECUTION:
bash filtermaf.sh



********FOR GENE EXPRESSION FILES********

FIRST COMMANDLINE EXECUTION
bash transform.sh 

SECOND COMMANLINE EXECUTION
bash transform2.sh 

THIRD COMMANDLINE EXECUTION
bash transpose3.sh 

FOURTH COMMANDLINE EXECUTION
bash geneExpression.sh

FIFTH COMMANDLINE EXECUTION
bash filtergeneexp.sh




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
