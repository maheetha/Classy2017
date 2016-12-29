#!/usr/bin/perl

use strict;
use warnings;

use FileHandle;
use Getopt::Long;
use IO::File;


#############################################################################
#############################################################################
####### The purpose of this script is to parse the intersected ##############
####### copy number file, which has been intersected with a bed file. #######
####### The output file will contain genes and the number of samples ########
####### within various ranges of copy number values. ########################
#############################################################################
#############################################################################

# Set Up USAGE 

my $USAGE =<<USAGE;

     Usage:

         perl cna_script --cna_file=intersected_file --maf_sample=maf_sample_file

     Description of options:

           cna_file = this is a file that should ideally have 9 fields, 
		      that matches chromosome, and region, to a gene and sample. 
		      Make sure the file is grouped by gene. 

			FILE FORMAT:
			
			CHR	START	END	GENE	CHR	START	END	NUM_PROBES	Segment_Mean	Sample
			1	2985731	3355185	PRDM16	1	3218610	5323429	1374	-0.3028	TCGA-XX-XXXX-01A-11D-XXXX-01
			1	2985731	3355185	PRDM16	1	3218610	5738460	1678	0.0028	TCGA-XX-XXXX-11A-01D-XXXX-01
		     
	   maf_sample = this is a file that should have been output when running the scripts to parse the MAF data. 
			This is the output of the script createMAFinput, and should end with *.samples. It contains
			three fields: sample, gene, and type of mutation

			FILE FORMAT:

			Tumor_Sample	Hugo_Symbol	Variant_Classification
			TCGA-XX-XXXX	PTCHD2	Silent
			TCGA-XX-XXXX	VPS13D	Silent
			TCGA-XX-XXXX	PHC2	Silent
			TCGA-XX-XXXX	LCE1B	In_Frame_Del

           help:  Prints out this helpful message

USAGE

## Initialize all 
my $cna_file;
my $samp_file;
my $help = 0;
my $result;

### Initializing maf data hash ###
my %mafdata = ();

### Initializing options #####
$result = GetOptions(
	"cna_file=s" => \$cna_file,
	"maf_sample=s" => \$samp_file,	
	"help" => \$help,
);

if ($help){
   print "$USAGE\n";
   exit 0;	
}


##### Initializing all variables #####
my @t;
my @array;
my $chr;
my $gene;
my $sample;
my $ploidy;

#Initialize all variables that we will be using
my $original = "NULL";
my $mutated_samples = IO::File->new("$samp_file") or die "Could not open  file. $!";

# Store Mutated Samples for Each Gene, derived from MAF Sample File #
while (my $line = $mutated_samples->getline){
	next if ($line =~ m/^(#|Tumor_Sample)/);
	chomp($line);
        my @array = split(/\t/, $line);
        my ($sample, $gene, $mutation_type) = @array[0,1,2];
	if ($mutation_type =~m/^(In_Frame_Ins|Missense_Mutation|Frame_Shift_Del|In_Frame_Del|Frame_Shift_Ins|Nonsense_Mutation|Splice_Site)/){
		$mafdata{$gene}{$sample}++;
	}
}


### Initialize all variables of storage ####
my $intersect_fh = IO::File->new("$cna_file") or die "Could not open  file. $!";
my $cnaoutput = IO::File->new("$cna_file.output", ">" ) or die "Could not create output file. $!";
my $num_0_05 = 0;
my $num_05_1 = 0;
my $num_1_15 = 0;
my $num_15_2 = 0;
my $num_2_25 = 0;
my $num_25_3 = 0;
my $num_3_35 = 0;
my $num_35_4 = 0;
print $cnaoutput "Gene\t0\t0.5\t1\t1.5\t2\t2.5\t3\t3.5\n";


#######################################################################################################
##### For each line in the cna_intersected file, get the chromosome, gene, sample, and ploidy. ########
##### Calculate the number of samples in each ploidy range category, and store. #######################
#######################################################################################################

while (my $line = $intersect_fh->getline) {

	#### store the line in an array, and extract the gene, sample, and ploidy
	chomp($line);
	@array = split(/\t/, $line);
	($gene, $sample, $ploidy) = @array[3,9,8];
	
	# confirming that this particular sample is mutated for this gene, else we move onto the next sample
	if (exists $mafdata{$gene}){
	    my @keys = keys %{$mafdata{$gene}};
	    my $cut_sample = substr($sample, 0, 12);
	    if ( grep( /^$cut_sample$/, @keys ) ) {
		 } else { next;}
	   } else { next;}       

	#If the current gene isn't what it was previously, then we know we're done with a gene, and we have to update the numbers we have already stored for it.
	if ($gene ne $original){
		
		#store data
		if ($original ne 'NULL'){

			print $cnaoutput "$gene\t$num_0_05\t$num_05_1\t$num_1_15\t$num_15_2\t$num_2_25\t$num_25_3\t$num_3_35\t$num_35_4\n";
		}

		#re initialize
		$original = $gene;
		$num_0_05 = 0;
		$num_05_1 = 0;
		$num_1_15 = 0;
		$num_15_2 = 0;
		$num_2_25 = 0;
		$num_25_3 = 0;
		$num_3_35 = 0;
		$num_35_4 = 0;
	}

	#############################################################################
	### However, if the gene is the same as it was before, we can keep adding ###
	### to the numbers we have now. We just add to the ploidy ranges we have ####
	#############################################################################
	
	$original = $gene;
	$ploidy = exp($ploidy*log(2.0)) * 2;
		if ($ploidy <= 0.5){
			$num_0_05++;
		}
		if ($ploidy > 0.5 and $ploidy <= 1.0){
			$num_05_1++;
		}
		if ($ploidy > 1.0 and $ploidy <= 1.5){
                        $num_1_15++;
                }
		if ($ploidy > 1.5 and $ploidy <= 2.0){
                        $num_15_2++;
                }
		if ($ploidy > 2.0 and $ploidy <= 2.5){
                        $num_2_25++;
                }
		if ($ploidy > 2.5 and $ploidy <= 3.0){
                        $num_25_3++;
                }
		if ($ploidy > 3.0 and $ploidy <= 3.5){
                        $num_3_35++;
                }
		if ($ploidy > 3.5 and $ploidy <= 4.0){
                        $num_35_4++;
                }
}








