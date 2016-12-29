#!/gsc/bin/perl

use strict;
use warnings;

use FileHandle;
use Getopt::Long;
use IO::File;
use PDL::Stats::Basic;
use Statistics::Descriptive;
use File::Basename;
use Statistics::R;
use Cwd;

#################################################################################################
#################################################################################################
##########The purpose of this script is to extract various information from a MAF file###########
##########From the MAF file, we hope to count the total number of inactivation mutations#########
##########as well as the number of recurrent missense mutations amino-acid specific##############
##########as well as at the domain-specific level. After this script, each MAF file##############
##########should have a *.trunc, *.recur, and a *.missensedomains file with it ##################
##########The .recur and .missensedomains file will be analyzed mergerecur.py ###################
########## which will result in counts for these features.#######################################
#################################################################################################
#################################################################################################

##### SET UP USAGE ######
my $USAGE =<<USAGE;

     Usage:

         perl maf_script.pl --maf_file=maf --gene_index==0 --mut_index=8 --samp_index=9 --codon_index=15 --domain_index=40

     Description of options:

           maf_file = taken from TCGA data portal. It contains information on mutations found from a particular cancer type.
		
	   gene_index = the field [0 based indexing] that contains the gene name. Usually the first field. 

	   mut_index= = the field index that contains the type of mutation. (integer)

	   samp_index = the field index that contains the sample_id. In the format TCGA-XX-XXXX-10X-01X-XXXX-XX (integer)

	   codon_index = the field index that contains the codon. The codon information usually comes in the form p.[correctAA][AAnumber][MutatedAA] EX: p.A329P (integer)

	   domain_index = the field index that contains the domain information containing "|"-separated PRUXXXXX domains. (integer)
			
           help:  Prints out this helpful message

USAGE

#### User input for the indices

my $maf_file;
my $help =  0;
my $result;
my $gene_index;
my $mut_index;
my $sample_index;
my $codon_index;
my $domain_index;


$result = GetOptions(
	"maf_file=s" => \$maf_file,
	"gene_index=i" => \$gene_index,
	"mut_index=i" => \$mut_index,
	"sample_index=i" => \$sample_index,
	"codon_index=i" => \$codon_index,
	"domain_index=i" => \$domain_index,
	"help" => \$help,	
	);

if ($help){
   print "$USAGE\n";
   exit 0;	
}

my $output = IO::File->new("$maf_file.output", ">" ) or die "Could not create output file. $!";

#### Keeps track of total mutations ####
my %totmuts = ();

#### Keeps track of truncating mutations ####
my %trunc = ();

### Keeps track of recurrent mutations #####
my %recur = ();

#### Keeps track of recurrent domain mutations #####
my %domains = ();



#EXTRACTING APPROPRIATE DATA FROM THE MAF FILE

my $maf_fh = IO::File->new($maf_file) or die "Couldn't open $maf_file. $!";

while (my $maf = $maf_fh->getline) {
	
	#read line by line
	next if ($maf =~ m/^(#|Hugo_Symbol)/);
	my @array = split(/\t/, $maf);
	my ($gene,  $mutclass, $sample, $codon, $domain) = @array[$gene_index,$mut_index,$sample_index,$codon_index,$domain_index];

	#skip unnecessary types of mutations
	next if ($mutclass =~ m/^Silent/ || $gene =~ m/^(-|ENSG\d+|LOC\d+|OR\d+\w+\d+\w*|uc\d+.*)$/);
	next if($mutclass =~ m/^(Intron|RNA|3'Flank|3'UTR|5'Flank|5'UTR|IGR|Targeted_Region|Silent)$/);

	# extract the codon
	$codon = substr($codon, 0, -1); 
	
	# if mutclass is loss of function or truncating, store appropriately
	if ($mutclass =~ m/^(Splice_Site|Frame_Shift_Del|Frame_Shift_Ins|Nonstop_Mutation|Nonsense_Mutation)$/){ $trunc{$gene}++;}

	# store total number of mutations
	$totmuts{$gene}++;


	# if missense, note down the codon and domain
	if ($mutclass =~ m/^(Missense_Mutation)$/){	

		# store total missense and the codon
		
		$recur{$gene}{$codon}++;
		
		if ($domain ne ''){
			my @current = split(/:/, $domain);
			foreach my $dom (@current){
				next if ($dom !~ m/^PRU/);
				my $truncated = substr($dom, 0, 8);
				$domains{$gene}{$truncated}++;
			}
		}
	}
}
		

print $output "Gene\tTrunc\tRecur_Pos\tRecur_Dom\tTotal\n";


##### Retrieve stored information for each gene and output to a file ######
foreach my $gene ( keys %totmuts){
	my $count_trunc = 0;
	my $count_recur_aa = 0;
	my $count_recur_domain = 0;

	my $count_total = $totmuts{$gene};
	if (defined $trunc{$gene}){
		$count_trunc = $trunc{$gene};
	} 

	if (defined $recur{$gene}){
		my @codons = keys %{$recur{$gene}};
		foreach my $cod (@codons){
			$count_recur_aa += $recur{$gene}{$cod} - 1;		
		}
	} 

	if (defined $domains{$gene}){
		my @doms = keys %{$domains{$gene}};
		foreach my $dom (@doms){
			$count_recur_domain += $domains{$gene}{$dom} - 1;		
		}
	}
	
	
	print $output "$gene\t$count_trunc\t$count_recur_aa\t$count_recur_domain\t$count_total\n";
	
}
	



