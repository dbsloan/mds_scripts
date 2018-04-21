#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use sloan;

my $usage = 
"\nUsage: perl $0 [options/arguments]
   
   This script takes the read fams text file from MDS_process_families.pl and
   summarizes variant calls.
   
   REQUIRED ARGUMENTS
   
   Input File
   
         --input
         Read fams output text file from MDS_process_families.pl

   MDS Marker Information File
   
         --marker_info_file
         A tab-delimited text file that reports marker names, primers, and 
         reference sequences


   OPTIONAL ARGUMENTS

   Output Name
   
         --output [default: input file name]
         Base name for all output files (additional extensions will be added)

   MDS Marker Name
   
         --marker_name [default: first name in Marker Information File]
         Name of MDS marker found within Marker Information File

   Minimum Agreement for Basecalls within Read Families

         --min_base_agreement [default: 0.8]   
         Minimum proportion of reads with the same basecall within a read family.
         Bases below this level will be excluded from further analysis.
 
   Minimum Read Family Size to Output Length Variants to Fasta File
         
         --min_len_fasta [default: 3]
         Read families with consensus length different from reference will be
         output to fasta file if they have a family size larger than this value.   

   Minimum Read Family Size to Output SNP Variants to Fasta File
         
         --min_snp_fasta [default: 3]
         Read families with SNP variant calls will be output to fasta file if 
         they have a family size larger than this value.   

   Minimum Read Family Size to Output MultiSNP Variants to Fasta File
         
         --min_multisnp_fasta [default: 2]
         Read families with more than one SNP variant (but with the expected reference
         length) will be output to a fasta file if they have at least this many 
         SNPs. This file should be checked for non-specific amplicons that happen
         to be the correct length (can greatly inflate variant call rates).
         Note that this option does not affect the tabulation of counts in the \"single\"
         output files. These treat any sequence with >1 SNPs as multivariant.

   Suppress Fasta Output Files
    
    	--suppress_fasta_output
    	Add this flag to avoid outputing fasta files summarizing with variants seqs.

   Exclude Read Families with Any Ns in Consensus Sequence
    
    	--exclude_n_families
    	Add this flag to exclude the entire sequence from any read families that have
    	an N at any position in their consensus sequences.

   Filter Known Contaminant Sequences
    
    	--filter_seq_file
    	Simple text file with one sequence per line (no headers or other content).
    	Read families with consensus sequences that are exact, full-length matches to any
    	of the sequences in this file will be excluded from statistical analysis. This 
    	can be used to filter known contaminants like NUMTs.

   EXAMPLE
   
         perl $0 \\
         	--marker_info_file=MDS_markers.txt --marker_name=At_mt1 \\
         	--min_base_agreement=1 \\
         	--input=INPUT_READ_FAMS_FILE.txt --output=my_output
       	
\n\n";


our $INPUT;
our $MIN_BASE_AGREEMENT = 0.8;
our $MIN_LEN_FASTA = 3;
our $MIN_SNP_FASTA = 3;
our $MIN_MULTISNP_FASTA = 2;
our $SUPPRESS_FASTA_OUTPUT;
our $EXCLUDE_N_FAMILIES;
our $MARKER_INFO_FILE;
our $MARKER_NAME;
our $FILTER_SEQ_FILE;
our $OUTPUT;


GetOptions(
    'INPUT=s'  => \$INPUT,
	'min_base_agreement=i'  => \$MIN_BASE_AGREEMENT,
	'min_len_fasta=i'  => \$MIN_LEN_FASTA,
	'min_snp_fasta=i'  => \$MIN_SNP_FASTA,
	'min_multisnp_fasta=i'  => \$MIN_MULTISNP_FASTA,
	'suppress_fasta_output'  => \$SUPPRESS_FASTA_OUTPUT,
	'exclude_n_families'  => \$EXCLUDE_N_FAMILIES,
    'marker_info_file=s'  => \$MARKER_INFO_FILE,
    'marker_name=s'  => \$MARKER_NAME,
    'filter_seq_file=s'  => \$FILTER_SEQ_FILE,
    'output=s'  => \$OUTPUT
);



#Check for proper inputs and return usage statement if there are errors
$INPUT or die ("\n$usage\n\nERROR: Must provide a read family input filename (from MDS_process_families.pl) with --input.\n\n");
$OUTPUT or $OUTPUT = $INPUT;
($MIN_BASE_AGREEMENT >= 0 and $MIN_BASE_AGREEMENT <= 1) or die ("\n$usage\n\nERROR: --min_base_agreement must be between from 0 to 1.\n\n");
(int($MIN_LEN_FASTA) == $MIN_LEN_FASTA and $MIN_LEN_FASTA >= 0) or die ("\n$usage\n\nERROR: --min_len_fasta must be a non-negative integer.\n\n");
(int($MIN_SNP_FASTA) == $MIN_SNP_FASTA and $MIN_SNP_FASTA >= 0) or die ("\n$usage\n\nERROR: --min_snp_fasta must be a non-negative integer.\n\n");
(int($MIN_MULTISNP_FASTA) == $MIN_MULTISNP_FASTA and $MIN_MULTISNP_FASTA >= 0) or die ("\n$usage\n\nERROR: --min_multisnp_fasta must be a non-negative integer.\n\n");


#Make sure information is provided for name of marker and reference sequence
my @marker_info;
my $first_marker_name;
my %markers_HoA;
$MARKER_INFO_FILE or die ("\n$usage\n\nERROR: A file must be specified with --marker_info_file that provides marker name, primers, and expected sequence.\n\n");

#Read in information in $MARKER_INFO_FILE. Record first marker name and store everything in hash of arrays
@marker_info = file_to_array ($MARKER_INFO_FILE);
my @split_marker_line = split(/\t/, $marker_info[0]);
$first_marker_name = $split_marker_line[0];
foreach (@marker_info){
	chomp $_;
	my @sl = split (/\t/, $_);
	my $name = shift (@sl);
	$markers_HoA{$name}=\@sl;
}

#Check if marker name was provided. If not set it to the name of the first marker in $MARKER_INFO_FILE
unless ($MARKER_NAME){
	$MARKER_NAME = $first_marker_name;
	print STDERR "\nWARNING: No marker name was provided. Assuming the marker is the first one listed in $MARKER_INFO_FILE ($MARKER_NAME\). If this is not correct, specify actual marker name with --marker_name.\n\n";
}
exists ($markers_HoA{$MARKER_NAME}) or die ("\n$usage\n\nERROR: Could not find the marker name $MARKER_NAME in $MARKER_INFO_FILE\.\n\n");

my $ref_seq = $markers_HoA{$MARKER_NAME}[2];
my $ref_len = length($ref_seq);

#if a file of sequences to filter out was provided save those in a hash to apply below
my %filter_seq_hash;
if ($FILTER_SEQ_FILE){
	my @bad_seqs = file_to_array ($FILTER_SEQ_FILE);
	foreach (@bad_seqs){
		chomp $_;
		$filter_seq_hash{$_} = 1;
	}
}


my $FH = open_file ($INPUT);

my $first_line = <$FH>;

my $input_fam_count = 0;
my %length_mismatch_count_hash;
my %N_count;
my $N_fam_excludes = 0;
my $total_len_mismatch_count = 0;
my %hetero_count_HoA;
my @variant_AoHoH;
my @variant_single_AoHoH;
my %hetero_above_threshold_HoA;
my %fasta_snp_hash;
my %fasta_multisnp_hash;
my %fasta_indel_hash;
my %counts_by_size_HoA;
my %multisnp_hash;
my %filter_seq_counts;
my $total_filter_seq_count = 0;

print "\n" . (localtime) . "\nAnalyzing read families in $INPUT to identify variants relative to the $MARKER_NAME reference sequence\n\n";

#Loop through each data line in the read family file
while (<$FH>){
	chomp $_;
	my @sl = split (/\t/, $_);
	
	++$input_fam_count;
	
	##If option is specified, exclude any families with a consensus seq that contains N
	if ($EXCLUDE_N_FAMILIES and $sl[2] =~ /N/){
		++$N_fam_excludes;
		next;
	}
	
	#if a file of "bad" (contaminating) sequences was provided, exclude families that have one of those exact sequences
	if ($FILTER_SEQ_FILE and exists ($filter_seq_hash{$sl[2]})){
		++$filter_seq_counts{$sl[2]};
		++$total_filter_seq_count;
		next;
	}
	
	##Record consensus seqs that do not match reference length and exclude from further SNP processing
	unless (length ($sl[2]) == $ref_len){ 
		++$length_mismatch_count_hash{$sl[1]};
		++$total_len_mismatch_count;
		unless ($SUPPRESS_FASTA_OUTPUT or $sl[1] < $MIN_LEN_FASTA){
			$fasta_indel_hash{"$sl[0] $sl[1]"} = $sl[2];
		}
		next;
	}
	
	#split the strings with match and mismatch counts by site.
	my @sl_match = split (/\ /, $sl[3]);
	my @sl_mismatch = split (/\ /, $sl[4]);
	
	scalar (@sl_match) == scalar (@sl_mismatch) or die ("\nERROR: match and mismatch strings do not have same number of scores\n\n");
	scalar (@sl_match) == $ref_len or die ("\nERROR: match and mismatch strings do not match reference length\n\n");

	
	my $fasta_snp_out = 0;
	my $read_snp_count = 0;
	my $solo_snp_fam_size = 0;
	
	#use this data structure to store info on all sites in read. Only store in final single_SNP dataset if read is not multiSNP
	my @read_variant_AoHoH;

	
	#loop over each site in sequence
	for (my $i=0; $i < $ref_len; ++$i){
		
		my $total_count = $sl_match[$i] + $sl_mismatch[$i];
		
		#record Ns and skip to next site
		my $nuc = substr($sl[2], $i, 1);
		if ($nuc eq "N"){
			++$N_count{$sl[1]};
			next;
		} 
		my $ref_nuc = substr($ref_seq, $i, 1);
		
		#check for disagreement within read family and check if agreement is below specified thereshold
		if ($sl_match[$i] / $total_count < $MIN_BASE_AGREEMENT){
			if ($nuc eq $ref_nuc){
				++$hetero_count_HoA{$total_count}[0];
			}else{
				++$hetero_count_HoA{$total_count}[1];
			}
			next;		
		}
		
		#Count up variants or matches in a two-nucleotide format (AA, AC, GG, etc.)
		++$variant_AoHoH[$i]->{"$ref_nuc$nuc"}->{$total_count};
		
		#add values to the single dataset for now. If the read ends up as a multiSNP, sustract them later (should be rare so I think this is faster)
		++$variant_single_AoHoH[$i]->{"$ref_nuc$nuc"}->{$total_count};
		
		#will use this data structure to substract if multisnp
		++$read_variant_AoHoH[$i]->{"$ref_nuc$nuc"}->{$total_count};
		
		#Record SNPs
		if ($ref_nuc ne $nuc){
			++$read_snp_count;
			++$counts_by_size_HoA{$total_count}[1];
			$solo_snp_fam_size = $total_count;
			unless ($SUPPRESS_FASTA_OUTPUT or $total_count < $MIN_LEN_FASTA){
				$fasta_snp_out = 1;
			}
		}else{
			++$counts_by_size_HoA{$total_count}[0];
		}
		
		#Record cases where agreement for site was not perfect but above min threshold
		unless($sl_match[$i] / $total_count == 1){
			if ($nuc eq $ref_nuc){
				++$hetero_above_threshold_HoA{$total_count}[0];
			}else{
				++$hetero_above_threshold_HoA{$total_count}[1];
			}
		}
	}
	
	if ($read_snp_count > 1){ 
		++$multisnp_hash{$sl[1]};
		#subtract counts form singleSNPs because this was a multiSNP read
		for (my $nuc_array_pos = 0; $nuc_array_pos < scalar(@read_variant_AoHoH); ++$nuc_array_pos){
			$read_variant_AoHoH[$nuc_array_pos] or next;
			foreach my $subtype_key (keys %{$read_variant_AoHoH[$nuc_array_pos]}){
				foreach my $famsize_key (keys %{$read_variant_AoHoH[$nuc_array_pos]->{$subtype_key}}){
					--$variant_single_AoHoH[$nuc_array_pos]->{$subtype_key}->{$famsize_key};
				}
			} 
		}
	}
	$read_snp_count == 1 and ++$counts_by_size_HoA{$solo_snp_fam_size}[2];
	
	unless ($SUPPRESS_FASTA_OUTPUT){
		$fasta_snp_out and $fasta_snp_hash{"$sl[0] $sl[1]"} = $sl[2];
		$read_snp_count > 1 and $sl[1] >= $MIN_MULTISNP_FASTA and $fasta_multisnp_hash{"$sl[0] $sl[1]"} = $sl[2];
	}
}


print "\n" . (localtime) . "\nAnalysis Complete; Printing Output\n\n";

print "$input_fam_count\tRead families analyzed\n";
$EXCLUDE_N_FAMILIES and print "$N_fam_excludes\tRead families excluded because consensus contains one or more Ns and --exclude_n_families was specified\n";
print "$total_len_mismatch_count\tRead families excluded because consensus length did not match reference sequence\n\n";
$FILTER_SEQ_FILE and print "$total_filter_seq_count\tRead families excluded because they match sequences provided with --filter_seq_file in $FILTER_SEQ_FILE\n";
if ($total_filter_seq_count){
	foreach (sort keys %filter_seq_counts){
		print "$filter_seq_counts{$_}\t$_\n";
	}
}

my $FHO_VARSUMM = open_output("$OUTPUT\.var_summ.txt");

print $FHO_VARSUMM "FamilySize\tMatch\tMismatch\tMismatch_singles\tMultiSNP_Fams\tLengthDiffs\tN_Count\tHetero_BelowThresh_MajorityRef\tHetero_BelowThresh_MinorityRef\tHetero_AboveThresh_MajorityRef\tHetero_AboveThresh_MinorityRef\n";

#Fill in zeros for undefined values in data structures
foreach (sort {$a <=> $b} keys %counts_by_size_HoA){
	#set all undefined values to 0;
	unless ($counts_by_size_HoA{$_}[0]){
		$counts_by_size_HoA{$_}[0] = 0;
	}
	unless ($counts_by_size_HoA{$_}[1]){
		$counts_by_size_HoA{$_}[1] = 0;
	}
	unless ($counts_by_size_HoA{$_}[2]){
		$counts_by_size_HoA{$_}[2] = 0;
	}
	unless($multisnp_hash{$_}){
		$multisnp_hash{$_} = 0;
	}
	unless ($length_mismatch_count_hash{$_}){
		$length_mismatch_count_hash{$_} = 0;
	}
	unless ($N_count{$_}){
		$N_count{$_}=0;
	}
	unless ($hetero_count_HoA{$_}[0]){
		$hetero_count_HoA{$_}[0] = 0;
	}
	unless ( $hetero_count_HoA{$_}[1]){
		$hetero_count_HoA{$_}[1] = 0;
	}
	unless ($hetero_above_threshold_HoA{$_}[0]){
		$hetero_above_threshold_HoA{$_}[0] = 0;
	}
	unless ($hetero_above_threshold_HoA{$_}[1]){
		$hetero_above_threshold_HoA{$_}[1] = 0;
	}	

	print $FHO_VARSUMM "$_\t$counts_by_size_HoA{$_}[0]\t$counts_by_size_HoA{$_}[1]\t$counts_by_size_HoA{$_}[2]\t";	
	print $FHO_VARSUMM "$multisnp_hash{$_}\t$length_mismatch_count_hash{$_}\t$N_count{$_}\t";	
	print $FHO_VARSUMM "$hetero_count_HoA{$_}[0]\t$hetero_count_HoA{$_}[1]\t$hetero_above_threshold_HoA{$_}[0]\t$hetero_above_threshold_HoA{$_}[1]\n";	
}


unless ($SUPPRESS_FASTA_OUTPUT){
	
	my $FH_SNP_FAS = open_output("$OUTPUT\.snp_family.fas");
	my $FH_INDEL_FAS = open_output("$OUTPUT\.indel_family.fas");
	my $FH_MULTISNP_FAS = open_output("$OUTPUT\.multisnp_family.fas");
	
	foreach (sort keys %fasta_snp_hash){
		print $FH_SNP_FAS ">$_\n$fasta_snp_hash{$_}\n"; 
	}
	foreach (sort keys %fasta_indel_hash){
		print $FH_INDEL_FAS ">$_\n$fasta_indel_hash{$_}\n"; 
	}
	foreach (sort keys %fasta_multisnp_hash){
		print $FH_MULTISNP_FAS ">$_\n$fasta_multisnp_hash{$_}\n"; 
	}
}

my $FH_VAR = open_output("$OUTPUT\.var.txt");

print $FH_VAR "Position\tType\tReference\tRead\tFamilySize\tCount\n";

foreach (my $pos = 0; $pos < scalar (@variant_AoHoH); ++$pos){
	
	my $HoH_ref = $variant_AoHoH[$pos];
	
	foreach my $dinuc (sort keys %{$HoH_ref}){
		my $split_dinuc_ref = substr ($dinuc, 0, 1);
		my $split_dinuc_alt = substr ($dinuc, 1, 1);
		my $match = "Mismatch";
		if ($split_dinuc_ref eq $split_dinuc_alt){
			$match = "Match";
		}
		foreach my $fam_size (sort {$a <=> $b} keys %{$HoH_ref->{$dinuc}}){
			print $FH_VAR $pos + 1, "\t$match\t$split_dinuc_ref\t$split_dinuc_alt\t$fam_size\t$HoH_ref->{$dinuc}->{$fam_size}\n";
		}
	}	
}


my $FH_VAR_SINGLE = open_output("$OUTPUT\.var.single.txt");

print $FH_VAR_SINGLE "Position\tType\tReference\tRead\tFamilySize\tCount\n";

foreach (my $pos = 0; $pos < scalar (@variant_single_AoHoH); ++$pos){
	
	my $HoH_ref = $variant_single_AoHoH[$pos];
	
	foreach my $dinuc (sort keys %{$HoH_ref}){
		my $split_dinuc_ref = substr ($dinuc, 0, 1);
		my $split_dinuc_alt = substr ($dinuc, 1, 1);
		my $match = "Mismatch";
		if ($split_dinuc_ref eq $split_dinuc_alt){
			$match = "Match";
		}
		foreach my $fam_size (sort {$a <=> $b} keys %{$HoH_ref->{$dinuc}}){
			print $FH_VAR_SINGLE $pos + 1, "\t$match\t$split_dinuc_ref\t$split_dinuc_alt\t$fam_size\t$HoH_ref->{$dinuc}->{$fam_size}\n";
		}
	}	
}