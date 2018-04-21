#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use sloan;

my $usage = 
"\nUsage: perl $0 [options/arguments]
   
   This script takes the trimmed fastq file from MDS_process_reads.pl and
   summarizes read families output in a tab delimited text file 
   (OUTPUT_NAME.read_fams.txt). Filters are applied and filtering metrics 
   and relevant sequences are reported in additional output files and STDOUT.
   
   REQUIRED ARGUMENTS
   
   Input File
   
         --input
         Trimmed fastq output from MDS_process_reads.pl


   OPTIONAL ARGUMENTS

   Output Name
   
         --output [default: fastq input name]
         Base name for all output files (additional extensions will be added)

   Minimum Agreement for Read Length within Read Families

         --min_len_agreement [default: 0.8]   
         Minimum proportion of reads with the same length within a read family.
         Families below this level will be excluded from further analysis.
 
   Suppress Output on Length Variation
    
    	--suppress_var_len_output
    	Add this flag to avoid outputing fasta file summarizing read families
    	with variation in read length.

   Minimum Base Quality
   
    	--min_base_qual [default = 20]
    	Change any basecalls with a quality score below this to N. These bases
    	will be excluded from summary calculations.

   Exclude Entire Reads that Contain Any Ns
   
    	--exclude_n_reads
    	Add this flag to remove the entire read if it contains any Ns.
    	         
   Illumina Phred Quality Score Version
    
    	--phred_offset [default: 33]
    	ASCII offset value for quality score encoding. 33 is standard for
    	current Illumina runs.


   EXAMPLE
   
         perl $0 \\
         	--exclude_n_reads \\
         	--input=Trimmed_Input.fq --output=my_output
       	
\n\n";

our $INPUT;
our $MIN_LEN_AGREEMENT = 0.8;
our $SUPPRESS_VAR_LEN_OUTPUT;
our $EXCLUDE_N_READS;
our $MIN_BASE_QUAL = 20;
our $PHRED_OFFSET = 33;
our $OUTPUT;


GetOptions(
    'input=s'  => \$INPUT,
	'min_len_agreement=i'  => \$MIN_LEN_AGREEMENT,
	'suppress_var_len_output'  => \$SUPPRESS_VAR_LEN_OUTPUT,
	'exclude_n_reads'  => \$EXCLUDE_N_READS,
	'min_base_qual=i'  => \$MIN_BASE_QUAL,
    'phred_offset=i'  => \$PHRED_OFFSET,
    'output=s'  => \$OUTPUT
);

#Check for proper inputs and return usage statement if there are errors
$INPUT or die ("\n$usage\n\nERROR: Must provide a fastq input filename with --input.\n\n");
$OUTPUT or $OUTPUT = $INPUT;
($MIN_LEN_AGREEMENT >= 0 and $MIN_LEN_AGREEMENT <= 1) or die ("\n$usage\n\nERROR: --min_len_agreement must be between from 0 to 1.\n\n");
(int($MIN_BASE_QUAL) == $MIN_BASE_QUAL and $MIN_BASE_QUAL >= 0) or die ("\n$usage\n\nERROR: min_base_qual must be a non-negative integer\.\n\n");


my $FH = open_file ($INPUT);

my %read_families_HoA;

my $n_read_count = 0;
my $total_read_count = 0;
my $reads_with_n_added_count = 0;
my $var_len_familes_count = 0;
my $tied_len_familes_count = 0;
my $zero_len_families_count = 0;
my $var_len_retained_count = 0;


print "\n" . (localtime) . "\nReading in sequences from $INPUT\, applying read filters, and grouping into read families.\n\n";


while (my $header = <$FH>){

	++$total_read_count;
	my $seq = <$FH>;
	my $drop_line = <$FH>;
	my $qual_line = <$FH>;
	chomp $header;
	chomp $seq;
	chomp $qual_line;
	
	my @sh = split (/\s/, $header);
	my $barcode = $sh[-1];
	
	my $n_added = 0;
	if ($MIN_BASE_QUAL){
		for (my $i = 0; $i < length ($seq); ++$i){
			ord(substr($qual_line, $i, 1)) - $PHRED_OFFSET < $MIN_BASE_QUAL and substr ($seq, $i, 1, "N") and $n_added = 1;
		}
	}
	$n_added and ++$reads_with_n_added_count;
	
	if ($seq =~ /N/){
		++$n_read_count;
		$EXCLUDE_N_READS and next;
	}
	
	push (@{$read_families_HoA{$barcode}}, $seq);
	
}


my %family_size_hash;
my %family_majority_read_length_hash;
my $total_read_families = 0;

my $FHO = open_output ("$OUTPUT\.read_fams.txt");
$SUPPRESS_VAR_LEN_OUTPUT or my $FHO_VARFAS = open_output ("$OUTPUT\.var_len_fams.fas");

print "\n" . (localtime) . "\nDeterming consensus sequences for each read family and filtering based on length variation within families.\n\n";

print $FHO "UMI\tFamily Size\tConsensus\tMatches\tMismatches\n";

foreach (sort keys %read_families_HoA){

	my @read_array = @{$read_families_HoA{$_}};
	my $fam_num = scalar(@read_array);
	my %read_seq_len_hash;
	foreach my $read (@read_array){
		++$read_seq_len_hash{length($read)};
	}
	
	my $majority = 0;
	my $majority_seq_len = 0;
	my $len_tie = 0;
	foreach my $read_length (keys %read_seq_len_hash){
		if ($read_seq_len_hash{$read_length} > $majority){
			$majority = $read_seq_len_hash{$read_length};
			$len_tie = 0;
			$majority_seq_len = $read_length;
		}elsif ($read_seq_len_hash{$read_length} == $majority){
			$len_tie = 1;
		}
	}

	$majority_seq_len or (++$zero_len_families_count and next);

	my $var_len_filter_fam = 0;
	if ($majority/$fam_num < $MIN_LEN_AGREEMENT){
		$var_len_filter_fam = 1;
		++$var_len_familes_count;
	}elsif($len_tie){
		$var_len_filter_fam = 1;
		++$tied_len_familes_count;
	}

	unless ($SUPPRESS_VAR_LEN_OUTPUT){
		if (scalar (keys %read_seq_len_hash) > 1){
			my $seq_within_var_fam_count = 0;
			foreach my $print_var_seq (@read_array){
				++$seq_within_var_fam_count;
				print $FHO_VARFAS ">$_\-$seq_within_var_fam_count\n$print_var_seq\n";
			}
			print $FHO_VARFAS "\n";
		}
	}

	$var_len_filter_fam and next;
	
	my @new_read_array;
	if (scalar (keys %read_seq_len_hash) > 1){
		++$var_len_retained_count;
		foreach my $read (@read_array){
			length($read) == $majority_seq_len and push (@new_read_array, $read);
		}
		@read_array = @new_read_array;
		$fam_num = scalar(@read_array);
	}

	my @consensus_counts_AoH;
	foreach my $read (@read_array){
		for (my $i = 0; $i < $majority_seq_len; ++$i){
			unless (substr ($read, $i, 1) eq "N"){
				++$consensus_counts_AoH[$i]{substr ($read, $i, 1)};
			}
		}
	}
		
		
	my $consensus_string;
	my $consensus_match;
	my $consensus_mismatch;
	for (my $i = 0; $i < $majority_seq_len; ++$i){
		if (exists $consensus_counts_AoH[$i]){
			my $href = $consensus_counts_AoH[$i];
			my $total = 0;
			my $max = 0;
			my $max_nuc;
				
			foreach my $nuc (sort keys %{$href}){
				$total += $href->{$nuc};
				if ($href->{$nuc} > $max){
					$max = $href->{$nuc};
					$max_nuc = $nuc;
				}elsif ($href->{$nuc} == $max){
					$max_nuc = "N";
				}
			}
			$consensus_string .= "$max_nuc";
			$consensus_match .= "$max\ ";
			my $minority = $total - $max;
			$consensus_mismatch .= "$minority\ ";	
		}else{
			$consensus_string .= "N";
			$consensus_match .= "0\ ";
			$consensus_mismatch .= "0\ ";
		}	
	}
	
	++$family_size_hash{$fam_num};
	++$family_majority_read_length_hash{$majority_seq_len};
	++$total_read_families;
	
	$consensus_match = substr($consensus_match, 0, -1);
	$consensus_mismatch = substr($consensus_mismatch, 0, -1);
	
	print $FHO "$_\t$fam_num\t$consensus_string\t$consensus_match\t$consensus_mismatch\n";
	
}

print "\n" . (localtime) . "\nAnalysis Complete\n\n";


print "$total_read_count\tTotal input read count\n";
$MIN_BASE_QUAL and print "$reads_with_n_added_count\tReads with at least one added N because base quality was below $MIN_BASE_QUAL\n";
print "$n_read_count\tTotal reads with at least one N";
if ($EXCLUDE_N_READS){
	print " (removed from analysis because --exclude_n_reads option was specified)\n";
}else{
	print " (retained in analysis because --exclude_n_reads option was not specified)\n";
}
$zero_len_families_count and print "$zero_len_families_count\tFamilies removed because all reads had length 0\n";
$MIN_LEN_AGREEMENT > 0 and print "$var_len_familes_count\tFamilies removed because reads with majority length were < $MIN_LEN_AGREEMENT of the read family\n";
$MIN_LEN_AGREEMENT < 1 and print "$var_len_retained_count\tFamilies retained because reads with majority length were >= $MIN_LEN_AGREEMENT of the read family\n";
$MIN_LEN_AGREEMENT <= 0.5 and print "$tied_len_familes_count\tFamilies removed because there were an equal number of reads with two or more different lengths within family\n";

print "\n$total_read_families\tTotal read families after filtering\n";

print "\nFamilySize\tCount\n";

foreach (sort {$a <=> $b} keys %family_size_hash){
	print "$_\t$family_size_hash{$_}\n";
}


print "\nMajoritySequenceLength\tFamilyCount\n";

foreach (sort {$a <=> $b} keys %family_majority_read_length_hash){
	print "$_\t$family_majority_read_length_hash{$_}\n";
}