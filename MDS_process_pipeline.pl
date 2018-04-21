#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use sloan;

my $usage = 
"\nUsage: perl $0 [options] --libraries_file=INPUT_FILE
   
   This is a wrapper script that calls a series for three scripts to process
   MDS datasets.
   
   REQUIRED ARGUMENTS
   
   File summarizing input libraries
   
         --libraries_file
         A tab-delimited text file that provides output name, R1 fastq file,
         marker info file, and optionally: R2 fastq file, marker name, filter
         seq file. See sample template (make sure new line characters can be 
         parsed by Perl in a UNIX environment).


   OPTIONAL ARGUMENTS (no defaults)

   File providing non-default options for each of the three processing scripts

         --options_file      
         A tab-delimited text file. First column is \"reads\", \"families\", 
         or \"variants\" to indicate the which script the option applies to.
         The second column is the option. See sample template.

   Scripts paths/names 

         --reads_script      
         Full path to MDS_process_reads.pl if not in working directory

         --families_script      
         Full path to MDS_process_families.pl if not in working directory

         --variants_script      
         Full path to MDS_process_variants.pl if not in working directory

   
   EXAMPLE
   
         perl $0 \\
         	--libraries_file=INPUT_FILE --options_file=OPTIONS_FILE \\
           	--reads_script=/PATH/TO/SCRIPTS/MDS_process_reads.pl \\
           	--families_script=/PATH/TO/SCRIPTS/MDS_process_families.pl \\
           	--variants_script=/PATH/TO/SCRIPTS/MDS_process_variants.pl
       	
\n\n";


our $LIBRARIES_FILE;
our $OPTIONS_FILE;
our $READS_SCRIPT = "MDS_process_reads.pl";
our $FAMILIES_SCRIPT = "MDS_process_families.pl";
our $VARIANTS_SCRIPT = "MDS_process_variants.pl";


GetOptions(
    'libraries_file=s'  => \$LIBRARIES_FILE,
	'options_file=s'  => \$OPTIONS_FILE,
	'reads_script=s'  => \$READS_SCRIPT,
	'families_script=s'  => \$FAMILIES_SCRIPT,
	'variants_script=s'  => \$VARIANTS_SCRIPT
);


#Check for proper inputs and return usage statement if there are errors
$LIBRARIES_FILE or die ("\n$usage\n\nERROR: Must provide an input file summarizing the libraries to process with --libraries_file.\n\n");
-e $READS_SCRIPT or die ("\n$usage\n\nERROR: Cannot find MDS_process_reads.pl. Provide full path to the script with --reads_script.\n\n");
-e $FAMILIES_SCRIPT or die ("\n$usage\n\nERROR: Cannot find MDS_process_families.pl. Provide full path to the script with --families_script.\n\n");
-e $VARIANTS_SCRIPT or die ("\n$usage\n\nERROR: Cannot find MDS_process_variants.pl. Provide full path to the script with --variants_script.\n\n");

my @reads_options;
my @families_options;
my @variants_options;

#parse options file
if ($OPTIONS_FILE){
	my @options = file_to_array ($OPTIONS_FILE);
	foreach (@options){
		chomp $_;
		my @sl =split(/\t/, $_);
		if ($sl[0] eq "reads"){
			push (@reads_options, $sl[1]);
		}elsif($sl[0] eq "families"){
			push (@families_options, $sl[1]);
		}elsif($sl[0] eq "variants"){
			push (@variants_options, $sl[1]);		
		}else{
			die ("\nERROR: could not parse the following line from $OPTIONS_FILE\. Must be tab delimited and the first field on each line must be \"reads\", \"families\", or \"variants\":\n\n$_\n\n");
		}
	}
}

my @libraries = file_to_array ($LIBRARIES_FILE);
my $first_line = shift @libraries;

#loop through lines in library input file. Generate command line strings for each of the three scripts and run them with system function.
foreach (@libraries){
	$_ =~ /^\s*$/ and next;
	chomp $_;
	my @sl = split (/\t/, $_);
	my $output = $sl[0];
	$output or die ("\nERROR: could not parse the following line from $LIBRARIES_FILE\. Must be tab delimited and the first field (output name) is mandatory:\n\n$_\n\n");
	my $reads1 = $sl[1];
	$reads1 or die ("\nERROR: could not parse the following line from $LIBRARIES_FILE\. Must be tab delimited and the second field (read1 fastq) is mandatory:\n\n$_\n\n");
	my $reads2;
	$sl[2] and $reads2 = $sl[2];
	my $marker_file = $sl[3];
	$marker_file or die ("\nERROR: could not parse the following line from $LIBRARIES_FILE\. Must be tab delimited and the fourth field (marker info file) is mandatory:\n\n$_\n\n");
	my $marker_name;
	$sl[4] and $marker_name = $sl[4];
	my $filter_file;
	$sl[5] and $filter_file = $sl[5];
	
	my $reads_string = "perl $READS_SCRIPT --r1_fastq=$reads1";
	$reads2 and $reads_string .= " --r2_fastq=$reads2";
	$reads_string .= " --marker_info_file=$marker_file";
	$marker_name and $reads_string .= " --marker_name=$marker_name";
	foreach my $ro (@reads_options){
		$reads_string .= " --$ro";
	}
	$reads_string .= " --output=$output";
	
	my $families_string = "perl $FAMILIES_SCRIPT --input=$output\.trimmed.fq";
	foreach my $fo (@families_options){
		$families_string .= " --$fo";
	}
	$families_string .= " --output=$output";

	my $variants_string = "perl $VARIANTS_SCRIPT --input=$output\.read_fams.txt --marker_info_file=$marker_file";
	$marker_name and $variants_string .= " --marker_name=$marker_name";
	$filter_file and $variants_string .= " --filter_seq_file=$filter_file";
	foreach my $vo (@variants_options){
		$variants_string .= " --$vo";
	}
	$variants_string .= " --output=$output";

	print "\n" . (localtime) . "\nRunning $READS_SCRIPT with the following command line:\n$reads_string\n\n";
	system ($reads_string);

	print "\n" . (localtime) . "\nRunning $FAMILIES_SCRIPT with the following command line:\n$families_string\n\n";
	system ($families_string);
	
	print "\n" . (localtime) . "\nRunning $VARIANTS_SCRIPT with the following command line:\n$variants_string\n\n";
	system ($variants_string);
	
}