#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use IPC::Cmd qw(can_run);
use sloan;

my $usage = 
"\nUsage: perl $0 [options/arguments]
   
   This script takes raw MDS Illumina fastq input files and returns a primer-
   trimmed fastq with the MDS barcode added to the header of each read 
   (OUTPUT_NAME.trimmed.fq). A number of filters are applied and filtering 
   metrics are reported in additional output files and STDOUT.
   
   REQUIRED ARGUMENTS
   
   R1 Fastq File
   
         --r1_fastq
         File containing Illumina read 1 MDS sequences (raw sequences). Can 
         be gzipped.

   MDS Marker Information File
   
         --marker_info_file
         A tab-delimited text file that reports marker names, primers, and 
         reference sequences


   OPTIONAL ARGUMENTS

   R2 Fastq File
   
         --r2_fastq
         File containing Illumina read 2 MDS sequences (raw sequences) if
         paired-end sequencing was performed. If provided, bbmerge.sh will
         be run to collapse fwd and rev reads. Can be gzipped.

   MDS Marker Name
   
         --marker_name [default: first name in Marker Information File]
         Name of MDS marker found within Marker Information File

   Output Name
   
         --output [default: r1_fastq file name]
         Base name for all output files (additional extensions will be added)

   BBmerge Executable

         --bbmerge_exe [default: bbmerge.sh]      
         Full path to bbmerge.sh if not in default PATH

   BBmerge Quality Score for Trimming

         --bbmerge_trimq [default: 20]   
         Quality score used for trimq parameter in bbmerge.sh
         
   BBmerge Minimum Overlap

         --bbmerge_minoverlap [default: 30]   
         Sequence length for minoverlap parameter in bbmerge.sh

   BBmerge Mismatches

         --bbmerge_mismatches [default: 5]   
         Maximum number of mismathces for mismatches parameter in bbmerge.sh

   Turn Off Cutadapt Trimming

         --disable_cutadapt 
         Add this flag to turn off primer trimming with cutadapt.

   Cutadapt Executable

         --cutadapt_exe [default: cutadapt]      
         Full path to cutadapt if not in default PATH

   Cutadapt Error Tolerance

         --cutadapt_err_tol [default: 0.15]   
         Error tolerance for cutadapt trimming

   Cutadapt Cores

         --cutadapt_cores [default: 4]   
         Number of cores to request for cutadapt step. This may need to be set
         to 1 if using a cutadapt version with Pyhton v2.

   Cutadapt Minimum Length

         --cutadapt_min_len [default: 1]   
         Value for minimum_length parameter in cutadapt

   Report Untrimmed Cutadapt Reads

         --cutadapt_report_untrimmed   
         Add this flag to export a fastq file with reads that are not primer-
         trimmed by cutadapt (regardless, these reads are not used in
         subsequent analysis steps).

   MDS Random Barcode Length

         --barcode_len [default: 14] 
         Number of Ns in the random barcode (unique molecular identifier) at
         beginning of each MDS reads.

   Disable Repetitive Barcode Filter
    
    	--disable_rep_filter
    	Add this flag to turn off the default filter that excludes reads with
    	MDS barcodes that are just one long homopolymer (e.g. AAAAAAAAAAAAAA).

   Filtering Optical Duplicates
    
    	--min_optical_dist [default: 0]
    	Minimum pixel distance for retaining a read with the same barcode as 
    	another read on the same tile within the same Illumina flow cell.
    	Specifying non-zero values here will help eliminate optical/clustering
    	duplicates (MDS ideally only compares amplification duplicates).
    	
   Filtering Tile Edge Reads
    
    	--min_tile_edge_dist [default: 0]
    	Minimum pixel distance for X and Y coordinates on the Illumina tile.
    	Specifying a non-zero value will help eliminate tile-edge duplicates
    	(MDS ideally only compares amplification duplicates). Note that this 
    	option will discard lots of non-duplicate reads. It just discards all
    	reads within that distance of left/bottom tile edge.
  	
   Exclude Reads with Low Barcode Quality
    
    	--min_barcode_qual [default: 20]
    	Discard reads that have basecalls in barcode with quality value lower
    	than this.

   Illumina Phred Quality Score Version
    
    	--phred_offset [default: 33]
    	ASCII offset value for quality score encoding. 33 is standard for
    	current Illumina runs.

   Disable Filtering of Barcodes with Ns
    
    	--disable_n_filter
    	Add this flag to turn off the default filter that excludes reads with
    	MDS barcodes that contain an N.

   Delete Fastq Files
    
    	--delete_intermediate_fastqs
    	Add this flag to delete the large intermediate fastq files and only
    	the final merged/trimmed file.


   EXAMPLE
   
         perl $0 \\
         	--r1_fastq=READ1.fastq.gz -r2_fastq=READ2.fastq.gz \\
         	--marker_info_file=MDS_markers.txt --marker_name=At_mt1 \\
         	--min_optical_dist=5000 --min_tile_edge_dist=2000 \\
         	--delete_intermediate_fastqs \\
         	--output=my_output
       	
\n\n";

our $R1_FASTQ;
our $R2_FASTQ;
our $BBMERGE_FASTQ;
our $BBMERGE_EXE = "bbmerge.sh";
our $BBMERGE_TRIMQ = 20;
our $BBMERGE_MINOVERLAP = 30;
our $BBMERGE_MISMATCHES = 5;
our $DISABLE_CUTADAPT;
our $CUTADAPT_EXE = "cutadapt";
our $CUTADAPT_ERR_TOL = 0.15;
our $CUTADAPT_CORES = 4;
our $CUTADAPT_MIN_LEN = 1;
our $CUTADAPT_REPORT_UNTRIMMED;
our $BARCODE_LEN = 14;
our $DISABLE_REP_FILTER;
our $MIN_OPTICAL_DIST = 0;
our $MIN_TILE_EDGE_DIST = 0;
our $MIN_BAROCODE_QUAL = 20;
our $PHRED_OFFSET = 33;
our $DISABLE_N_FILTER;
our $MARKER_INFO_FILE;
our $MARKER_NAME;
our $DELETE_INTERMEDIATE_FASTQS;
our $OUTPUT;

GetOptions(
    'r1_fastq=s'  => \$R1_FASTQ,
    'r2_fastq=s'  => \$R2_FASTQ,
    'bbmerge_exe=s'  => \$BBMERGE_EXE,
    'bbmerge_trimq=i'  => \$BBMERGE_TRIMQ,
    'bbmerge_minoverlap=i'  => \$BBMERGE_MINOVERLAP,
    'bbmerge_mismatches=i'  => \$BBMERGE_MISMATCHES,
    'disable_cutadapt'  => \$DISABLE_CUTADAPT,
    'cutadapt_exe=s'  => \$CUTADAPT_EXE,
    'cutadapt_err_tol=i'  => \$CUTADAPT_ERR_TOL,
    'cutadapt_cores=i'  => \$CUTADAPT_CORES,
    'cutadapt_min_len=i'  => \$CUTADAPT_MIN_LEN,
    'cutadapt_report_untrimmed'  => \$CUTADAPT_REPORT_UNTRIMMED,
    'barcode_len=i'  => \$BARCODE_LEN,
    'disable_rep_filter'  => \$DISABLE_REP_FILTER,
	'min_optical_dist=i'  => \$MIN_OPTICAL_DIST,
    'min_tile_edge_dist=i'  => \$MIN_TILE_EDGE_DIST,
    'min_barcode_qual=i'  => \$MIN_BAROCODE_QUAL,
    'phred_offset=i'  => \$PHRED_OFFSET,
    'disable_n_filter'  => \$DISABLE_N_FILTER,
    'marker_info_file=s'  => \$MARKER_INFO_FILE,
    'marker_name=s'  => \$MARKER_NAME,
    'delete_intermediate_fastqs'  => \$DELETE_INTERMEDIATE_FASTQS,
    'output=s'  => \$OUTPUT
);

#Check for proper inputs and return usage statement if there are errors
(int($BBMERGE_MINOVERLAP) == $BBMERGE_MINOVERLAP and $BBMERGE_MINOVERLAP >= 0) or die ("\n$usage\n\nERROR: bbmerge_minoverlap must be a non-negative integer\.\n\n");
(int($BBMERGE_TRIMQ) == $BBMERGE_TRIMQ and $BBMERGE_TRIMQ >= 0) or die ("\n$usage\n\nERROR: bbmerge_trimq must be a non-negative integer\.\n\n");
(int($BBMERGE_MISMATCHES) == $BBMERGE_MISMATCHES and $BBMERGE_MISMATCHES >= 0) or die ("\n$usage\n\nERROR: bbmerge_mismatches must be a non-negative integer\.\n\n");
(int($BARCODE_LEN) == $BARCODE_LEN and $BARCODE_LEN > 0) or die ("\n$usage\n\nERROR: barcode_len must be a positive integer\.\n\n");
(int($MIN_OPTICAL_DIST) == $MIN_OPTICAL_DIST and $MIN_OPTICAL_DIST >= 0) or die ("\n$usage\n\nERROR: min_optical_dist must be a non-negative integer\.\n\n");
(int($MIN_TILE_EDGE_DIST) == $MIN_TILE_EDGE_DIST and $MIN_TILE_EDGE_DIST >= 0) or die ("\n$usage\n\nERROR: min_tile_edge_dist must be a non-negative integer\.\n\n");
(int($MIN_BAROCODE_QUAL) == $MIN_BAROCODE_QUAL and $MIN_BAROCODE_QUAL >= 0) or die ("\n$usage\n\nERROR: min_barcode_qual must be a non-negative integer\.\n\n");
(int($PHRED_OFFSET) == $PHRED_OFFSET and $PHRED_OFFSET >= 0) or die ("\n$usage\n\nERROR: phred_offset must be a non-negative integer\.\n\n");
($CUTADAPT_ERR_TOL >= 0 and $CUTADAPT_ERR_TOL <= 1) or die ("\n$usage\n\nERROR: cutadapt_err_tol must be a number from 0 to 1\.\n\n");
(int($CUTADAPT_CORES) == $CUTADAPT_CORES and $CUTADAPT_CORES >= 1) or die ("\n$usage\n\nERROR: cutadapt_cores must be a positive integer\.\n\n");
(int($CUTADAPT_MIN_LEN) == $CUTADAPT_MIN_LEN and $CUTADAPT_MIN_LEN >= 0) or die ("\n$usage\n\nERROR: cutadapt_min_len must be a non-negative integer\.\n\n");

if ($CUTADAPT_REPORT_UNTRIMMED and $CUTADAPT_CORES > 1){
	$CUTADAPT_CORES = 1;
	print STDERR "\nWARNING: Changing number of cores for cutadapt to 1. Multithreading not supported when outputting untrimmed reads to separate file. To allow multithreading do not use --cutadapt_report_untrimmed and set number of cores with --cutadapt_cores.\n\n";
}

unless ($DISABLE_CUTADAPT){
	can_run ($CUTADAPT_EXE) or die ("\n$usage\n\nERROR: cutadapt is enabled, but could not find following executable in PATH: $CUTADAPT_EXE. Provide a path and executable name with --cutadapt_exe or disable with --disable_cutadapt.\n\n");
}
$R1_FASTQ or die ("\n$usage\n\nERROR: Must provide a fastq input filename with --r1_fastq.\n\n");
-e $R1_FASTQ or die ("\n$usage\n\nERROR: File specified with --r1_fastq does not exist: $R1_FASTQ.\n\n");
$OUTPUT or $OUTPUT = $R1_FASTQ;

#Make sure information is provided for name of marker and primers (only needed if cutadapt will be run)
my @marker_info;
my $first_marker_name;
my %markers_HoA;
unless ($DISABLE_CUTADAPT){
	$MARKER_INFO_FILE or die ("\n$usage\n\nERROR: If cutadapt will be run, a file must be specified with --marker_info_file that provides marker name, primers, and expected sequence. Cutadapt can be disabled with --disable_cutadapt.\n\n");

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
}

#Check if marker name was provided. If not set it to the name of the first marker in $MARKER_INFO_FILE
unless ($DISABLE_CUTADAPT){
	unless ($MARKER_NAME){
		$MARKER_NAME = $first_marker_name;
		print STDERR "\nWARNING: No marker name was provided. Assuming the marker is the first one listed in $MARKER_INFO_FILE ($MARKER_NAME\). If this is not correct, specify actual marker name with --marker_name.\n\n";
	}
	exists ($markers_HoA{$MARKER_NAME}) or die ("\n$usage\n\nERROR: Could not find the marker name $MARKER_NAME in $MARKER_INFO_FILE\.\n\n");
}

my $total_unfiltered_count = 0;
my $n_count = 0;
my $rep_count = 0;
my $low_barcode_qual_count = 0;
my $edge_count = 0;
my $opt_count = 0;
my $no_trim_count = 0;
my $too_short_count = 0;
my $bbmerge_filter_count = 0;

my $starting_read_count;

#If a read2 file is provided (paired-end reads), run bbmerge.sh to collapse read 1 and read 2.
my $unfiltered_fastq_file;
if ($R2_FASTQ){
	-e $R2_FASTQ or die ("\n$usage\n\nERROR: File specified with --r2_fastq does not exist: $R2_FASTQ.\n\n");
	can_run ($BBMERGE_EXE) or die ("\n$usage\n\nERROR: A read2 file was specified with --r2_fastq, so bbmerge.sh is needed to collapse reads. Could not find following executable in PATH: $BBMERGE_EXE. Provide a path and executable name with --bbmerge_exe.\n\n");
	$unfiltered_fastq_file = "$OUTPUT\.bbmerge.fq";
	
	my $bbmerge_string = "$BBMERGE_EXE in1=$R1_FASTQ in2=$R2_FASTQ out=$unfiltered_fastq_file outu=$OUTPUT\.bbmerge.unmerged_R1.fq outu2=$OUTPUT\.bbmerge.unmerged_R2.fq ihist=$OUTPUT\.bbmerge_hist.txt ordered=t qtrim=r trimq=$BBMERGE_TRIMQ minoverlap=$BBMERGE_MINOVERLAP mismatches=$BBMERGE_MISMATCHES 2> $OUTPUT\.bbmerge_log.txt";
	print "\n" . (localtime) . "\nRead1 and Read2 files both provided. Running $BBMERGE_EXE to collapse reads with the following command:\n\n$bbmerge_string\n\n";
	system ($bbmerge_string);
	
	my @bbmerge_log = file_to_array("$OUTPUT\.bbmerge_log.txt");
	my $bbmerge_join_count;
	foreach (@bbmerge_log){
		if ($_ =~ /^Pairs\:\s+(\d+)\s+$/){
			$starting_read_count = $1;
		}elsif($_ =~ /^Joined\:\s+(\d+)\s+/){
			$bbmerge_join_count = $1;
		}
	}
	$bbmerge_join_count == 0 and die ("\nERROR: $BBMERGE_EXE failed to join any read pairs. Cannot proceed.\n\n");
	($starting_read_count and $bbmerge_join_count) or die ("\nERROR: Could not parse bbmerge.sh log file: $OUTPUT\.bbmerge_log.txt.\n\n");
	$bbmerge_filter_count = $starting_read_count - $bbmerge_join_count;
}else{
	print "\n" . (localtime) . "\nOnly Read 1 file provide (no paired-end). Skipping bbmerge.sh step.\n\n";
}

#if only read1 was provided, use that fastq as input for next step.
$unfiltered_fastq_file or $unfiltered_fastq_file = $R1_FASTQ;
my $FH = open_file ($unfiltered_fastq_file);

#main data structure (hash of arrays) to store array of all reads by barcode. Each read contains 3 concatenated elements (header, seq, quality score) delimited by newlines.
my %barcode_read_HoA;

print "\n" . (localtime) . "\nApplying individual read filters to reads from $unfiltered_fastq_file\.\n\n";

#loop through each read in fastq file
while (my $header = <$FH>){
	++$total_unfiltered_count;
	my $seq = <$FH>;
	my $drop_line = <$FH>;
	my $qual_line = <$FH>;
	
	#extract MDS barcode from the beginning of the read
	my $barcode = substr ($seq, 0, $BARCODE_LEN);
	
	#Exclude if barcode region contains N
	unless ($DISABLE_N_FILTER ){
		$barcode =~ /N/ and ++$n_count and next;
	}

	#Exclude if barcode is all the same character
	unless ($DISABLE_REP_FILTER){
		$barcode =~ /^(.)\1+$/ and ++$rep_count and	next;
	}

	#Exclude if barcode contains any sites with quality score below threshold
	my $low_qual = 0;
	if ($MIN_BAROCODE_QUAL){
		for (my $i=0; $i < $BARCODE_LEN; ++$i){
			ord(substr($qual_line, $i, 1)) - $PHRED_OFFSET < $MIN_BAROCODE_QUAL and $low_qual=1 and ++$low_barcode_qual_count and last;
		}
		$low_qual and next;
	}	

	#Exclude if cluster is too close to low-X or low-Y side of tile.
	if ($MIN_TILE_EDGE_DIST){	
		my @sh1 = split (/\s/, $header);
		my @sh2 = split (/\:/, $sh1[0]);
		unless ($sh2[5] >= $MIN_TILE_EDGE_DIST && $sh2[6] >= $MIN_TILE_EDGE_DIST){
			++$edge_count;
			next;
		}	
	}
	
	#If passed above filters, store read in hash of arrays based on barcode.
	chomp $qual_line;
	push (@{$barcode_read_HoA{$barcode}}, $header . $seq . $qual_line);
}

#set the total count from this file as the starting count if starting count was not already defined (i.e., bbmerge.sh was skipped because only R1 input file).
$starting_read_count or $starting_read_count = $total_unfiltered_count;

if ($R2_FASTQ and $DELETE_INTERMEDIATE_FASTQS){
	unlink("$OUTPUT\.bbmerge.fq");
	unlink("$OUTPUT\.bbmerge.unmerged_R1.fq");
	unlink("$OUTPUT\.bbmerge.unmerged_R2.fq");
}

if ($MIN_OPTICAL_DIST){
	print "\n" . (localtime) . "\nApplying optical distance filter and outputting filtered reads and barcode tags to $OUTPUT\.filt.fq\n\n";
}else{
	print "\n" . (localtime) . "\nOutputting filtered reads and barcode tags to $OUTPUT\.filt.fq\n\n";
}

my $FH_FILT = open_output ("$OUTPUT\.filt.fq");
##loop through the read families defined by each barcode
foreach my $barcode (keys %barcode_read_HoA){
	
	my @read_array = @{$barcode_read_HoA{$barcode}};
	
	#start list of reads to exclude within this read family
	my %exclude_hash;
	
	#if minimum optical distance is specified, loop through all pairs in a family and exclude duplicates that are too close 	
	if($MIN_OPTICAL_DIST){
		for (my $i = 0; $i < scalar (@read_array) - 1; ++$i){
			my @sh1 = split (/\s/, $read_array[$i]);
			my @split_read1 = split (/\:/, $sh1[0]);
		
			for (my $j = $i+1; $j < scalar(@read_array); ++$j){
	
				exists ($exclude_hash{$j}) and next;
			
				my @sh2 = split (/\s/, $read_array[$j]);
				my @split_read2 = split (/\:/, $sh2[0]);

				my $tile1 = $split_read1[2] . "_" . $split_read1[3] . "_" . $split_read1[4];
				my $tile2 = $split_read2[2] . "_" . $split_read2[3] . "_" . $split_read2[4];
			
				#exclude the second read in the pair if they are on the same tile of same lane of same flow cell and distance is less than specified minimum				
				if ($tile1 eq $tile2){
					if (sqrt(abs($split_read1[5] - $split_read2[5])**2 + abs($split_read1[6] - $split_read2[6])**2) < $MIN_OPTICAL_DIST){
						$exclude_hash{$j} = 1;
					}
				}			
			} 
		}	
	}

	#print out reads in fastq format
	for (my $i = 0; $i < scalar (@read_array); ++$i){
		exists ($exclude_hash{$i}) and ++$opt_count and next;
		my @sr = split (/\n/, $read_array[$i]);	
		print $FH_FILT "$sr[0] $barcode\n" . substr($sr[1], $BARCODE_LEN) . "\n+\n" . substr($sr[2], $BARCODE_LEN) . "\n";		
	}
}

#release the barcode-read hash from memory (as much as possible in Perl)
undef %barcode_read_HoA;

#run cutadapt unless this option is turned off
unless ($DISABLE_CUTADAPT){
	#use primer sequences as the adapters to trim. Must reverse complement the reverse primer
	my $adapter1 = $markers_HoA{$MARKER_NAME}[0];
	my $adapter2 = revcom ($markers_HoA{$MARKER_NAME}[1]);
	my $cutadapt_string = "$CUTADAPT_EXE -a $adapter1\...$adapter2\$";
	
	if ($CUTADAPT_REPORT_UNTRIMMED){
		$cutadapt_string .= " --untrimmed-output $OUTPUT\.untrimmed.fq"; 
	}else{
		$cutadapt_string .= " --discard-untrimmed"
	}
	$cutadapt_string .= " --minimum-length $CUTADAPT_MIN_LEN --cores $CUTADAPT_CORES -e $CUTADAPT_ERR_TOL -o $OUTPUT\.trimmed.fq $OUTPUT\.filt.fq > $OUTPUT\.cutadapt_log.txt";

	print "\n" . (localtime) . "\nRunning $CUTADAPT_EXE on $OUTPUT\.filt.fq to trim locus-specific MDS primers. Untrimmed reads are being discarded. Trimmed reads are being written to $OUTPUT\.trimmed.fq, using the following command:\n\n$cutadapt_string\n\n";
	system ($cutadapt_string);

	my @cutadapt_log = file_to_array("$OUTPUT\.cutadapt_log.txt");
	my $cutadapt_start_count;
	my $cutadapt_trim_count;
	foreach (@cutadapt_log){
		if ($_ =~ /^Total\ reads\ processed\:\s+([\d\,]+)\s+$/){
			$cutadapt_start_count = $1;
		}elsif($_ =~ /^Reads\ with\ adapters\:\s+([\d\,]+)\s+/){
			$cutadapt_trim_count = $1;
		}elsif($_ =~ /^Reads\ that\ were\ too\ short\:\s+([\d\,]+)\s+/){
			$too_short_count = $1;
		}
	}
	($cutadapt_start_count and $cutadapt_trim_count) or die ("\nERROR: Could not parse cutadapt log file: $OUTPUT\.cutadapt_log.txt.\n\n");
	$cutadapt_start_count =~ s/\,//g;
	$cutadapt_trim_count =~ s/\,//g;
	$no_trim_count = $cutadapt_start_count - $cutadapt_trim_count;
}

if ($DELETE_INTERMEDIATE_FASTQS){
	$DISABLE_CUTADAPT or unlink ("$OUTPUT\.filt.fq");
}
print "\n" . (localtime) . "\nAnalysis Complete";
$DELETE_INTERMEDIATE_FASTQS and print ". Only the final fastq file has been retained because --delete_intermediate_fastqs was specified.";
print "\n\n";


print "Read counts and filtering:\n";
print "$starting_read_count\tReads in $R1_FASTQ\n";
$R2_FASTQ and print "$bbmerge_filter_count\tFiltered reads that were not merged by bbmerge.sh\n";
$DISABLE_N_FILTER or print "$n_count\tFiltered reads that contained one or more Ns in barcode\n";
$DISABLE_REP_FILTER or print "$rep_count\tFiltered reads that contained the same repetitive base for the their entire barcode\n";
$MIN_BAROCODE_QUAL and print "$low_barcode_qual_count\tFiltered reads that contained one or more base quality scores below $MIN_BAROCODE_QUAL in barcode\n";
$MIN_TILE_EDGE_DIST and print "$edge_count\tFiltered reads that were within $MIN_TILE_EDGE_DIST pixels of left and/or lower border of tile\n";
$MIN_OPTICAL_DIST and print "$opt_count\tFiltered reads that were within $MIN_OPTICAL_DIST pixels of another read with the same barcode\n";
$DISABLE_CUTADAPT or print "$no_trim_count\tFiltered reads for which cutadapt failed to identify a forward and/or reverse primer\n";
$DISABLE_CUTADAPT or ($CUTADAPT_MIN_LEN and print "$too_short_count\tFiltered reads with length less than $CUTADAPT_MIN_LEN after cutadapt trimming\n");