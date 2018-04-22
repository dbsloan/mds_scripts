#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use sloan;

my $usage = 
"\nUsage: perl $0 [options/arguments]
   
   This script takes variant output files from one or more libraries analyzed with 
   MDS_process_variants.pl and summarizes various statistical measures, including 
   variant frequencies by read-family size, mutation spectrum, and position specific
   variant counts.
   
   REQUIRED ARGUMENTS
   
   Input Libraries. Input files must be specified in one of two ways. The output 
   basename (but not full ouptut filenames) from MDS_process_families.pl

         --input_lib_list
         Inputs can be provided with in a simple text file (one library per line) with
         this option.
         
         Or they can be provided as a space-delimited list at the end of the command
         line call.   
 
    Family Size Threshold for Variant Frequencies
   
         --fam_size_freqs
         A cutoff for read family size. All families of this size and larger will be 
         pooled for reporting variant frequencies.
        
    Family Size Threshold for Variant Sprectrum and Position-Specific Information
   
         --fam_size_spect
         A cutoff for read family size. All families below this size will be discarded
         for variant-spectrum and position-distribution summaries.

   Output Name
   
         --output [default: input file name]
         Base name for all output files (additional extensions will be added)


   OPTIONAL ARGUMENTS

   Include multi-SNP families
   
         --include_multisnps
         Add this flag to include read families that have multiple variants relative to 
         the reference. By default, only families with a single SNP are used. If this 
         option is specified, mutliSNP families should have been manually inspected to
         confirm that they are genuine.



   EXAMPLES
   
         perl $0 \\
         	--input_lib_list=INPUT_LIB_LIST_FILE \\
         	--fam_size_freqs=5 --fam_size_spect=5 \\
         	--output=my_output 

         perl $0 \\
         	--fam_size_freqs=5 --fam_size_spect=5 \\
         	--output=my_output \\
         	Lib1 Lib2 Lib3 ...

       	
\n\n";


our $INPUT_LIB_LIST;
our $FAM_SIZE_FREQS;
our $FAM_SIZE_SPECT;
our $OUTPUT;
our $INCLUDE_MULTISNPS;

GetOptions(
    'input_lib_list=s'  => \$INPUT_LIB_LIST,
	'fam_size_freqs=i'  => \$FAM_SIZE_FREQS,
	'fam_size_spect=i'  => \$FAM_SIZE_SPECT,
    'output=s'  => \$OUTPUT,
    'include_multisnps'  => \$INCLUDE_MULTISNPS
);



#Check for proper inputs and return usage statement if there are errors
$OUTPUT or die ("\n$usage\n\nERROR: Must provide a basename for output files with --output.\n\n");
$INPUT_LIB_LIST or @ARGV or die ("\n$usage\n\nERROR: Must provide one or more input library names either with --input_lib_list or as a space-delimited list at the end of the command line.\n\n");
(int($FAM_SIZE_FREQS) == $FAM_SIZE_FREQS and $FAM_SIZE_FREQS >= 0) or die ("\n$usage\n\nERROR: --fam_size_freqs is required and must be a non-negative integer.\n\n");
(int($FAM_SIZE_SPECT) == $FAM_SIZE_SPECT and $FAM_SIZE_SPECT >= 0) or die ("\n$usage\n\nERROR: --fam_size_spect is required and must be a non-negative integer.\n\n");

my @libs;

if ($INPUT_LIB_LIST){
	@ARGV and die ("\n$usage\n\nERROR: Input libraries provided with --input_lib_list and at end of command line. Only one of these methods sould be used to input library names.\n\n");
	@libs = file_to_array ($INPUT_LIB_LIST);
	for (my $i =0; $i < scalar (@libs); ++$i){
		chomp $libs[$i];
	}
}else{
	while (my $lib = shift){
		push (@libs, $lib);
	}
}

#single mismatches are stored in col index 3, unless multisnp variants are included (col index 2)
my $col_num = 3;
$INCLUDE_MULTISNPS and $col_num = 2;

#store variant counts by family size
my %matches_HoH;
my %mismatches_HoH;
foreach (@libs){
	my @var_summ_array = file_to_array ("$_\.var_summ.txt");
	my $first_line = shift @var_summ_array;
	foreach my $line (@var_summ_array){
		chomp $_;
		my @sl = split (/\t/, $line);
		if ($sl[0] >= $FAM_SIZE_FREQS){
			$matches_HoH{$_}->{$FAM_SIZE_FREQS} += $sl[1];
			$mismatches_HoH{$_}->{$FAM_SIZE_FREQS} += $sl[$col_num];
		}else{
			$matches_HoH{$_}->{$sl[0]} = $sl[1];
			$mismatches_HoH{$_}->{$sl[0]} = $sl[$col_num];
		}
	}
}


#store spectrum and position information
my %spect_HoH;
my %pos_HoAoA;
my $max_len = 0;
my %nuc_pos_hash;
foreach (@libs){
	my $file = "$_\.var.single.txt";
	$INCLUDE_MULTISNPS and $file = 	"$_\.var.txt";
	my @var_array = file_to_array ($file);
	my $first_line = shift @var_array;
	foreach my $line (@var_array){
		chomp $line;
		my @sl = split (/\t/, $line);
		$sl[4] < $FAM_SIZE_SPECT and next;
		if (exists ($nuc_pos_hash{$sl[0]})){
			$nuc_pos_hash{$sl[0]} eq $sl[2] or die ("\nERROR: Different reference nucleotides detected at position $sl[0]\. All libraries must use the same marker.\n\n");
		}else{
			$nuc_pos_hash{$sl[0]} = $sl[2];
		}
		$spect_HoH{$_}->{"$sl[2]$sl[3]"} += $sl[5];
		if ($sl[1] eq "Match"){
			$pos_HoAoA{$_}[$sl[0]][0] += $sl[5];
		}elsif ($sl[1] eq "Mismatch"){
			$pos_HoAoA{$_}[$sl[0]][1] += $sl[5];			
		}else{
			die ("\nERROR: Could not parse match/mismatch status from following line in $file:\n$line\n\n");
		}
		$sl[0] > $max_len and $max_len = $sl[0];
	}
}

#print freqs output
my $FHF = open_output ("$OUTPUT\.var_freqs.txt");

print $FHF "FamilySize";
foreach (@libs){
	print $FHF "\t$_\_Matches";
}
foreach (@libs){
	print $FHF "\t$_\_Mismatches";
}
foreach (@libs){
	print $FHF "\t$_\_Freqs";
}
print $FHF "\n";

for (my $i=1; $i <= $FAM_SIZE_FREQS; ++$i){
	if ($i == $FAM_SIZE_FREQS){
		print $FHF "$i\+";
	}else{
		print $FHF $i;		
	}
	
	foreach (@libs){
		print $FHF "\t$matches_HoH{$_}->{$i}";
	}
	foreach (@libs){
		print $FHF "\t$mismatches_HoH{$_}->{$i}";
	}
	foreach (@libs){
		my $freq = $mismatches_HoH{$_}->{$i} / ($mismatches_HoH{$_}->{$i} + $matches_HoH{$_}->{$i});
		print $FHF "\t$freq";
	}
	print $FHF "\n";
}

close $FHF;


#print spect output
my $FHS = open_output ("$OUTPUT\.var_spect.txt");
my @nucs = ("A", "C", "G", "T");

print $FHS "SubType";
foreach (@libs){
	print $FHS "\t$_\_Count";
}
print $FHS "\n";

foreach my $nuc1 (@nucs){
	foreach my $nuc2 (@nucs){
		my $dinuc = $nuc1 . $nuc2;
		print $FHS $dinuc;
		foreach my $lib (@libs){
			my $count = 0;
			exists ($spect_HoH{$lib}->{"$dinuc"}) and $count = $spect_HoH{$lib}->{"$dinuc"};
			print $FHS "\t$count";
		}
		print $FHS "\n";
	}
}

close $FHS;

#print pos output
my $FHP = open_output ("$OUTPUT\.var_pos.txt");

print $FHP "Position\tReference";
foreach (@libs){
	print $FHP "\t$_\_Matches";
}
foreach (@libs){
	print $FHP "\t$_\_Mismatches";
}
foreach (@libs){
	print $FHP "\t$_\_Freqs";
}
print $FHP "\n";

for (my $i=1; $i <= $max_len; ++$i){
	print $FHP "$i\t$nuc_pos_hash{$i}";
	
	foreach (@libs){
		exists ($pos_HoAoA{$_}[$i][0]) or $pos_HoAoA{$_}[$i][0] = 0;
		print $FHP "\t$pos_HoAoA{$_}[$i][0]";
	}
	foreach (@libs){
		exists ($pos_HoAoA{$_}[$i][1]) or $pos_HoAoA{$_}[$i][1] = 0;
		print $FHP "\t$pos_HoAoA{$_}[$i][1]";
	}
	foreach (@libs){
		my $freq = $pos_HoAoA{$_}[$i][1] / ($pos_HoAoA{$_}[$i][0] + $pos_HoAoA{$_}[$i][0]);
		print $FHP "\t$freq";
	}
	print $FHP "\n";
}

close $FHP;