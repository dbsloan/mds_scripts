# mds_scripts
Scripts for analyzing MDS (maximum depths sequencing) datasets

These scripts are used to analyze Illumina sequencing datasets generated with methods described by [Jee et al. (2016 Nature)](https://www.nature.com/articles/nature18313). The initial input files are raw Illumina fastq reads (gzipped is OK). Note that unlike, the original MDS protocol, these scripts assume that libraries have been demultiplexed using a standard i5/i7 demultiplexing strategy and that the reads in the R1 file start with a random barcode (unique molecular identified) followed by a locus-specific fwd primer, the region of interest and potentially the locus-specific reverse primer and some of the i7 Illumina adapter.


The following three scripts can be called in order.

- `MDS_process_reads.pl`: Merges and trims fastq reads to region of interest
- `MDS_process_families.pl`: Summarizes read families from trimmed fastq
- `MDS_process_variants.pl`: Produces variant calls from read family summaries



For processing multiple libraries in batches, the following wrapper script can be called to run the above three scripts in series for each library.

- `MDS_process_pipeline.pl`
 


## MDS_process_pipeline.pl

Usage: `perl MDS_process_pipeline.pl [options] --libraries_file=INPUT_FILE`
   
   This is a wrapper script that calls a series for three scripts to process
   MDS datasets.
   
   REQUIRED ARGUMENTS
   
   File summarizing input libraries. See [example file](sample_config_files/input_libraries_file.txt).

   
         --libraries_file
         A tab-delimited text file that provides output name, R1 fastq file,
         marker info file, and optionally: R2 fastq file, marker name, filter
         seq file. See sample template (make sure new line characters can be 
         parsed by Perl in a UNIX environment).




   OPTIONAL ARGUMENTS (no defaults)

   File providing non-default options for each of the three processing scripts. See [example file](sample_config_files/options.txt)

         --options_file      
         A tab-delimited text file. First column is "reads", "families", 
         or "variants" to indicate the which script the option applies to.
         The second column is the option. See sample template.


   Scripts paths/names 

         --reads_script      
         Full path to MDS_process_reads.pl if not in working directory

         --families_script      
         Full path to MDS_process_families.pl if not in working directory

         --variants_script      
         Full path to MDS_process_variants.pl if not in working directory

   
   EXAMPLE
   
         perl MDS_process_pipeline.pl \
         	--libraries_file=INPUT_FILE --options_file=OPTIONS_FILE \
           	--reads_script=/PATH/TO/SCRIPTS/MDS_process_reads.pl \
           	--families_script=/PATH/TO/SCRIPTS/MDS_process_families.pl \
           	--variants_script=/PATH/TO/SCRIPTS/MDS_process_variants.pl
       	


## MDS_process_reads.pl

Usage: `perl MDS_process_reads.pl [options/arguments]`
   
   This script takes raw MDS Illumina fastq input files and returns a primer-
   trimmed fastq with the MDS barcode added to the header of each read 
   (OUTPUT_NAME.trimmed.fq). A number of filters are applied and filtering 
   metrics are reported in additional output files and STDOUT.
   
   REQUIRED ARGUMENTS
   
   R1 Fastq File
   
         --r1_fastq
         File containing Illumina read 1 MDS sequences (raw sequences). Can 
         be gzipped.

   MDS Marker Information File.  See [example file](sample_config_files/MDS_markers.txt).
   
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
   
         perl MDS_process_reads.pl \
         	--r1_fastq=READ1.fastq.gz -r2_fastq=READ2.fastq.gz \
         	--marker_info_file=MDS_markers.txt --marker_name=At_mt1 \
         	--min_optical_dist=5000 --min_tile_edge_dist=2000 \
         	--delete_intermediate_fastqs \
         	--output=my_output



## MDS_process_families.pl


Usage: `perl MDS_process_families.pl [options/arguments]`
   
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
   
         perl MDS_process_families.pl \
         	--exclude_n_reads \
         	--input=Trimmed_Input.fq --output=my_output


## MDS_process_variants.pl

Usage: `perl MDS_process_variants.pl [options/arguments]`
   
   This script takes the read fams text file from MDS_process_families.pl and
   summarizes variant calls.
   
   REQUIRED ARGUMENTS
   
   Input File
   
         --input
         Read fams output text file from MDS_process_families.pl

   MDS Marker Information File. See [example file](sample_config_files/MDS_markers.txt)
   
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
         Note that this option does not affect the tabulation of counts in the "single"
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
   
         perl MDS_process_variants.pl \
         	--marker_info_file=MDS_markers.txt --marker_name=At_mt1 \
         	--min_base_agreement=1 \
         	--input=INPUT_READ_FAMS_FILE.txt --output=my_output
       	

