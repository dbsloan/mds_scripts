# mds_scripts
Scripts for analyzing MDS (maximum depths sequencing) datasets

These scripts are used to analyze Illumina sequencing datasets generated with methods described by [Jee et al. (2016 Nature)](https://www.nature.com/articles/nature18313). The initial input files are raw Illumina fastq reads (gzipped is OK). Note that unlike, the original MDS protocol, these scripts assume that libraries have been demultiplexed using a standard i5/i7 demultiplexing strategy and that the reads in the R1 file start with a random barcode (unique molecular identified) followed by a locus-specific fwd primer, the region of interest and potentially the locus-specific reverse primer and some of the i7 Illumina adapter.


\


The following three scripts can be called in order.

- `MDS_process_reads.pl`: Merges and trims fastq reads to region of interest
- `MDS_process_families.pl`: Summarizes read families from trimmed fastq
- `MDS_process_variants.pl`: Produces variant calls from read family summaries

\

For processing multiple libraries in batches, the following wrapper script can be called to run the above three scripts in series for each library.

- `MDS_process_pipeline.pl`: Merges and trims fastq reads to region of interest
 


## MDS_process_pipeline.pl

Usage: perl MDS_process_pipeline.pl [options] --libraries_file=INPUT_FILE
   
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
   
         perl MDS_process_pipeline.pl \\
         	--libraries_file=INPUT_FILE --options_file=OPTIONS_FILE \\
           	--reads_script=/PATH/TO/SCRIPTS/MDS_process_reads.pl \\
           	--families_script=/PATH/TO/SCRIPTS/MDS_process_families.pl \\
           	--variants_script=/PATH/TO/SCRIPTS/MDS_process_variants.pl
       	

