# NetSeekR

A networks analysis pipeline for RNASeq time series data.

NetSeekR is a network analysis R package that includes the capacity to analyze time series of RNASeq data, perform correlation and regulatory network inferences and use network analysis methods to summarize the results of a comparative genomics study. 

Authors: Himangi Srivastava, Drew Ferrell, and George V. Popescu.

The NetSeekR code requires specific versions for packages that are used.  
| Package        | Version     |
| -------------- |-------------| 
| pacman         |0.5.1|
| BiocManager    |1.30.10|
| magrittr       |1.5|
| readr          |1.3.1|
| purrr          |0.4.2|
| stringr        |0.3.3|
| ggplot2        |1.4.0|
| devtools       |3.2.1|
| flashClust     |2.2.1|
| tidyr          |1.01-2|
| networkD3      |1.0.0|
| igraph         |0.4|
| limma          |1.2.4.2|
| edgeR          |3.42.0|
| topGO          |2.37.0|
| WGCNA          |1.68|
| biomaRt        |2.42.0|
| Rgraphviz      |2.30.0|
| dplyr          |0.8.3|

Below are the steps to run NetSeekR.

1. Set the working directory to the NetSeekR path. 

` setwd(<<path/to/NetSeekR>>)`

2. Unzip the NetSeekR file.  

`unzip('NetSeekR.zip')`  

3. Load packages and source functions for NetSeekR.  

`source('scripts/NetSeekR.R')`  

4. Edit configuration file and sample comparison matrix.

 * **note** Below is a template configuration file which needs to be edited per usage.
 
 
 | analysis_type	 | custom tag | 
 | --------------  |-------------| 
 | design_matrix	 | path to experimental design matrix | 
 | edger_adjustment_method	 | edgeR: p-value adjustment method | 
 | edger_lfc	 | limma: minimum log2-fold-change that is considered scientifically meaningful | 
 | edger_method	 | NOT USED | 
 | feature_counts_path	 | path to feature counts software | 
 | kallisto	 | Boolean value for Kallisto execution decision | 
 | kallisto_bias	 | sequence based bias correction | 
 | kallisto_bootstrap_samples	 | bootstrap sample number | 
 | kallisto_chromosomes	 | tab separated file with chromosome names and lengths | 
 | kallisto_fasta_files	 | path to genome annotation file | 
 | kallisto_fastq_files	 | reads to quantify | 
 | kallisto_fr_stranded	 | strand specific reads, first read forward | 
 | kallisto_fragment_length	 | estimated average fragment length | 
 | kallisto_fusion	 | search for fusions for Pizzly | 
 | kallisto_genomebam	 | project pseudoalignments to genome sorted BAM file | 
 | kallisto_gtf	 | GTF file for transcriptome information (required for --genomebam) | 
 | kallisto_index	 | location to write genome index from Kallisto (required for Kallisto alignment) | 
 | kallisto_kmer_size	 | k-mer (odd) length (defaut: 31, max value: 31 | 
 | kallisto_make_unique	 | replace repeated target names with unique names | 
 | kallisto_output_dir	 | directory to write quantification output to | 
 | kallisto_path	 | path to Kallisto software | 
 | kallisto_plaintext	 | output plaintext instead of HDF5 | 
 | kallisto_pseudobam	 | save pseudoalignments to transcriptome to BAM file | 
 | kallisto_rf_stranded	 | strand specific reads, first read reverse | 
 | kallisto_sd	 | estimated standard deviation of fragment length (default: -l, -s values are estimated from paired end data, but are required when using --single) | 
 | kallisto_seed	 | seed for the bootstrap sampling (default: 42) | 
 | kallisto_single	 | quantify single-end reads | 
 | kallisto_single_overhang	 | include reads where unobserved rest of fragment is predicted to lie outside a transcript | 
 | kallisto_threads	 | number of threads to use (default: 1) | 
 | sample_comparisons_file	 | path to the sample comparison file for differential gene expression testing | 
 | sample_covariates	 | experimental design matrix column names to be used as covariates with Sleuth | 
 | significance_cutoff	 | a cutoff value for determining significance | 
 | sleuth_gene_mode	 | Boolean value for Sleuth gene mode execution decision | 
 | sleuth_transcript_mode	 | Boolean value for Sleuth transcript mode execution decision | 
 | star	 | Boolean value for STAR execution decision | 
 | star_genomeDir	 | path to the directory where the genome indices are stored | 
 | star_genomeFastaFiles	 | path to a FASTA file with the genome reference sequences | 
 | star_path	 | path to STAR software | 
 | star_readFilesIn	 | path to the folder containing the sequences to be mapped | 
 | star_runThreadN	 | number of threads to be used for genome generation, it has to be set to the number of available cores on the server node | 
 | star_sjdbGTFfile	 | path to the file with annotated transcripts in the standard GTF format | 
 | star_sjdbOverhang	 | length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database> | 



5. Align batches of reads.

`alignment_results <- implement_alignment(arguments_file = <<path/to/configuration file>>)`

6. Test for differential gene expression.

 * **note** An edited sample comparison matrix needs to be supplied for differential gene expression testing. An example of an edited sample comparison matrix file is below. A sample comparison matrix does not need headers, only sample identifiers. Which samples to compare in differential testing should be written row-wise. 
 
|Example | | | | | |
| --------------  |-------------| --------------  |-------------| --------------  |-------------|
|SL209924    |SL209925    |SL209926    |SL209921    |SL209922    |SL209923|
|SL209927    |SL209928    |SL209937    |SL209921    |SL209922    |SL209923|
|SL209938    |SL209939    |SL209940    |SL209921    |SL209922    |SL209923|
|SL209944    |SL209945    |SL209946    |SL209941    |SL209942    |SL209943|
|SL209947    |SL209948    |SL209949    |SL209941    |SL209942    |SL209943|
|SL209950    |SL209951    |SL209952    |SL209941    |SL209942    |SL209943|

 
`implement_differential_gene_expression(alignment_results)`
 
 7. Network analysis is then conducted assuming sets of differentially expressed genes are available.
 
 `implement_network_analysis(alignment_tool = 'star', alignment_results = alignment_results)`
 
 The network analysis function above will produce an image of the network. The image can be saved from the ‘Export’ tab at the top of the window.
