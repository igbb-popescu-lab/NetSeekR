# NetSeekR
A networks analysis pipeline for RNASeq time series data

NetSeekR is a network analysis R package that includes the capacity to analyze time series of RNASeq data, perform correlation and regulatory network inferences and use network analysis methods to summarize the results of a comparative genomics study. 

Authors: Himangi Srivastava, Drew Ferrell, and George V. Popescu.

1. Set the working directory to the NetSeekR path.  
setwd(<<path/to/NetSeekR>>)  
2. Unzip the NetSeekR file.  
unzip('NetSeekR.zip')  
3. Load packages and source functions for NetSeekR.  
source('scripts/NetSeekR.R')  
4. Edit configuration file and sample comparison matrix as per NetSeekR.html
5. Execute pipeline with the three main functions: implement_alignment.R, implement_differential_expression.R, and implement_network_analysis.R, as per NetSeekR.html
