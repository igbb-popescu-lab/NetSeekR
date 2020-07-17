# NetSeekR
A networks analysis pipeline for RNASeq time series data

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


1. Set the working directory to the NetSeekR path.  
setwd(<<path/to/NetSeekR>>)  
2. Unzip the NetSeekR file.  
unzip('NetSeekR.zip')  
3. Load packages and source functions for NetSeekR.  
source('scripts/NetSeekR.R')  
4. Edit configuration file and sample comparison matrix as per NetSeekR.html
5. Execute pipeline with the three main functions: implement_alignment.R, implement_differential_expression.R, and implement_network_analysis.R, as per NetSeekR.html
