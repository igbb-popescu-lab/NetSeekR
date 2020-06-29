# Load all packages and source all functions.
if (!'pacman' %in% installed.packages()){
  install.packages('pacman')
}

library(pacman)

p_load(kable, BiocManager, magrittr, readr, rlang, purrr, stringr, ggplot2, devtools,
       flashClust, tidyr, networkD3, igraph, scales, reshape, tibble, ggraph, tidygraph)

# Install packages from Bioconductor.
bioconductor_packages <- c('limma', 'edgeR', 'topGO', 'sleuth','WGCNA', 'biomaRt', 'Rgraphviz', 'STRINGdb', 'BiocManager')
    
for (package in bioconductor_packages){
  if (!package %in% installed.packages()){
    BiocManager::install(package, update = F)
  }
}

p_load(char = bioconductor_packages, install = F)

rm(bioconductor_packages)

# Prevent conflict of 'select' function in AnnotationDBi.
p_load(dplyr)

# Source all necessary scripts.
list.files('./scripts/', pattern = 'implement.*R$', full.names = TRUE) %>% 
  walk(source)