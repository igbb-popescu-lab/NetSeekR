# Function name: implement_network_analysis
# Purpose: Performs network visualization and find hubbed genes 
#          in the network using Go enrichment analysis also.
# Input: A string specifying which differential gene expression 
#        was used upstream.
# Output: Network visualization and GO terms associated with 
#         significant genes. 
implement_network_analysis <- function(alignment_tool, alignment_results, exectute){
  alignment_decision <- alignment_results %>% 
    dplyr::first()
  pipeline_input <- alignment_results %>% 
    dplyr::last() 
  
  # WGCNA
  if(alignment_tool %>% str_detect('kallisto|KALLISTO|Kallisto')){
    wgcna_input <- alignment_decision %>% 
      filter(key == 'kallisto') %>% 
      select(directories) %>% 
      unlist(use.names = FALSE) %>% 
      first() %>% 
      str_replace('Kallisto', paste0('Sleuth/', pipeline_input %>% extract2('analysis_type'), '/gene_level/gene_level_results'))
  }
  
  if(alignment_tool %>% str_detect('STAR|star|Star')){
    wgcna_input <- pipeline_input %>% 
      extract2('star_genomeDir') %>% 
      str_replace('STAR.*', paste0('edgeR/', pipeline_input$analysis_type, '/')) %>% 
      str_replace(pattern = '//', replacement = '/')
    
    
    # Remove tags after gene names in edgeR results.
    deg_files <- wgcna_input %>% 
      list.files(path = ., pattern = '.tsv$', full.names = TRUE, recursive = FALSE) %>% 
      tibble(deg_file = .)
    
    test_for_extraneous <- deg_files %>% 
      use_series(deg_file) %>% 
      first() %>% 
      read_tsv(col_types = cols()) %>%
      use_series(target_id) %>% 
      first() %>% 
      unlist(use.names = FALSE)
    
    # str matching specific to test data.
    if(test_for_extraneous %>% str_detect('\\..*')){
      deg_files_remove <- deg_files %>% 
        mutate(
          deg_data = map(deg_file,
                         ~read_tsv(.x, col_types = cols()) %>%
                           mutate(
                             target_id = map_chr(target_id, 
                                                 ~gsub(pattern = '\\..*', replacement = '', x = .x)
                             )
                           )
          ),
          save_deg_data = map2(deg_data, 
                               deg_file,
                               write_tsv
          )
        )
    }
  }
  
  message('Running WGCNA.')
  
  expr2 <- wgcna_input_data(wgcna_input, pipeline_input) 
  
  significant_hits <- expr2 %>% 
    select(original_filter) %>% 
    unnest(cols = original_filter) %>% 
    distinct()
  
  expr2 <- expr2 %>% 
    wgcna_plot_sample_tree() %>% 
    wgcna_plot_power_results() %>% 
    wgcna_plot_power_histogram() %>% 
    wgcna_clustering()
  
  
  message('Performing GO enrichment.')
  GO_results <- implement_GO_enrichment(deg_tool = wgcna_input, alignment_results = alignment_results) %>% 
    post_process_GO_results()
  
  expr2 <- expr2 %>%
    mutate(
      GO_results = GO_results
    )
  
  # DREM
  message('Writing executable files for DREM.')
  time_series_count_data <- DREM_main(pipeline_input = pipeline_input, wgcna_input = wgcna_input, significant_hits, execute)
  # Overlap differentially expressed genes with network
  DREM_network_overlap(pipeline_input, deg_files)
  
  message('Conducting network analysis.')
  network_analysis_results <- mapping_network_analysis(expr2)
  
  adjacency_matrices <- expr2 %>% 
    mutate(
      comparison_adj_mat = map(clustering, ~.x$adjacency)
    ) %>% 
    select(comparison, comparison_adj_mat) %>% 
    mutate(
      gene_tf = map(comparison_adj_mat, map_gene_tf_network_analysis, pipeline_input)
    )
  
  
}






map_gene_tf_network_analysis <- function(adj_mat, pipeline_input){
  
  setwd(paste0(getwd(), '/data/DREM/', pipeline_input$analysis_type))
  
  paths <- getwd() %>% 
    list.files(pattern = 'path', full.names = T) %>% 
    tibble(path_file = .) %>% 
    dplyr::mutate(
      path_data = map(path_file, vroom, col_types = cols()),
      path_data = map(path_data, ~select(.x, !contains('SPOT') & !starts_with('H'))),
      path_data = map(path_data, ~set_colnames(.x, gsub("\\|.*","", colnames(.x))))
    ) 
  
  tf_list <- '../../../Ath_TF_list.txt'
  
  paths <- paths %>% 
    mutate(
      tmp = purrr::map(path_data, 
                       gene_tf_network_analysis, 
                       adj_mat, 
                       tf_list
      )
    )
  
  
}

  
  


gene_tf_network_analysis <- function(drem_gene_tf, adj_mat, tf_list){
  
  #drem_gene_tf=read.table("gene_tf_table.txt",header = TRUE)
 
  
  
  
  ids = drem_gene_tf %>% 
    select(1) 
  
  drem_gene_tf <- as.data.frame(drem_gene_tf)
  
  
  rownames(drem_gene_tf)=drem_gene_tf[,1]
  drem_gene_tf<- drem_gene_tf[,-1]
  rr=drem_gene_tf[rowSums(drem_gene_tf)>2,]
  rr=rr[,colSums(rr)>2]
  
  
  g <- graph.incidence(rr)
  V(g)$color <- V(g)$type
  V(g)$color=gsub("FALSE","red",V(g)$color)
  V(g)$color=gsub("TRUE","blue",V(g)$color)
  tkp.id<-tkplot(g, edge.color="gray30", layout=layout_as_bipartite)
  
  tk_center(tkp.id)
  #tk_fit(tkp.id, width = 467, height = 567)
  #tk_rotate(tkp.id, degree = -90, rad = NULL)
  
  
  types_drem <- V(g)$type                 ## getting each vertex `type` let's us sort easily
  degree_drem <- igraph::degree(g)
  betweenness_drem<- betweenness(g)
  closeness_drem <- closeness(g)
  eig_drem <- eigen_centrality(g)$vector
  
  centrality_drem <- data.frame(types_drem , degree_drem, betweenness_drem, closeness_drem, eig_drem)
  
  centrality_drem[order(centrality_drem$types_drem, decreasing = TRUE),]
  
  #adj_m <- adj_mat[(rownames(adj_mat)%in% ids),]
  
  adjacancy_matrix_of_drem_TF <- adj_mat[(rownames(adj_mat) %in% ids),]
  
  #adjacency mtrix for drem TF data that matches with genes
  adjacancy_matrix_of_drem_TF = adj_mat[,colnames(adj_mat)[colnames(adj_mat) %in% colnames(drem_gene_tf)]]
  
  
  
  ath_tf_list=read.table(tf_list, header = TRUE)
  ath_tf_list=gsub("\\..*","",ath_tf_list$TF_ID)
  ath_tf_list=as.data.frame(unique(ath_tf_list))
  #colnames(ath_tf_list)="TFs"
  #adjacancy matrix of database TF data that matches with genes
  adjacancy_matrix_of_database_TF=adj_mat[,colnames(adj_mat)[colnames(adj_mat) %in% ath_tf_list]]
  
  
  data_tf_gene_corr <- unique(cbind(adjacancy_matrix_of_database_TF,adjacancy_matrix_of_drem_TF))
  
  
  m <- data_tf_gene_corr
  source_node=c()
  target_node=c()
  correlation<-c()
  genes_names <- rownames(data_tf_gene_corr)
  tf_names<-colnames(data_tf_gene_corr)
  i<-genes_names[1:dim(m)[1]]
  j<-tf_names[1:dim(m)[2]]

  for(gene in i)
  {
    for(gen in j)
    {
      if(m[gene,gen]>0.5){
        source_node<-c(source_node,gene)
        target_node<-c(target_node,gen)
        correlation<-c(correlation,m[gene,gen])
      }
    }
  }
  
  
  
  NetworkData <- data.frame(source_node, target_node, correlation)
  
  
  net=NetworkData %>% filter(correlation > 0.5)
  net=net %>% filter(correlation != 1)
  
  net=unique(net)
  
  g <- graph.empty(directed = F)
  node.out <- unique(net$target_node) #stringsAsFactor = F in data frame
  node.in <- unique(net$source_node) #stringsAsFactor = F in data frame
  g <- graph.data.frame(net, directed = F)
  V(g)$type <- V(g)$name %in% net[,2] #the second column of edges is TRUE type
  E(g)$weight <- as.numeric(net[,3])
  g
  
  rr=get.incidence(g,attr = "weight")
  
  
  
  V(g)$color <- V(g)$type
  V(g)$color=gsub("FALSE","red",V(g)$color)
  V(g)$shape=gsub("FALSE","square",V(g)$shape)
  V(g)$color=gsub("TRUE","blue",V(g)$color)
  tkplot(g, edge.color="gray30",edge.width=E(g)$weight, layout=layout_as_bipartite)
  
  
  types <- V(g)$type    
  deg <- igraph::degree(g)
  bet <- betweenness(g)
  clos <- closeness(g)
  eig <- eigen_centrality(g)$vector
  
  cent_df <- data.frame(types, deg, bet, clos, eig)
  
  cent_df[order(cent_df$type, decreasing = TRUE),]
}
  
  #quantile foreach gene -TF relationship
  # q <- vector()
  # for (i in 1:length(data_tf_gene_corr)){
  #   q[i]<- quantile(data_tf_gene_corr[,i],  probs = 0.9)
  #   
  # }
  # 
  
  
  
  
  # adj_mat=read.table("gene-gene_adj.txt",header = TRUE)
  # row.names(adj_mat)=adj_mat$Column1
  # adj_mat=adj_mat[,-1]
  # 
  # 
  # row.names(adj_mat)=adj_mat$Column1
  # adj_mat=adj_mat[,-1]
  # 
  #adj<-as.data.frame(adjacency)
  #adj<-adj_mat
  
  #since DREM filters some data so we will keep the genes that are not filtered by DREM
  #drem_gene_tf=read.table("ILK_path1.txt",header = TRUE)
  #drem_gene_tf = select_if(!contains('SPOT', 'H[[:digit:]]+'))
    
    
    #drem_gene_tf[,-(2:6)]
  
  #colnames(paths$path_data[[1]]) = gsub("\\|.*","", colnames(paths$path_data[[1]]))
  
  






DREM_network_overlap <- function(p, d){
  # p: pipeline input
  # d: DEG file tibble
  
  
  
  # point to network directory
  nets <- p %>% 
    extract2('DREM') %>% 
    paste0('TFInput') 
  # list potential networks to access
  
  net_full <- nets %>% 
    list.files(full.names = T)
  # files without source target filtered away.
  net_full_tib <- net_full %>% 
    tibble(f = .) %>% 
    mutate(
      x = map_lgl(f, ~readLines(.x, n = 1) %>% 
                    str_detect('TF\tGene\tInput')
      )
    ) %>% 
    filter(x)
  
  nets <- net_full_tib %>% 
    select(f) %>% 
    unlist(use.names = F)
  
  # invite user to select network
  message('Please pick an item (corresponding to a network) to overlap differentially expressed genes with: ')
  # format menu to select from
  net_select <- nets %>% 
    as_tibble(.) %>% 
    rownames_to_column() %>% 
    dplyr::rename('item' = rowname, 'network' = value)
  # print selection menu
  print(net_select)
  # input selection item
  chosen_item <- readline('Type item number: ')
  
  # extract network file name
  net <- net_select %>% 
    spread(item, network) %>% 
    select(all_of(chosen_item))
  
  # extract full path to network chosen
  path_to_chosen_net <- net_full[str_detect(net_full, net %>% unlist(use.names = F))]
  # check which tool
  edgeR_tool <- d %>% 
    unlist(use.names = F) %>% 
    first() %>% 
    str_detect(pattern = 'edgeR')
  if (edgeR_tool){
    pval_filter <- 'adj.P.Val'
    summarise_by <- 'logFC'
    file_name_pattern <- 'clean_counts'
    
    # Select the directory which contains the count data from STAR.
    counts_data <- d %>% 
      dplyr::slice(1) %>% 
      first() %>% 
      str_replace(pattern = 'edgeR.*', replacement = 'STAR/Feature_counts/')
    
    sk <- 1
    
  }
  if (!edgeR_tool){
    message('Unable to analyze Sleuth results.')
    #pval_filter <- 'qval'
    #summarise_by <- 'b'
    
  }
  # column id for removing extraneous chars.
  gene_col <- 'target_id'
  # source and target definitions as in network files.
  target_source <- 'TF'
  target <- 'Gene'
  clean_dm <-  p %>% 
    extract2('design_matrix') %>% 
    clean_design_matrix(., p)
  # load deg files
  summarised_differentially_expressed_gene_sets <- d %>%
    mutate(
      
      # Read in differentially expressed gene sets; filter based on p-val cut-off
      de_gene_sets = map(deg_file, ~read_tsv(.x, col_types = cols()) %>% 
                           drop_na() %>% 
                           filter(!!sym(pval_filter) <= p %>% extract2('significance_cutoff'))
                        
      ),
      summarised_de_gene_sets = map(de_gene_sets, 
                                    ~summarise_differential_expression(differential_expression_gene_set = .x, gene_column = gene_col, mean_summerize = summarise_by)
      ),
      extract_title = map_chr(deg_file, extract_title_from_file_name)
    ) %>% 
    select(-de_gene_sets)
  # unzip if needed
  # if (str_detect(string = path_to_chosen_net, '.gz$')){
  #   unzip(zipfile = path_to_chosen_net)
  #   path_to_chosen_net <- path_to_chosen_net %>% str_remove('.gz')
  # }
  
  network <- read.delim(path_to_chosen_net) %>% 
    select(all_of(target_source), all_of(target))
  # Make network edges unique.
  edge_list <- network %>% 
    mutate(
      !!target_source := !!sym(target_source) %>% map(extract_edges)
    ) %>%
    unnest(cols = c(target_source)) %>% 
    mutate(
      !!target := !!sym(target) %>% map_chr(extract_edges)
    ) 
  
  # Extract unique source nodes. 
  unique_sources <- edge_list %>% 
    unique_and_relabel('label_name' = target_source)
  
  # Extract unique target nodes. 
  unique_targets <- edge_list %>% 
    unique_and_relabel(label_name = target)
  
  # Combine all unique sources and all unique targets and label each a unique identifier.                                  
  all_nodes <- full_join(unique_sources, unique_targets, by = 'label') %>% 
    rowid_to_column('id')
  
  # Associate unique gene ids with edges. 
  edge_list_ids <- edge_list %>%
    left_join(all_nodes, by = setNames(nm = target_source, 'label')) %>%
    dplyr::rename(from = id) %>% 
    left_join(all_nodes, by = setNames(nm = target, 'label')) %>%
    dplyr::rename(to = id)
  
  
  # Number of columns the extracted title can be divided into.
  column_number <- extract_condition_column_number(extracted_title_dataset = summarised_differentially_expressed_gene_sets)
  
  # Group comparison set elements.
  comparison_elements <- associate_comparison_elements(extracted_title_dataset = summarised_differentially_expressed_gene_sets, number_of_columns = column_number)
  
  # Associate expression data (counts) replicate sets with conditions.
  associate_expression_replicates <- associate_replicate_sets_to_conditions(expression_data_location = counts_data, clean_dm, sk, p)
  
  # Associate replicate sets to differential gene expression comparison sets. 
  associate_expression_values_to_comparison_sets <- associate_expression_to_comparison_elements(replicate_set_association = associate_expression_replicates, comparisons = comparison_elements)
  
  # Uniquely name the column in the large tibble to contain overlapped nodes. 
  overlap_sources_and_targets_column <- paste(target_source, target, sep = '_and_')
  
  
  # Associate read count data (expression data) with comparison sets (comparison elements). 
  comparison_set_expression <- associate_expression_values_to_comparison_sets %>% 
    mutate(
      
      # Paste groups together to join with extract_title column in node_overlap variable.
      extract_title = map2_chr(group1, group2, paste_groups_to_join_with_title_extract),
      
      # Take the mean expression value for replicates, and the mean expression value across each comparison.
      gene_expression = map2(group1_replicate_data, group2_replicate_data, mean_replicate_est_counts)
    )
  
  
  # Overlap differentially expressed genes with network edges.
  node_overlap <- summarised_differentially_expressed_gene_sets %>%
    
    # remove '.tsv'
    
    mutate(
      extract_title = str_remove(extract_title, pattern = '\\..*') %>% 
        str_remove('_vs')
    ) %>% 
    
    dplyr::right_join(comparison_set_expression, by = 'extract_title') %>% 
    
    select(-contains('group')) %>% 
    
    mutate(
      
      
      # Extract read count data to associate along side differential expression measurement. 
      summarised_de_gene_sets = map2(summarised_de_gene_sets, gene_expression, ~inner_join(.x, .y, by = gene_col)),
      
      
      # Collect target sources found in each differentially expressed gene set. 
      !!overlap_sources_and_targets_column := map(summarised_de_gene_sets, ~overlap_sources(differential_gene_expression_set = .x, edges = edge_list_ids, node_type = target_source, summarise_col = summarise_by)),
      
      
      # Collect target measurements similarly.
      !!overlap_sources_and_targets_column := map2(summarised_de_gene_sets, !!sym(overlap_sources_and_targets_column), ~overlap_targets(differential_gene_expression_set = .x, target_sources = .y, node_type = target, summarise_col = summarise_by)),
      
      
      # Adjacency matrix for each network. 
      adjacency_matrix = map(!!sym(overlap_sources_and_targets_column), get_adjacency_matrix),
      
      
      # Collect edge information.
      edges = map(!!sym(overlap_sources_and_targets_column), ~select(.x, matches('from|to'))),
      
      # ID the nodes within each set.
      relabel_nodes = map(!!sym(overlap_sources_and_targets_column), ~relabel_local_nodes(local_edges = .x, target_sources = target_source, targets = target)),
      
      
      # Map beta values to unique genes.
      map_betas_to_nodes = map(!!sym(overlap_sources_and_targets_column), ~beta_and_read_mapping_to_nodes(target_source_set = .x, node_type_1 = target_source, node_type_2 = target, summarise_col = summarise_by)),
      
      
      # Label each identifier as either target or source.
      label_type = map(map_betas_to_nodes, ~label_and_rescale_mapped_values(mapped_values_set = .x, edges = edge_list, first_label = target_source, second_label = target))
      
      
      # Extract in-degree and out-degree for each node in a network.
      #node_degree = map2(label_type, adjacency_matrix, ~map_degrees_to_nodes(node_information = .x, adj_matrix = .y, node_universe = all_nodes))
      
    )
  
  tmp <- getwd() %>% paste0('/data/', p$analysis_type, '_DEG_networks')
  dir.create(tmp, recursive = T)
  
  expression_networks <- node_overlap %>%
    mutate(
      #
      graph_data = pmap(list(map_betas_to_nodes, edges, label_type, extract_title), graph_differential_gene_expression_network),
      
      #
      write_graph_data = map2(extract_title, graph_data, ~paste_and_save(output_directory = tmp, file_name_to_paste = .x, graph = .y))
    )
}







graph_differential_gene_expression_network <- function(nodes, edges, label_types, titles_extracted) {
  
  as_tbl_graph(nodes, edges %>% unlist()) %>%
    activate(nodes) %>%
    left_join(label_types, by = 'name') %>%
    dplyr::rename('Identifier' = name) %>%
    dplyr::rename('Class' = label) %>% 
    ggraph(layout = 'kk') +
    geom_edge_link(arrow = arrow(length = unit(3, 'mm'))) +
    geom_node_point(aes(alpha = scaled_beta, colour = Identifier, size = reads, shape = Beta_coefficient)) +
    geom_node_text(aes(label = Class), color = 'black', size = 2, vjust = 2, show.legend = FALSE) +
    theme_graph() +
    ggtitle(titles_extracted) + 
    labs(alpha = 'Scaled LFC value', size = 'Mean of counts', shape = 'LFC') +
    guides(color = FALSE)
}



paste_and_save <- function(output_directory, file_name_to_paste, graph){
  output_directory <- output_directory %>% 
    paste0('/', file_name_to_paste, '.png') 
  
  graph %>% 
    ggsave(filename = output_directory, device = 'png', dpi = 320, width = 10.00, height = 10.00, units = 'in')
  
}



label_and_rescale_mapped_values <- function(mapped_values_set, edges, first_label, second_label){
  
  # Extract value range for scaling. 
  rescale_set <- mapped_values_set %>%
    select(beta) %>%
    range()
  
  
  mapped_values_set <- mapped_values_set %>%
    mutate(
      
      # Rescale values for graphing with alpha. 
      scaled_beta = map_dbl(beta, ~rescale(x = .x, from = rescale_set, to = c(0,1))),
      
      # Label beta coefficient.
      Beta_coefficient = map_chr(beta, ~if_else(condition = (.x < 0), 
                                                true = '-', 
                                                false = '+'))
      
    )
  
  mapped_values_set %>%
    left_join(x = ., y = edges, by = c('name' = first_label)) %>%
    dplyr::rename(!!first_label := second_label) %>%
    
    
    mutate_at(first_label, ~if_else(condition = is.na(.),
                                    true = replace(x = ., values = second_label),
                                    false = replace(x = ., values = first_label)
    )
    ) %>%
    dplyr::rename(label := !!first_label) %>%
    distinct()
}








beta_and_read_mapping_to_nodes <- function(target_source_set, node_type_1, node_type_2, summarise_col){
  node_type_1_set <- c(node_type_1, paste0(node_type_1, '_', summarise_col), paste0(node_type_1, '_mean_counts'))
  node_type_2_set <- c(node_type_2, paste0(node_type_2, '_', summarise_col), paste0(node_type_2, '_mean_counts'))
  
  target_source_set %>% 
    nest(data = c(node_type_1_set, node_type_2_set)) %>% 
    dplyr::rename(tmp := data) %>% 
    dplyr::select(tmp) %>% 
    mutate(
      tmp = map(tmp, gather) 
    ) %>% 
    unnest(cols = c(tmp)) %>% 
    select(value) %>% 
    mutate(
      ind = rep(c(1,2,3), length.out = n())
    ) %>% 
    group_by(ind) %>% 
    mutate(
      id = row_number()
    ) %>%
    spread(ind, value) %>%
    dplyr::select(-id) %>%
    dplyr::rename(name = '1', beta = '2', reads = '3') %>%
    mutate_at(vars(beta, reads), as.numeric) %>% 
    distinct()
}



relabel_local_nodes <- function(local_edges, target_sources, targets){
  local_edges %>% 
    select(target_sources, targets) %>% 
    unlist(use.names = FALSE) %>% 
    tibble(name = .) %>% 
    distinct() %>% 
    rowid_to_column('id')
}


get_adjacency_matrix <- function(overlap_sources_and_targets_column){
  overlap_sources_and_targets_column %>%
    select(to, from) %>%
    as.data.frame() %>%
    graph.data.frame() %>%
    get.adjacency() %>%
    as.matrix()
}



overlap_targets <- function(differential_gene_expression_set, target_sources, node_type, summarise_col){
  differential_gene_expression_set %>% 
    dplyr::rename(!!node_type := target_id) %>% 
    right_join(target_sources, by = node_type) %>% 
    drop_na() %>% 
    dplyr::rename(!!paste0(node_type, '_', summarise_col) := paste0(summarise_col, '_mean')) %>% 
    dplyr::rename(!!paste0(node_type, '_mean_counts') := mean_counts)
}




overlap_sources <- function(differential_gene_expression_set, edges, node_type, summarise_col){
  differential_gene_expression_set %>% 
    dplyr::rename(!!node_type := target_id) %>%
    left_join(edges, by = node_type) %>% 
    drop_na() %>% 
    dplyr::rename(!!paste0(node_type, '_', summarise_col) := paste0(summarise_col, '_mean')) %>% 
    dplyr::rename(!!paste0(node_type, '_mean_counts') := mean_counts)
}




average_counts_across_comparison_sets <- function(est_count_set){
  est_count_set %>% 
    dplyr::rename(target_id = 1, counts = 2, counts1 = 4) %>% 
    select(ends_with('id'), contains('counts')) %>% 
    mutate(
      mean_counts = map2_dbl(counts, counts1, ~mean(x = c(.x, .y)))
    ) %>% 
    select(target_id, mean_counts)
}


average_est_counts_for_replicates <- function(group){
  group %>% 
    unnest(cols = c(count_data)) %>%
    group_by(target_id) %>% 
    summarise(
      counts = mean(counts, na.rm = TRUE)
    )
}

mean_replicate_est_counts <- function(column_one, column_two){
  
  column_one <- column_one %>% average_est_counts_for_replicates()
  column_two <- column_two %>% average_est_counts_for_replicates()
  bind_cols(column_one, column_two) %>% 
    average_counts_across_comparison_sets(est_count_set = .)
}



paste_groups_to_join_with_title_extract <- function(column_one, column_two){
  column_one %>% 
    paste(column_two, sep = ' ') %>%
    str_replace_all(pattern = ' ', replacement = '_')
}


associate_expression_to_comparison_elements <- function(replicate_set_association, comparisons){
  comparisons %>% 
    inner_join(replicate_set_association, by = c('group1' = 'condition')) %>%
    dplyr::rename(group1_replicate_data = replicate_data) %>% 
    
    left_join(replicate_set_association, by = c('group2' = 'condition')) %>%
    dplyr::rename(group2_replicate_data = replicate_data)
}




associate_comparison_elements <- function(extracted_title_dataset, number_of_columns){
  
  start <- (number_of_columns - number_of_columns) + 1
  mid <- number_of_columns/2
  mid_right <- round(mid) + 1
  
  extracted_title_dataset %>% 
    select(extract_title) %>% 
    tidyr::extract(col = extract_title, 
                   into = rep('id', times = number_of_columns) %>% paste(1:number_of_columns, sep = ''), 
                   regex = rep('(.*)', times = number_of_columns) %>% paste(collapse = '_')) %>% 
    
    unite(col = group1, rep('id', times = (mid)) %>% paste0(start:mid), sep = ' ') %>% 
    unite(col = group2, rep('id', times = (mid)) %>% paste0(mid_right:number_of_columns), sep = ' ') %>% 
    mutate(
      group2 = str_remove(group2, '.tsv')
    )
}


clean_design_matrix <- function(dm, p){
  #p: pipeline input
  covars <- p$sample_covariates
  dm %>% 
    read_csv(col_types = cols()) %>% 
    tidyr::extract(condition, 
            into = str_split(pattern = ', ', covars) %>% unlist(),
            regex = '(.*) (.*) ([[:digit:]].*)')
}


associate_replicate_sets_to_conditions <- function(expression_data_location, clean_dm, sk, p){  
  tmp <- str_split(pattern = ', ', p$sample_covariates) %>% unlist()
  expression_data <- expression_data_location %>% 
    tibble(files = list.files(path = ., full.names = TRUE, recursive = TRUE)) %>% 
    
    mutate(
      
      files = map_chr(files, ~str_replace(.x, pattern = '//', replacement = '/')),
      count_data = map(files, 
                       read_tsv_filter_extraneous,
                       sk
      ),
      sample = map_chr(files, ~basename(.x) %>% str_remove('\\_.*'))
    ) %>% 
    
    
    # Join by sample identifiers.
    left_join(y = clean_dm, by = 'sample') %>% 
    dplyr::select(files, count_data, sample, all_of(tmp)) %>%
    tidyr::unite(col = condition, tmp, sep = ' ') %>% 
    select(count_data, condition) %>% 
    group_by(condition) %>% 
    nest() %>% 
    dplyr::rename(replicate_data = data)
}



read_tsv_filter_extraneous <- function(current, s){
  # current: current counts file
  # s: skip lines
  current %>% 
    read_tsv(file = ., col_names = T, skip = s, col_types = cols()) %>% 
    dplyr::rename(target_id = 1, est_counts = 2) %>% 
    mutate(
      target_id = str_remove(string = target_id, '\\..*')
    ) %>% 
    group_by(target_id) %>% 
    summarise(counts = mean(est_counts)) %>% 
    ungroup()
}



associate_expression_to_comparison_elements <- function(replicate_set_association, comparisons){
  comparisons %>% 
    inner_join(replicate_set_association, by = c('group1' = 'condition')) %>%
    dplyr::rename(group1_replicate_data = replicate_data) %>% 
    
    left_join(replicate_set_association, by = c('group2' = 'condition')) %>%
    dplyr::rename(group2_replicate_data = replicate_data)
}



extract_condition_column_number <- function(extracted_title_dataset){
  extracted_title_dataset %>% 
    select(extract_title) %>% 
    mutate(
      column_number = map_dbl(extract_title, ~str_count(.x, '_') %>% add(1) %>% as.numeric())
    ) %>% 
    select(column_number) %>% 
    unlist(use.names = FALSE) %>% 
    unique()
}



unique_and_relabel <- function(edges, label_name){ 
  edges %>%
    distinct(!!sym(label_name)) %>% 
    dplyr::rename(label = !!label_name)
}

extract_edges <- function(field){
  field %>% 
    str_split(pattern = '\\|') %>% 
    unlist() %>% 
    unique()
}

# remove extraneous chars in gene column. summarise chosen column by mean of values in column.
summarise_differential_expression <- function(differential_expression_gene_set, gene_column, mean_summerize){
  differential_expression_gene_set %>% 
    mutate(
      !!gene_column := map_chr(!!sym(gene_column), ~str_remove(.x, pattern = '\\..*'))
    ) %>% 
    group_by(!!sym(gene_column)) %>% 
    summarize(!!sym(paste0(mean_summerize, '_mean')) := mean(!!sym(mean_summerize), na.rm = T))
}

extract_title_from_file_name <- function(file_name){
  file_name %>% 
    basename() %>% 
    str_remove('_results[[:punct:]].*')
}

# Function name: wgcna_input_data
# Purpose: Takes the input data and convert it into 
#          differentially expreesed expression data for wgcna analysis.
# Input: path to differentially expressed genes and pipeline input. 
# Output: filtered count data.
wgcna_input_data <- function(differentialy_expressed_genes, pipeline_input){
  
  # Load design matrix into environment.
  design_matrix <- pipeline_input$design_matrix %>% 
    read_csv(col_types = cols())
  
  # Detect whether analyzing differentially expressed genes from
  # edgeR or Sleuth.
  edgeR_tool <- differentialy_expressed_genes %>% 
    str_detect(pattern = 'edgeR')
  
  # Select the column which contains the p-values or q-values
  # depending on which differential expression tool was selected. 
  if(edgeR_tool){
    
    # The edgeR analysis output has 'adj.P.Val' as the p-value
    # column name.
    filter_value <- 'adj.P.Val'
    
    # Counts data from feature counts contains the 'clean_counts'
    # string in the nane. 
    file_name_pattern <- 'clean_counts'
    
    # Select the directory which contains the count data from STAR.
    counts_data <- differentialy_expressed_genes %>%
      str_replace('edgeR.*', 'STAR/Feature_counts/')
    
  } 
  else if(!edgeR_tool){
    
    # Sleuth has q-values in the 'qval' column.
    filter_value <- 'qval'
    
    # Counts data are in TSV files from Kallisto.
    file_name_pattern <- '^abundance.tsv$'
    
    counts_data <- differentialy_expressed_genes %>% 
      str_replace('Sleuth.*', 'Kallisto/Kallisto_quantifications/')
  }
  
  # Load read count data and filter lowly expressed genes. 
  keep_counts <- counts_data %>%
    
    # Read count data into environment; use explicit patterns
    # which allow distinctions between counts datasets.
    counts_keep(file_name_pattern, recurse_subdirectories = !edgeR_tool) %>% 
    
    
    # Map count data to sample names from design matrix.
    #associate_counts_with_sample_names(design_matrix) %>%
    
    # Filter lowly expressed genes from count data.
    filter_counts_data_wgcna(counts = ., design_matrix = design_matrix) %>% 
    
    mutate(
      target_id = map(target_id, 
                      ~gsub(pattern = '\\..*', replacement = '', x = .x)
      ),
      
      target_id = target_id %>% 
        as.character()
    )
  
  
  
  
  # Load differential gene expression data.
  differentially_expressed_genes_keep(differentialy_expressed_genes, filter_value, pipeline_input) %>%
    mutate(
      original_filter = deg_data2,
      deg_data2 = map(deg_data2,
                      ~suppressMessages(semi_join(keep_counts, .x))
      ),
      preprocess_wgcna_input = map(deg_data2, 
                                   wgcna_data_processing
      )
    )
  

}




# Function name: counts_keep
# Purpose: Load the count data into the environment using decisions
#          from which count data type is being sourced.
# Input: Counts data, file name regex, Boolean value.
# Output: Counts data from files.
counts_keep <- function(feature_counts, file_name_pattern, recurse_subdirectories){
  
  # Skip lines for feature counts data, but not for Kallisto data. 
  if(recurse_subdirectories){
    nrow_skip <- 0
  } else{
    nrow_skip <-  1
  }
  
  # Load count data into environment based on the file name pattern.  
  list.files(feature_counts, pattern = file_name_pattern, full.names = TRUE, recursive = recurse_subdirectories) %>% 
    tibble(files = .) %>%
    
    # Allow column names to be passed from data source and skip a number of rows. 
    mutate(
      sample_count_data = map(files, read_tsv, skip = nrow_skip, col_names = TRUE, col_types = cols())
    )
}




# Function name: filter_counts_data_wgcna
# Purpose: Filter the lowly expressed genes from the data 
#          obtained from feature counts.
# Input: Counts data and the experimental design matrix. 
# Output: Counts mapped to sample in a dataframe. 
filter_counts_data_wgcna <- function(counts, design_matrix){
  
  # Check for Kallisto data.
  is_kallisto <- counts %>% 
    slice(1) %>% 
    select(files) %>% 
    str_detect('Kallisto_quantifications')
  
  # Cut extra statistics, leaving only the count data
  # for genes. 
  if(is_kallisto){
    gene_names <- counts %>% 
      slice(1) %>% 
      select(sample_count_data) %>% 
      unnest(sample_count_data) %>% 
      select(gene)
    
    counts <- counts %>% 
      mutate(
        sample_count_data = map(sample_count_data, reset_splice_variants, gene_names = gene_names)
      )
  }
  
  samplename <- design_matrix
  
  data <- data.frame(counts$sample_count_data)
  
  # Storing data as data frame.
  total_samples <- nrow(counts)*2
  ans <- seq(2,total_samples,2)
  data = data[c(1,ans)]
  
  #(optional) naming the samples
  #d=colnames(sleuth_table)[1]
  d = "target_id"
  x <- samplename$sample
  a <- c(d,x)
  colnames(data) <- a
  rowname <- data[c(1)]
  
  #filterng the data
  keep <- rowSums(cpm(data[,2:length(data)]) > 0.5) >= 2
  
  data[keep,] 
}




# Function name: reset_splice_variants
# Purpose: Remove splice variant notation.
# Input: Count data and unique gene names.
# Output: Gene names without spliced notation.
reset_splice_variants <- function(sample_counts, gene_names){
  sample_counts %>% 
    mutate(
      target_id = gene_names %>% 
        unlist(use.names = F)
    ) %>% 
    select(target_id, est_counts)
}





# Function name: wgcna_data_processing
# Purpose: Preprocess the input data for wgcna.
# Input: Count data specific to a comparison set from 
#        differential gene expression testing.  
# Output: Transposed count data. 
wgcna_data_processing <- function(comparison_set_counts){
  
  Expression_data0 <- comparison_set_counts[, -c(1)]
  Expression_data0 <- as.data.frame(t(Expression_data0))
  names(Expression_data0) = comparison_set_counts$target_id
  goodsamples = goodSamplesGenes(Expression_data0, verbose = 0)
  if (!goodsamples$allOK)
  {
    if (sum(!goodsamples$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(Expression_data0)[!goodsamples$goodGenes], collapse = ", ")));
    if (sum(!goodsamples$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(Expression_data0)[!goodsamples$goodSamples], collapse = ", ")));
    Expression_data0 = Expression_data0[goodsamples$goodSamples, goodsamples$goodGenes]
  }
  Expression_data0 <- as.data.frame(Expression_data0)
  return(Expression_data0)
}




# Function name: differentially_expressed_genes_keep
# Purpose: Find the count data of differentially expressed 
#          genes from TSV files of EdgeR/ Sleuth results 
#          regardless of which tool was used.
# Input: DGE data directory, p-value, and pipeline input structure.
# Output: Genes filtered with the cut-off (p-value).
differentially_expressed_genes_keep <- function(deg_directory, filter_value, pipeline_input){
  # Count data files are TSV files regardless of which tool was used.
  tibble(files = list.files(deg_directory, pattern = '.tsv', full.names = TRUE)) %>%
    mutate(
      comparison = map_chr(files, 
                           extract_comparison
      ),
      deg_data2 = map(files, 
                      read_and_filter, 
                      filter_value, 
                      pipeline_input
      )
    )
}





# Function name: extract_comparison
# Purpose: Find the comparision file name and remove .tsv string pattern
# Input: Filename for a differential gene expression test result.
# Output: Input without .tsv extension.
extract_comparison <- function(comparison_filename){
  comparison_filename %>% 
    basename() %>% 
    str_remove('.tsv')
}



# Function name: read_and_filter
# Purpose: Reads the data from Sleuth/ EdgeR and 
#          filters out the differentially expressed genes.
# Input: Set of differentially expressed genes, cut-off value,
#        and pipeline input structure.
# Output: The input set filtered on the cut-off value.
read_and_filter <- function(differentialy_expressed_genes_sets, filter_value, pipeline_input){
  differentialy_expressed_genes_sets %>% 
    read_tsv(col_types = cols()) %>% 
    drop_na() %>% 
    filter(!!sym(filter_value) <= pipeline_input$significance_cutoff) %>% 
    select(target_id)
}






# Function name: wgcna_plot_sample_tree
# Purpose: Plot to find any outlier sample in the data 
#          for all the comparisions. 
# Input: Expression_data0.
# Output: Plotted sample tree.
wgcna_plot_sample_tree <- function(Expression_data0){
  
  Expression_data0 %>%
    mutate(
      sample_tree = map(preprocess_wgcna_input,
                        hclust_distance_matrix,
                        'average'
      ),
      
      # Plot the sample tree: 
      sample_clustering = pmap(
        list(sample_tree,
             preprocess_wgcna_input,
             files,
             comparison),
        plot_sample_tree
      )
    )
}




# Function name: hclust_distance_matrix
# Purpose: Find the distance matrix to perform clustering.
# Input: Preprocessed WGCNA data and the 'average' method. 
# Output: Analyzed hierarchical cluster.
hclust_distance_matrix <- function(Expression_data0, method){
  Expression_data0 %>% 
    dist() %>% 
    hclust(method)
}




# Function name: plot_sample_tree
# Purpose: Plot to find any outlier sample in the data .
# Input: hclust result, preprocessed WGCNA data, file name from 
#        differential gene expression test, and the names of
#        comparisons made. 
# Output: List: datExp, sft, power, and k values. 
plot_sample_tree <- function(sample_tree, Expression_data0, file_name, comparison){
  sizeGrWindow(12,9)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  
  p.plot <- sample_tree %>% plot(main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  
  choose_line <- readline(prompt = "Enter the number at which the limit that shoud be cut or remove the outlier Hint:the sample that seems to be the outlier: ")
  choose_line <- as.integer(choose_line)
  
  if (is.na(choose_line)){
    print("Enter a valid number.")
  }
  
  
  sample_clust <- sample_tree %>% 
    cutreeStatic(cutHeight = choose_line, minSize = 10)
  
  keepSamples <- (sample_clust==1)
  
  datExpr <- Expression_data0[keepSamples, ]
  powers <-  c(c(1:10), seq(from = 12, to=100, by=2))
  invisible(capture.output(sft <-  pickSoftThreshold(datExpr, powerVector = powers, verbose = 0, networkType = "signed")))
  power <- sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq==max(sft$fitIndices$SFT.R.sq))]
  k <- softConnectivity(datE = datExpr, power = power, verbose = 0)
  
  return(list(datExpr = datExpr, sft = sft, power = power, k = k))
}



# Function name: wgcna_plot_power_results
# Purpose: Histogram plot to analyse and choose 
#          correct power value for all the comparisions 
# Input: Column which contains the sft values. 
# Output: Plot in Rstudio Plots pane.
wgcna_plot_power_results <- function(contains_sft){
  contains_sft %>% 
    mutate(
      power_results = map2(sample_clustering, 
                           comparison,
                           plot_power_results)
    )
}




# Function name: plot_power_results
# Purpose: Histogram plot to analyse and choose correct power value.
# Input: datExp, sft, power, and k values; comparison names. 
# Output: Analysis for scale-free topology and mean connectivity.
plot_power_results <- function(sample_clustering, comparison){
  comparison <- comparison %>% 
    str_replace_all('_', ' ') %>% 
    paste('Analyzing', .)
  
  print(comparison)
  
  sft <- sample_clustering %>% 
    use_series(sft)
  
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  powers = c(c(1:10), seq(from = 12, to=100, by=2))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  
  abline(h=0.90,col="red")
  
  p1.plot<-plot(sft$fitIndices[,1], sft$fitIndices[,5],
                xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
                main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  readline(prompt="Press [enter] to see next figure. ")
}





# Function name: wgcna_plot_power_histogram
# Purpose: Plot to choose the correct power 
#          from scale free topology for all the comparisions.
# Input: Column with k values. 
# Output: Plot in Rstudio Plots pane.
wgcna_plot_power_histogram <- function(contains_k){
  contains_k %>% 
    mutate(
      power_histogram_results = map2(sample_clustering,
                                     comparison,
                                     plot_power_histogram)
    )
}




# Function name: plot_power_histogram
# Purpose: Plot to choose the correct power 
#          from sclae free topology.
# Input: datExp, sft, power, and k values; comparison names.
# Output: Scale free topology plot.
plot_power_histogram <- function(sample_clustering, comparison){
  
  comparison <- comparison %>% 
    str_replace_all('_', ' ') %>% 
    paste('Analyzing', .)
  
  print(comparison)
  
  k <- sample_clustering %>% 
    use_series(k)
  
  
  sizeGrWindow(10,5)
  par(mfrow=c(1,2))
  p3.plot<-hist(k)
  p4.plot<-scaleFreePlot(k, main="Check scale free topology\n")
  
  readline(prompt = "Press [enter] to see next figure. ")
  
  return(list(p3.plot, p4.plot))
  
}



# Function name: wgcna_clustering
# Purpose: Perform clustering and saves result in their respective directory.
# Input: datExp, sft, power, and k values.
# Output: Plot in Rstudio Plots pane.
wgcna_clustering <- function(sample_clustering){
  sample_clustering %>% 
    mutate(
      clustering = map2(sample_clustering,
                        files, 
                        clustering)
    )
  
}





# Function name: wgcna_clustering
# Purpose: Perform clustering and saves result in their respective directory.
# Input: datExp, sft, power, and k values; file path to write dendrograms to.
# Output: datExp, dynamic_Modules, adjacency, and dynamic_Colors values.
clustering <- function(sample_clustering, file_path){
  
  cluster_dend_dir <- file_path %>% 
    dirname() %>% 
    paste0('/cluster_dendrograms/') %>% 
    str_replace(pattern = '//', replacement = '/')
  
  dir.create(cluster_dend_dir, showWarnings = FALSE, recursive = TRUE)
  
  file_path <- file_path %>% 
    basename() %>% 
    paste0(cluster_dend_dir, .) %>% 
    str_replace('.tsv$', '.pdf')
  
  datExpr <- sample_clustering %>% 
    use_series(datExp)
  
  sft <- sample_clustering %>% 
    use_series(sft)
  
  softPower <- sample_clustering %>% 
    use_series(softPower)
  
  k <- sample_clustering %>% 
    use_series(k)
  
  
  adjacency = adjacency(datExpr, power = 12)
  
  # topological overlap matrix
  TOM = TOMsimilarity(adjacency, verbose = 0)
  dissTOM = 1-TOM
  
  # hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = 30
  dynamic_Modules = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                  deepSplit = 2, pamRespectsDendro = FALSE,
                                  minClusterSize = minModuleSize, verbose = 0)
  #convert colors 
  
  dynamic_Colors = labels2colors(dynamic_Modules)
  
  
  Module_eigenList = moduleEigengenes(datExpr, colors = dynamic_Colors)
  Module_eigenvals = Module_eigenList$eigengenes
  Module_eigenDiss = 1-cor(Module_eigenvals)
  
  
  module_eigenTree = flashClust(as.dist(Module_eigenDiss), method = "average"); 
  Module_eigenDissThres = 0.25
  
  merge = mergeCloseModules(datExpr, dynamic_Colors, cutHeight = Module_eigenDissThres, verbose = 0)
  
  
  mergedColors = merge$colors;
  
  sizeGrWindow(12, 9)
  pdf(file = file_path, wi = 9, he = 6)
  p.plot <- plotDendroAndColors(geneTree, cbind(dynamic_Colors, mergedColors),
                                c("Dynamic Tree Cut", "Merged dynamic"),
                                dendroLabels = FALSE, hang = 0.03,
                                addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  
  mergedMEs = merge$newMEs
  
  return(list(datExpr = datExpr, 
              dynamic_Modules = dynamic_Modules,
              adjacency = adjacency,
              dynamic_Colors = dynamic_Colors))
  
}



# Function name: implement_GO_enrichment
# Purpose: Checks the DEG data on which GO analysis to be performed 
#          and performs GO analysis.
# Input: String specifying which tool was used; alignment results
#        structure.
# Output: GO enrichment results.
implement_GO_enrichment <- function(deg_tool, alignment_results){
  
  alignment_decision <- alignment_results %>% 
    dplyr::first()
  
  pipeline_input <- alignment_results %>% 
    dplyr::last()
  
  if(deg_tool %>% str_detect('Sleuth')){
    deg_input <- alignment_decision %>% 
      filter(key == 'kallisto') %>% 
      select(directories) %>% 
      unlist(use.names = FALSE) %>% 
      first() %>% 
      str_replace('Kallisto', paste0('Sleuth/', pipeline_input$analysis_type, '/gene_level/gene_level_results'))
    
    pattern <- '.tsv$'
    
    filter_value <- 'qval'
    
    GO_enrichment <- enrichment(deg_directory = deg_input, pattern = pattern)
  }
  
  if(deg_tool %>% str_detect('edgeR')){
    deg_input <- pipeline_input$star_genomeDir %>% 
      str_replace('STAR.*', paste0('edgeR/', pipeline_input$analysis_type, '/')) %>% 
      str_replace(pattern = '//', replacement = '/')
    
    pattern <- '.tsv$'
    
    filter_value <- 'adj.P.Val'
    
    GO_enrichment <- enrichment(deg_directory = deg_input, pattern = pattern, pipeline_input, filter_value)
    
  }
  
}




# Function name: enrichment
# Purpose: Test for significantly enriched genes in the deg gene sets.
# Input: Directory to deg files, file name pattern, pipeline input, and
#        filter value. 
# Output: Go enrichment. 
enrichment <- function(deg_directory, pattern, pipeline_input, filter_value){
  
  
  deg_file_tib <- list.files(path = deg_directory, pattern = pattern, full.names = TRUE, recursive = TRUE) %>% 
    tibble(deg_file = .) %>% 
    mutate(
      deg_data = map(deg_file, 
                     read_tsv,
                     col_types = cols()
      )
    )
  
  
  mart_arabdopsis <- biomaRt::useMart(biomart = "plants_mart",
                                      dataset = "athaliana_eg_gene",
                                      host = 'plants.ensembl.org', 
                                      verbose = FALSE)
  
  
  Gene_go <- suppressMessages(biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id"), mart = mart_arabdopsis, verbose = F))
  
  geneID2GO <- by(Gene_go$go_id,
                  Gene_go$ensembl_gene_id,
                  function(x) as.character(x))
  
  filename_for_go <- extract_filename_go(deg_file_tib)
  go_input <- deg_file_tib$deg_data
  
  tmp <- list()
  gene_name <- list()
  result_GO_gene <- list()
  for (i in 1:length(go_input))
  {     
    tmp[[i]] <- go_input[[i]] %>% 
      filter(!!sym(filter_value) <= pipeline_input$significance_cutoff)
    
    gene_name[[i]] <- tmp[[i]] %>% 
      dplyr::select("target_id")
    
    if (filter_value == 'pval') {
      tmp[[i]] <- tmp[[i]] %>%
        dplyr::select( 'pval' ) 
    } else {
      tmp[[i]] <- tmp[[i]] %>%
        dplyr::select( 'adj.P.Val' ) 
    }
    
    
    geneList <- as.numeric(unlist(tmp[[i]]))
    names(geneList) <- as.character(unlist(gene_name[[i]]))
    
    # Create topGOData object
    GOdata <- suppressMessages(
      new("topGOdata",
          ontology = "BP",
          allGenes = geneList,
          geneSelectionFun = function(x)(x == 1),
          annot = annFUN.gene2GO, gene2GO = geneID2GO)
    )
    
    # Kolmogorov-Smirnov testing
    result_KS <- suppressMessages(
      runTest(GOdata, algorithm = "weight01", statistic = "ks")
    )
    
    GO_result_tab <- GenTable(GOdata, raw.p.value = result_KS, topNodes = length(result_KS@score), numChar = 120)
    
    par(cex = 1)
    
    showSigOfNodes(GOdata, score(result_KS), firstSigNodes = 10, useInfo = "def")
    
    print(head(GO_result_tab))
    
    printGraph(GOdata, result_KS, firstSigNodes = 10, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
    
    my_GO_term <- c(GO_result_tab$GO.ID[1])
    my_genes <- genesInTerm(GOdata, my_GO_term)
    
    
    for (j in 1:length(my_GO_term))
    {
      my_GO_term2 <- my_GO_term[j]
      my_genes_GO_term <- my_genes[my_GO_term2][[1]]
      my_genes_GO_term <- paste(my_genes_GO_term, collapse=',')
      print(paste("Term",my_GO_term,"genes:", my_genes_GO_term))
      
    }
    result_GO_gene[i] <- my_genes_GO_term
  }
  return(list(result_GO_gene = result_GO_gene, filename_for_go = filename_for_go))
}




# Function name: extract_filename_go
# Purpose: Extracts the comparision set for which 
#          GO enrichmen and network analysisi to be performed.
# Input: Tibble with differentially expressed genes. 
# Output: Edited file name. 
extract_filename_go <- function(deg_file_tib){
  filename <- deg_file_tib$deg_file
  filename<-filename %>% 
    str_remove('.tsv')
  filename<-sub(".*/", "", filename)
}





# Function name: post_process_GO_results
# Purpose: Assosiate file names with GO results.
# Input: GO results. 
# Output: Named list.
post_process_GO_results <- function(GO_results){
  file_names <- GO_results[length(GO_results)] %>% 
    unlist(use.names = FALSE)
  GO_results[length(GO_results)] <- NULL
  names(GO_results$result_GO_gene) <- file_names
  GO_results <- GO_results %>% 
    flatten()
}




# Function name: DREM_main
# Purpose: Extract required data for executing DREM.
# Input: Pipeline input structure, wgcna_input.
# Output: None
DREM_main <- function(pipeline_input, wgcna_input, sig, exec){
  
  analysis_type <- pipeline_input %>%
    extract2('analysis_type')
  
  if(wgcna_input %>% str_detect('Sleuth')){
    map_dir <- pipeline_input %>% 
      extract2('kallisto_output_dir') %>% 
      paste0('Kallisto_quantifications/')
    
    dataset <- 'kallisto'
  }
  
  if(wgcna_input %>% str_detect('edgeR')){
    map_dir <- pipeline_input %>% 
      extract2('star_genomeDir') %>% 
      str_replace('genome_Dir', 'Feature_counts/')
    
    dataset <- 'Feature_counts'
  }
  
  # Set paths to DREM input data to be created. 
  design_matrix <- pipeline_input %>% 
    extract2('design_matrix') 
  
  drem_time_series_input_path <- getwd() %>%
    paste0('/data/DREM/', analysis_type, '/')
  
  
  dir.create(drem_time_series_input_path, recursive = TRUE, showWarnings = FALSE)
  
  drem_defaults_template <- pipeline_input %>% 
    extract2('DREM') %>% 
    paste0('/defaults.txt') %>% 
    str_replace('//', '/')
  
  
  design_matrix <- design_matrix %>% 
    read_csv(col_types = cols()) %>%
    dplyr::select(sample, condition, which_replicate)
  
  
  loaded_read_data <- load_read_data(map_dir, design_matrix, dataset)
  
  tmp <- rearrange_count_data(loaded_read_data, pipeline_input$sample_covariates) %>% 
    # reduce by sig hits
    mutate(
      reads = map(data, right_join, sig, by = 'target_id')
    )
  
  time_series_count_data <- tmp %>% 
    write_DREM_time_series_data(., pipeline_input, drem_time_series_input_path) %>% 
    write_default_files(defaults_template = drem_defaults_template) %>% 
    ungroup() 
  
  
  DREM <- pipeline_input %>% 
    extract2('DREM') %>%
    list.files(pattern = '*.jar', full.names = T)
  
  # Execute DREM script.
  DREM_command_line <- paste('java -mx1024M -jar', DREM, '-b')
  
  # Write batch DREM to script. 
  DREM_execute <- time_series_count_data %>% 
    select(new_default_file_name) %>% 
    mutate(
      command = map_chr(new_default_file_name,
                        ~paste(DREM_command_line, .x)
      ),
      file = str_replace(command, '^.* (.*)$', '\\1') %>% 
        str_replace('.txt$', '_outfile.txt'),
      command_complete = paste(command, file)
    ) %>% 
    select(command_complete) %>% 
    unlist(use.names = F)
  
  drem_script_location <- getwd() %>% 
    paste0('/scripts/DREM_', pipeline_input %>% extract2('analysis_type'), '.sh')
  
  drem_script <- drem_script_location %>% 
    file()
  
  writeLines(DREM_execute, drem_script)
  
  close(drem_script)
  
  previous_dir <- getwd()
  setwd(pipeline_input$DREM)
  
  message('Currently running DREM. \n')
  
  for(i in seq_along(DREM_execute)){
    current <- DREM_execute[[i]] %>% 
      str_extract('-b .*defaults.txt') %>% 
      basename() %>% 
      paste0('Generated defaults file: ', .)
    
    out <- DREM_execute[[i]] %>% 
      stri_reverse() %>% 
      gsub(pattern = ' .*', replacement = '') %>% 
      stri_reverse() %>%
      basename()
    
    message(current)
    message('DREM configuration file: ', out)  
    system(DREM_execute[[i]], show.output.on.console = T)
  }
  
  for(i in seq_along(DREM_execute)){
    current <- DREM_execute[[i]] %>% 
      str_extract('-b .*defaults.txt') %>% 
      basename()
    system(paste0('java -mx1024M -jar drem.jar'))
    
  }
  
  
  #if (exec){
  #  system(drem_script_location)
  #}
  setwd(previous_dir)
  return(tmp)
  
}




# Function name: load_read_data
# Purpose: Load count data. 
# Input: Directory with reads, experimental design
#        matrix, and flag.
# Output: Loaded read data. 
load_read_data <- function(map_dir, dm, dataset){
  if(dataset == 'Feature_counts'){
    SKIP <- 1
    TARGETS <- c(1, 2)
  } else{
    SKIP = 0
    TARGETS <- c('target_id', 'est_counts')
  }
  
  
  # manipulate read data to be input to DREM.
  reads <- list.files(map_dir, full.names = TRUE, recursive = T, pattern = '.txt$') %>% 
    str_replace(pattern = '//', '/') %>% 
    tibble(file = .) %>% 
    mutate(
      reads = map(file, ~read_tsv(.x, skip = SKIP, col_types = cols()) %>% dplyr::select(!!!TARGETS)
      )
    ) %>% 
    bind_cols(dm)
}



clean_target_id <- function(data){
  data %>% 
    select(target_id = 1, counts = 2) %>% 
    mutate(
      target_id = str_remove(target_id, '\\..*')
    )
}



# Function name: rearrange_count_data
# Purpose: Create data sets formatted for DREM.
# Input: Counts. 
# Output: Time series count data. 
rearrange_count_data <- function(reads, covars){
  
  reads2 <- reads %>% 
    tidyr::extract(condition, 
            into = str_split(pattern = ', ', covars) %>% unlist(),
            regex = '(.*) (.*) ([[:digit:]].*)') %>% 
    mutate(
      reads = purrr::map(reads,
                         clean_target_id
      )
    )
  
  
  time_series_spread <- reads2 %T>%
    { 
      # extract and sort numeric time series points.
      hour_vector <<- reads2 %>% 
        use_series(hour) %>% 
        unique() %>% 
        str_remove(pattern = '[A-Za-z]') %>% 
        as.numeric() %>% 
        sort.int() %>% 
        paste0('H', .) 
      
    } %>%
    
    # Match hour vector element structures to the hour column in the data. 
    mutate(
      hour = hour %>% 
        str_remove('[a-zA-Z]') %>% 
        str_replace('^', 'H')
    ) %>% 
    select(-file, -sample) %>% 
    
    # Associate each condition with a time series from the expression data. 
    unnest(reads) %>% 
    group_by(genotype, 
             condition, 
             which_replicate, 
             target_id) %>%
    spread(hour, 
           counts) %>%
    ungroup() %>% 
    select(genotype,
           condition, 
           which_replicate, 
           target_id, 
           hour_vector) %>% 
    group_by(genotype,
             condition, 
             which_replicate) %>% 
    nest() 
  
  return(time_series_spread)
}



# Function name: write_DREM_time_series_data
# Purpose: Write DREM-formatted data sets to files. 
# Input: Reformatted counts, pipeline input structure,
#        path to each individual reformatted count data
#        set.
# Output: Reformatted counts.
write_DREM_time_series_data <- function(time_series_spread, pipeline_input, drem_time_series_input_path){
  # create path to time series data. 
  time_series_spread <- time_series_spread %>% 
    mutate(
      file_name = paste(genotype, condition, which_replicate, sep = '_') %>% 
        paste0(drem_time_series_input_path, ., '.tsv')
    )
  # write time series counts to tsv files. 
  time_series_spread %$% 
    walk2(data, 
          file_name, 
          write_tsv)
  
  return(time_series_spread)
}





# Function name: write_default_files
# Purpose: Create DREM configuration file and write data.
# Input: Reformatted counts, and a configuration file template. 
# Output: Reformatted counts nested on genotype.
write_default_files <- function(reformatted_counts, defaults_template){
  
  
  # load original defaults file from DREM2.
  defaults_file <- defaults_template %>%
    read.delim(row.names = NULL) %>% 
    as_tibble() %>% 
    mutate_if(is.factor, as.character)
  
  
  nest_on_genotype <- reformatted_counts %>%
    
    select(genotype, condition, which_replicate, file_name) %>% 
    group_by(genotype, condition) %>% 
    nest() %>% 
    select(write_data = data) %>% 
    ungroup() %>% 
 
    mutate(
      new_default_file = map(write_data,
                             insert_DREM_input_to_defaults,
                             defaults_file
      ),
      new_default_file_name = pmap_chr(
        list(genotype, condition, write_data),  
        write_new_default_file_name
      )
    )
  
  nest_on_genotype %$% 
    walk2(new_default_file,
          new_default_file_name,
          write_tsv,
          col_names = FALSE)
  return(nest_on_genotype)
  
}





# Function name: insert_DREM_input_to_defaults
# Purpose: Create default files to execute DREM
#          for samples. 
insert_DREM_input_to_defaults <- function(condition_data, defaults_file){
  
  # Max replicates.
  total_replicates <- condition_data %>%
    select(which_replicate) %>%
    max()
  
  
  replicate_file <- condition_data %>% 
    ungroup() %>% 
    select(file_name) 
  
  defaults_file <- defaults_file %>% 
    select(a = 1, b = 2)
  
  defaults_file$b[3] <- replicate_file[1, 1] %>% 
    pull()
  
  defaults_file$b[12] <- replicate_file[c(2:total_replicates), 1] %>% 
    pull() %>% 
    str_c(collapse = ',')
  
  return(defaults_file %>% unnest(b))
}



# Function name: write_new_default_file_name
# Purpose: Create a file name to write DREM defaults script.
# Input: Sample characteristics (genotype, condition).
# Output: Tibble for collecting DREM.
write_new_default_file_name <- function(genotype, condition, drem_data){
  
  file_name <- drem_data %>% 
    ungroup() %>% 
    select(file_name) %>% 
    first() %>% 
    str_replace('.tsv$', 'defaults.txt') %>% 
    first()
  
  return(file_name)
}





# Function name: mapping_network_analysis
# Purpose: Apply network analysis to each comparison dataset. 
# Input: expr2.
# Output: Top genes or intersection with top GO results. 
mapping_network_analysis <- function(expr2){
  expr2 %>% 
    mutate(
      network_analsis = map2(clustering, 
                             GO_results,
                             network_analysis)
    )
}




# Function name: network_analysis
# Purpose: finds the network of differentialy expressed genes and collect the hubbed genes of each comaprisions
# Input: adjacency matrix, expression data, Go results , Dynamic modules
# Output: network graph, most important genes
network_analysis <- function(expr2, GO_results){
  
  m <- expr2$adjacency
  source_node = c()
  target_node = c()
  correlation <- c()
  genes_names <- names(expr2$datExpr)
  i <- genes_names[1:dim(m)[1]]
  j <- i
  
  
  for(gene in i)
  {
    for(gen in j)
    {
      if(m[gene,gen] < 0.9999 & m[gene,gen] > 0.3){
        source_node <- c(source_node, gene)
        target_node <- c(target_node, gen)
        correlation <- c(correlation, m[gene,gen])
      }
    }
  }
  
  NetworkData <- data.frame(source_node, target_node, correlation)
  
  
  network_dataframe <- data.frame(source_node,target_node,correlation)
  
  
  # cluster membership info from WCGNA
  nodes <- as.data.frame(cbind(genes_names,  expr2$dynamic_Modules))
  rownames(nodes)<-NULL
  colnames(nodes) <- c("genes","cluster")
  
  # igraph_network dataset
  igraph_network_dataframe <- network_dataframe
  igraph_network_dataframe <- aggregate(igraph_network_dataframe[,3], igraph_network_dataframe[,-3], sum)
  igraph_network_dataframe <- igraph_network_dataframe[order(igraph_network_dataframe$source_node, igraph_network_dataframe$target_node),]
  colnames(igraph_network_dataframe)[3] <- "weight"
  rownames(igraph_network_dataframe) <- NULL
  
  # igraph object:
  net <- graph.data.frame(igraph_network_dataframe, directed=F)
  
  
  # Generate colors base on clustering done by WGCNA:
  nodeIndex<-c()
  for(name in V(net)$name){
    nodeIndex<-c(nodeIndex,which(nodes$genes==name))
  }
  colors <- expr2$dynamic_Colors[nodeIndex]
  
  
  V(net)$color <- colors
  
  
  degree_nodes <- igraph::degree(net, mode ="total")
  V(net)$size <- degree_nodes
  
  
  #ching arrow width , label color, edge width etc.
  V(net)$label.color <- "black"
  E(net)$width <- E(net)$weight*2
  
  E(net)$arrow.size <- 10
  edge_width <-(E(net)$weight-mean(E(net)$weight))*5+1
  
  #ploting the graph
  tkplot(net,edge.arrow.size=1,edge.curved=0,edge.width=edge_width,edge.color="gray80")
  
  #taking the cluster color
  cluster_colors<-unique(V(net)$color)
  
  #making an empty list of hubbed genes
  genes_hubbed<-list()
  
  #finding out the node centrality score of all the clusters and therefore finding out the hubbed genes 
  for(i in 1: length(unique(V(net)$color))){
    #finding nodes for the particlular cluster
    nodes_of_interest <- V(net)[which(V(net)$color == cluster_colors[i])]
    #finding the subgraph for the cluster nodes
    selgraph <- ego(net, order = 1, nodes =nodes_of_interest , mode = "all",
                    mindist = 0)
    subgraph_cluster <- induced_subgraph(net,unlist(selgraph))
    
    #find the centrality score of each node in that cluster
    centarlity_score <- hub_score(subgraph_cluster , scale = TRUE, weights = NULL,
                                  options = arpack_defaults)
    
    #find out the hubbed genes based on centrality score
    hubbed_genes <-names(centarlity_score$vector[which(centarlity_score$vector >= quantile(centarlity_score$vector,.95))])
    hubbed_genes<-as.data.frame(hubbed_genes)
    genes_hubbed[i]<-hubbed_genes
    
    
    rm(hubbed_genes)
    rm(nodes_of_interest)
    rm(selgraph)
    rm(selegoG)
    rm(centarlity_score)
    
  }
  
  hubbed_genes_in_total <- as.data.frame(levels(unlist(genes_hubbed)))
  colnames(hubbed_genes_in_total)<-"hubbed_genes"
  
  #selecting important genes on the basis of node degree
  top_node_dgree_genes<-as.data.frame(names(V(net)[which(V(net)$size>=quantile(V(net)$size,0.95))]))
  colnames(top_node_dgree_genes)<-"hubbed_genes"
  
  #merging the important genes from the node degree and hubbed genes and finding top genes from Network Analysis 
  top_genes<-merge(x = hubbed_genes_in_total, y = top_node_dgree_genes, by = "hubbed_genes", all = TRUE) %>% 
    unlist(use.names = FALSE)
  
  #finding the TOP hubbed genes frm the GO results and network anlysis 
  hub_go_intersection <- intersect(GO_results, top_genes)
  
  if(length(hub_go_intersection) > 0){
    print('Hub genes intersected with top GO genes.')
    return(hub_go_intersection)
  } else{
    print('Hub genes do not intersect with top GO genes.')
    return(top_genes)
  }
}