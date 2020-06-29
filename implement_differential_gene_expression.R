# Function name: implement_differential_gene_expression
# Purpose: Extract alignment function output and execute differential gene 
#          expression testing depending on the upstream alignment or pseudo-alignment 
#          software used.  
# Input: Alignment function results.   
# Output: Alignment function results and processed configuration file. 
implement_differential_gene_expression <- function(alignment_results){

  # Separate alignment results (tibble with tool commands 
  # and directories) from the pipeline input (configuration 
  # file exraction).
  alignment <- alignment_results %>% 
    first()
  
  pipeline_input <- alignment_results %>% 
    last()
  
  # Load experimental design matrix file. 
  design_matrix <- pipeline_input %>% 
    extract2('design_matrix') %>% 
    read_csv(col_types = cols())
  
  # Separate covariates for column selection.
  sample_covariates <- pipeline_input %>% 
    extract2('sample_covariates') %>% 
    split_and_unlist_conditions()
  
  # Map design matrix contents with sets of samples to 
  # compare in differential expression testing. 
  sample_comparisons <- pipeline_input %>% 
    extract2('sample_comparisons_file') %>% 
    read_tsv(col_names = FALSE, col_types = cols()) %>% 
    extract_sample_comparison_sets(design_matrix, sample_covariates)
  
  kallisto <- alignment %>% 
    extract2('key') %>% 
    str_detect('kallisto') %>% 
    any()
  
  STAR <- alignment %>% 
    extract2('key') %>% 
    str_detect('star') %>% 
    any()
  
  # Sleuth non-functional.
  #if(kallisto){
  #  implement_sleuth(alignment, pipeline_input, design_matrix, sample_comparisons, sample_covariates)
  #}
  
  if(STAR){
    implement_edgeR(alignment, pipeline_input, design_matrix, sample_comparisons)
  }
  
}



# Function name: split_and_unlist_conditions
# Purpose: Separate covariates from list of conditions. 
# Input: Condition names from experimental design matrix and
#        configuration file. 
# Output: Split list with spaces removed. 
split_and_unlist_conditions <- function(covariate_character_vector){
  covariate_character_vector %>% 
    str_split(pattern = ',') %>% 
    unlist(use.names = FALSE) %>% 
    str_remove('[:space:]')
}




# Function name: extract_sample_comparison_sets
# Purpose: Associate sample sets with sample comparisons.
# Input: Sample comparisons file, experimental design matrix,
#        and the split sample conditions. 
# Output: Comparison sets associated with file names. 
extract_sample_comparison_sets <- function(comparison_sets, design_matrix, sample_conditions){
  comparison_sets %>%
    
    # Combine sample identifiers in sample comparisons file with 
    # a pipe for searching the experimental design matrix with. 
    unite(col = samples, 
          sep = '|') %>%
    
    transmute(
      comparison_set = map(samples, 
                           ~filter(design_matrix, 
                                   str_detect(!!sym('sample'), .x)
                           )
      ),
      
      file_name = map_chr(comparison_set, 
                          ~write_differential_testing_file_name(comparison_set = .x, 
                                                                sample_conditions)
      )
    )
}




# Function name: write_differential_testing_file_name
# Purpose: Combine conditions from each sample comparison
#          set to create a file name with.
# Input: Comparison set, and split conditions. 
# Output: A file name with '_vs_' between conditions compared. 
write_differential_testing_file_name <- function(comparison_set, sample_conditions){
  comparison_set %>%
    select(condition) %>% 
    unique() %>% 
    arrange(desc(condition)) %>% 
    unlist(use.names = FALSE) %>% 
    str_c(collapse = '_vs_') %>% 
    str_replace_all(pattern = ' ', replacement = '_')
}




# Function name: implement_sleuth
# Purpose: Run Sleuth on one sample comparison set at a time. 
# Input: Alignment results data structure, processed configuration file, 
#        experimental design matrix, sample comparison sets, and split
#        sample comparisons. 
implement_sleuth <- function(alignment, pipeline_input, design_matrix, sample_comparisons, sample_covariates){
  
  # Provide access to Kallisto directories.
  kallisto_directories <- alignment %>% 
    filter(key == 'kallisto') %>% 
    select(directories) %>% 
    unlist(use.names = F) 
  
  # List sample quantification results.
  quantification_files <- kallisto_directories %>% 
    str_subset('Kallisto_quantifications') %>% 
    list.files(full.names = TRUE) %>% 
    tibble(path = .) %>% 
    mutate(
      # Ensure no double forward-slash.
      path = map_chr(path, 
                     str_replace,
                     '//',
                     '/'),
      # Give sample name. 
      sample = map_chr(path, 
                       str_remove,
                       '.*/'
      )
    )
  
  # Reset quantification file paths in design matrix and comparison sets 
  # from the original design matrix to the sample quantification results.
  design_matrix <- design_matrix %>% 
    mutate(
      path = quantification_files %>%
        select(path) %>% 
        unlist(use.names = F)
    )
  
  # 
  sample_comparisons <- sample_comparisons %>%
    unnest(cols = c(comparison_set)) %>% 
    left_join(quantification_files, by = 'sample') %>% 
    select(-path.x) %>% 
    dplyr::rename(path = path.y) %>% 
    nest(data = c(sample, genotype, condition, hour, which_replicate, testing_condition, path)) %>% 
    dplyr::rename(comparison_set = data) %>% 
    mutate(
      testing_condition = map(comparison_set, 
                              structure_covariates_for_DGE_testing,
                              sample_covariates
      ),
      comparison_set = map(comparison_set, 
                           separate_conditions,
                           sample_covariates
      )
    )
  
  # Determine which mode(s) to execute Sleuth in.
  gene_mode <- pipeline_input %>% 
    extract2('sleuth_gene_mode')
  
  transcript_mode <- pipeline_input %>% 
    extract2('sleuth_transcript_mode')
  
  sleuth_mode <- gene_mode %>% 
    select_sleuth_mode(transcript_mode)
  
  # Place Sleuth results next to Kallisto.
  sleuth_path <- kallisto_directories %>% 
    str_subset(pattern = 'Kallisto/$') %>% 
    str_replace('Kallisto', 'Sleuth')
  
  transcript_level <- str_detect(sleuth_mode, 'transcript')
  
  gene_level <- str_detect(sleuth_mode, 'gene')
  
  transcript_and_gene_level <- str_detect(sleuth_mode, 'transcript_and_gene_level')
  
  analysis_type <- pipeline_input %>% 
    extract2('analysis_type')
  
  sleuth_level_directories <- tibble(sleuth_path, transcript_level, gene_level, transcript_and_gene_level) %>% 
    mutate(
      transcript_level = case_when(transcript_level ~ paste0(sleuth_path, analysis_type, '/transcript_level/')
      ),
      gene_level = case_when(gene_level ~ paste0(sleuth_path, analysis_type, '/gene_level/')
      )
    ) %>% 
    gather() %>% 
    filter(
      str_detect(value, paste0('Sleuth/', analysis_type)
      )
    ) %>% 
    mutate(
      write_directory = map(value, 
                            dir.create,
                            recursive = TRUE,
                            showWarnings = FALSE
      )
    ) %>% 
    select(-write_directory) 
  
  if(transcript_and_gene_level | gene_level){
    target_mapping <- design_matrix %>%
      slice(1) %>% 
      gene_level_analysis_data_structure() %>% 
      select(ttg) %>%
      unnest()
  }
  
  sample_comparisons <- sample_comparisons %>%  
    mutate(
      testing_condition_formula = map(testing_condition, 
                                      extract_single_testing_condition
      )
    )
  
  if(transcript_and_gene_level){
    transcript_comparisons <- sample_comparisons %>% slice(1:2) %>% 
      transcript_prep() 
    
    gene_comparisons <- sample_comparisons %>% slice(1:2) %>% 
      gene_prep(target_mapping)
    
    comparisons <- bind_rows(transcript_comparisons, gene_comparisons)
  }
  
  if(gene_level & !(transcript_and_gene_level | transcript_level)){
    comparisons <- sample_comparisons %>% 
      gene_prep(target_mapping)
  }
  
  if(transcript_level & !(gene_level | transcript_and_gene_level)){
    comparisons <- sample_comparisons %>% 
      transcript_prep()
  }
  
  
  rm(sample_comparisons)
  
  comparisons <- comparisons %>% 
    left_join(sleuth_level_directories) %>% 
    dplyr::rename(level_directory = value) %>%
    select(-key)
  
  comparisons <- comparisons %>% 
    mutate(
      pca = map2(prep, 
                 testing_condition_formula, 
                 run_sleuth_pca
      ),
      pca_save = pmap(
        list(pca, level_directory, file_name), 
        save_sleuth_pca
      ),
      fit = map2(prep, 
                 testing_condition_formula, 
                 run_sleuth_fit
      ),
      lrt = map(fit, 
                run_sleuth_lrt
      ),
      lrt_results = map(lrt, 
                        extract_sleuth_lrt_results
      ),
      save_lrt_results = pmap(
        list(lrt_results, level_directory, file_name), 
        save_sleuth_results
      )
    ) 
}




# Function name: structure_covariates_for_DGE_testing
# Purpose: Search for conditions in the comparison sets which vary in a column.
# Input: Comparison set and split sample conditions.
# Output: Unique key, value pairs for conditions. 
structure_covariates_for_DGE_testing <- function(comparison_set, sample_conditions){
  condition_count <- sample_conditions %>% 
    length()
  
  comparison_set <- comparison_set %>% 
    extract(col = condition, 
            into = sample_conditions, 
            regex = rep('(.*)', times = condition_count) %>% 
              paste(collapse = ' ')
    ) %>% 
    select(sample_conditions) %>% 
    distinct()
  
  row_num <- comparison_set %>% 
    nrow()
  
  comparison_set %>% 
    gather() %>% 
    group_by(key) %>% 
    mutate(
      variant = value %>% n_distinct
    ) %>% 
    filter(variant > 1) %>% 
    mutate(
      test = map2_chr(key, 
                      value, 
                      paste0)
    ) %>% 
    select(key, 
           test
    ) %>% 
    distinct()  
}




# Function name: separate_conditions
# Purpose: Separate condition column into multiple columns  based on the covariates.
# Input: Sample comparison set and split condition list.
# Output: A sample comparison set with one column for each covariate. 
separate_conditions <- function(clean_design_matrix, conditions){
  condition_count <- conditions %>% 
    length()
  
  clean_design_matrix %>% 
    tidyr::extract(col = condition, 
                   into = conditions, 
                   regex = rep('(.*)', times = condition_count) %>% 
                     paste(collapse = ' ')
    )
}



# Function name: select_sleuth_mode 
# Purpose: Determine  at what level to run Sleuth: gene, transcript, or both. 
# Input: Boolean values for gene mode and transcript mode. 
# Output: A mode to run Sleuth in.
select_sleuth_mode <- function(gene_mode, transcript_mode){
  if_else(condition = (gene_mode == 'TRUE'), 
          true = if_else(condition = (transcript_mode == 'TRUE'),
                         true = 'transcript_and_gene_level',
                         false = 'gene_level'),
          false = if_else(condition = (transcript_mode == 'TRUE'),
                          true = 'transcript_level',
                          false = 'NULL')
  )
}




# Function name: gene_level_analysis_data_structure
# Purpose: Removes splicing notation from one quantification file result to use for gene level Sleuth analysis.
# Input: Experimental design matrix. 
# Output: Gene to transcript mapping. 
gene_level_analysis_data_structure <- function(design_matrix){
  design_matrix <- design_matrix %>%
    mutate(
      path = map_chr(path, 
                     list.files,
                     '^abundance.tsv$',
                     full.names = TRUE, 
                     recursive = T
      ),
      reads = map(path, 
                  read_tsv,
                  col_types = cols())
    ) %>% 
    unnest()
  
  # Check if the gene names have been extracted from 
  # transcripts for the first quantification dataset. 
  gene_column_existent <- design_matrix %>% 
    colnames() %>% 
    str_detect('gene') %>% 
    any()
  
  if(!gene_column_existent){
    design_matrix <- design_matrix %>% 
      group_by(target_id) %>% 
      nest() %>% 
      mutate(
        gene = target_id %>% 
          str_remove('[[:punct:]].*')
      ) %>% 
      unnest()
    
  } else{
    design_matrix <- design_matrix
  }
  
  design_matrix %>% 
    group_by(sample, 
             path, 
             condition, 
             which_replicate
    ) %>% 
    nest() %>% 
    mutate(
      data = map2(data, 
                  path, 
                  write_tsv
      ),
      ttg = map(data, 
                select,
                gene, 
                target_id
      ), 
      path = map_chr(path, 
                     dirname)
    ) %>% 
    select(-data)
}




# Function name: extract_single_testing_condition
# Purpose: Select a single covariate to test with by 
#          identifying the condition column with varying 
#          contents across samples.
# Input: Split sample conditions. 
# Output: Covariate which varies in its contents.
extract_single_testing_condition <- function(sample_conditions){
  
  count_comparisons <- sample_conditions %>%
    select(key) %>%
    plyr::count()

  key_no_variance <- count_comparisons %>%
    filter(freq < 2) %>%
    select(key)

  if(nrow(key_no_variance) >=1 ){
    no_variance_warning <- key_no_variance %>%
      unlist(use.names = FALSE) %>%
      str_c(collapse = ' and ') %>%
      paste0('Warning: Removing ', ., '. Sleuth needs two comparisons in each category.')
    
    message(no_variance_warning)
  }

  sample_conditions %>%
    anti_join(key_no_variance) %>%
    ungroup() %>%
    select(key) %>%
    unlist(use.names = FALSE) %>%
    paste0('~', .) %>% 
    unique()
}




# Function name: transcript_prep
# Purpose: A Sleuth_prep implementation specific
#          to transcript analaysis. 
# Input: A sample comparison set. 
# Output: Sleuth prep object. 
transcript_prep <- function(comparisons){
  comparisons %>% 
    mutate(
      prep = map2(comparison_set, 
                  testing_condition_formula,
                  ~sleuth::sleuth_prep(sample_to_covariates = .x, 
                                       full_model = .y %>% 
                                         unique() %>% 
                                         as.formula(), 
                                       gene_mode = FALSE)
      ),
      key = 'transcript_level'
    )
}





# Function name: gene_prep
# Purpose: A Sleuth_prep implementation specific to gene analaysis. 
# Input: A sample comparison set.
# Output: Sleuth prep object. 
gene_prep <- function(comparisons, mapped_targets){
  comparisons %>% 
    mutate(
      testing_condition_formula = map(testing_condition_formula,
                                      as.formula
      ),
      prep = map2(comparison_set, 
                  testing_condition_formula, 
                  sleuth::sleuth_prep,
                  target_mapping = mapped_targets, 
                  aggregation_column = 'gene', 
                  extra_bootstrap = TRUE, 
                  gene_mode = TRUE
      ),
      key = 'gene_level'
    )
}




# Function name: run_sleuth_pca
# Purpose: Plot PCA plots for a comparison.
# Input: Sleuth prep object and a covariate to 
#        color with. 
# Output: PCA plot for comparison. 
run_sleuth_pca <- function(prep, covariate_color){
  covariate_color <- covariate_color %>% 
    str_remove('^~')
  
  prep %>%
    sleuth::plot_pca(color_by = covariate_color, text_labels = FALSE)
}




# Function name: save_sleuth_pca
# Purpose: Save PCA plots for a comparison.
# Input: PCA plot to save, Parent directory, and file name.
# Output: ggsave object.
save_sleuth_pca <- function(pca_plot, level_path, file_name){
  level_path <- level_path %>% 
    basename() %>% 
    paste0(level_path, ., '_PCA')
  
  dir.create(level_path, showWarnings = FALSE, recursive = TRUE)
  
  pca_plot %>% 
    ggsave(filename = paste0(file_name, '.pdf'), path = level_path, device = 'pdf')
}




# Function name: run_sleuth_fit
# Purpose: Measurement error model fitting with Sleuth.
# Input: Sleuth prep object and testing formula. 
# Output: Sleuth fit object.
run_sleuth_fit <- function(prep, formula_string){
  prep %>% 
    sleuth::sleuth_fit(obj = ., 
                       formula = formula_string %>% 
                         as.formula, 
                       'full'
    ) %>% 
    sleuth::sleuth_fit(obj = ., 
                       formula = ~1, 
                       'reduced')
}




# Function name: run_sleuth_lrt
# Purpose: Perform Likelihood Ratio Test.
# Input: Sleuth fit object.
# Output: LRT results. 
run_sleuth_lrt <- function(fit){ 
  fit %>% 
    sleuth::sleuth_lrt(obj = ., 
                       'reduced', 
                       'full')
}




# Function name: extract_sleuth_lrt_results
# Purpose: Extract Likelihood Ratio test results from a sleuth object.
# Input: Sleuth object
# Output: Likelihood ratio test results.
extract_sleuth_lrt_results <- function(lrt){
  lrt %>% 
    sleuth::sleuth_results(test = 'reduced:full', 
                           test_type = 'lrt')
}




# Function name: save_sleuth_results
# Purpose: Save LRT results to a file for a comparison.
# Input: LRT results, analysis level, file name.
# Output: None.
save_sleuth_results <- function(results, level_path, file_name){
  level_path <- level_path %>% 
    basename() %>% 
    paste0(level_path, ., '_results/')
  dir.create(level_path, showWarnings = FALSE, recursive = TRUE)
  
  results %>% 
    write_tsv(path = paste0(level_path, file_name, '.tsv'))
}




# Function name: implement_edgeR
# Purpose: Conduct differential testing on STAR
#          mapped reads with edgeR.
# Input: Alignment decision data structure, processed configuration
#        file, experimental design matrix, and sample comparison sets.
implement_edgeR <- function(alignment, pipeline_input, design_matrix, sample_comparisons){
  STAR_counts_path <- alignment %>%  
    filter(key == 'star') %>%  
    select(directories) %>% 
    unlist(use.names = FALSE) %>% 
    str_subset(pattern = 'Feature_counts') %>% 
    list.files(full.names = TRUE) %>% 
    str_subset(pattern = 'clean_counts.txt$')
  
  design_matrix <- design_matrix %>%
    mutate(
      path = STAR_counts_path
    )
  
  no_return <- STAR_counts_path %>%
    tibble(files = .) %>%
    mutate(
      # Select the column one and column two from feature counts data; to be ran with submitted data.
      counts_data = map(files, ~read_tsv(.x, col_names = TRUE, skip = 1, col_types = cols()) %>% select(1, 2)
      )
    ) %>%
    filter_counts_data(design_matrix = design_matrix) %>%
    edgeR_preliminary(count_keeps = ., samplename = design_matrix, edgeR_directory = STAR_counts_path, comparisons = sample_comparisons, pipeline_input)
}




# Function name: filter_counts_data
# Purpose: filter the lowly expressed genes for edgeR analysis.
# Input: counts, experimental design matrix.
# Output: data frame of counts for each sample.
filter_counts_data <- function(counts, design_matrix){
  data <- data.frame(counts$counts_data)
  
  # storing data as data frame
  total_samples <- nrow(counts)*2
  ans <- seq(2,total_samples,2)
  data <- data[c(1, ans)]
  
  # Load the sample information and include the header. 
  samplename <- design_matrix
  
  #(optional) naming the samples
  #d=colnames(sleuth_table)[1]
  d="target_id"
  x <- samplename$`sample`
  a <- c(d,x)
  colnames(data) <- a
  rowname<-data[c(1)]
  
  #filterng the data 
  keep <- rowSums(cpm(data[,2:length(data)]) > 0.5) >= 2
  data[keep,]
}




# Function name: edgeR_preliminary
# Purpose: process input data for Edge R analysis followed by the analysis and set the directory where data has to be saved 
# Input: kept counts, a sample name, edgeR directory, 
#        comparison matrix, and the extracted configuration
#        data structure. 
edgeR_preliminary <- function(count_keeps, samplename, edgeR_directory, comparisons, pipeline_input){
  
  edgeR_output_directory <- edgeR_directory %>% 
    first() %>% 
    str_replace(pattern = 'STAR.*', replacement = 'edgeR/') %>% 
    paste0(pipeline_input$analysis_type, '/')
  
  dir.create(edgeR_output_directory, recursive = TRUE, showWarnings = FALSE)
  
  # grouping the samples
  edgeR_comparision_set1 <- gsub(" ","_",samplename$condition)
  group <- factor(edgeR_comparision_set1)
  y <- DGEList(count_keeps[,2:length(count_keeps)],group=group,genes=count_keeps[c(1)])
  logcounts <- cpm(y,log=TRUE)
  y <- calcNormFactors(y, method='TMM')
  logcounts <- cpm(y,log=TRUE)
  #########################################################
  # forming the design matrix from DGE list generated before
  design.mat <- model.matrix(~ 0 + y$samples$group)
  
  colnames(design.mat) <- levels(y$samples$group)
  v <- voom(y,design.mat)
  #fitting the model
  vfit <- lmFit(v, design.mat)
  ###############contrast matrix1#######################################
  #can be of any choice
  comparisons$file_name %>%
    tibble(C = .) %>%
    mutate(
      B = map_chr(C, ~str_remove(string = .x, pattern = '_vs.*')),
      A = map_chr(C, ~str_remove(string = .x, pattern = '.*vs_')),
      all = pmap(list(C, A, B), ~paste0(..1, ' = ', ..2, ' - ', ..3))
    ) %>% 
    mutate(
      make_contrasts = map(all, ~makeContrasts(contrasts = .x, levels = group)),
      vfit = map(make_contrasts, ~contrasts.fit(fit = vfit, contrasts = .x)),
      tfit = map(vfit, ~treat(fit = .x, lfc = pipeline_input$edger_lfc %>% as.numeric())),
      dt = map(tfit, ~decideTests(object = .x, method = pipeline_input$edger_method, adjust.method = pipeline_input$edger_adjustment_method, p.value = pipeline_input$significance_cutoff)),
      top_tables = map2(tfit, all, ~topTable(fit = .x, coef = .y, sort.by = 'p', n = 'Inf')),
      write_out = map2(top_tables, C, ~write_tsv(x = .x, path = paste0(edgeR_output_directory, .y, '.tsv')))
    )
}
