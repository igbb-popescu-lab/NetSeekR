# Function name: implement_alignment
# Purpose: Extract all arguments from the configuration file, decide which 
#          pipeline(s) to run, assemble arguments and implement each 
#          decided pipeline.
# Input: Arguments file path.
# Output: Alignment decision data structure and processed configuration file. 
implement_alignment <- function(arguments_file, execute_script){
  
  # Convert arguments from configuration file 
  # into a named list data structure.
  pipeline_input <- arguments_file %>% 
    extract_pipeline_input_from_configuration()
  
  # Determine which pipeline(s) to execute, 
  # and write directory trees for each tool
  # chosen.
  alignment_decision <- pipeline_input %>% 
    decide_alignment_tool() 
  
  # Extract paths to all written directories.
  alignment_directories <- alignment_decision %>% 
    select(directories) %>% 
    unlist(use.names = FALSE)
  
  # Map over elements and create each directory.
  alignment_directories %>% 
    map(dir.create, 
        recursive = TRUE, 
        showWarnings = FALSE)
  
  # Subset arguments from the configuration file
  # that are specific to an alignment tool.
  alignment_decision <- alignment_decision %>% 
    
    mutate(
      arguments = map(key, 
                      assemble_alignment_arguments,
                      pipeline_input
      )
    )
  
  # Set Boolean values for tools to be executed.
  Kallisto <- alignment_decision %>% 
    align_tool_bool('key', 'kallisto')
  
  STAR <- alignment_decision %>% 
    align_tool_bool('key', 'star')
  
  # Implement the selected alignment tool(s).
  
  # Kallisto 
  if(Kallisto){
    # Subset arguments specific to Kallisto.
    alignment_decision %>% 
      subset_tool_arguments('kallisto') %>% 
      implement_kallisto(Kallisto_arguments = ., directories = alignment_directories, execute_script)
  }
  
  # STAR
  if(STAR){
        
    # Subset arguments
    args <- alignment_decision %>% 
      filter(str_detect(key, 'star')) %>% 
      select(arguments) %>% 
      unnest(arguments)
    
    # Write bash scripts for STAR pipeline tools (STAR, Feature counts).
    star <- implement_STAR(STAR_arguments = args, directories = alignment_directories, execute_script)
    
    # feature counts
    feature_counts_output <- star %>% 
      str_replace('genome_Dir', 'Feature_counts')
    
    dir.create(feature_counts_output, showWarnings = FALSE)
    
    gtf <- pipeline_input %>% 
      extract2('star_sjdbGTFfile')
    
    path_to_feature_counts <- pipeline_input %>% 
      str_subset('.*bin/featureCounts.*')
    
    implement_feature_counts(output_dir = feature_counts_output, annotation = gtf, feature_counts = path_to_feature_counts, genome_Dir = star, execute_script)
  }
  return(list(alignment_decision, pipeline_input))
}




# Function name: extract_pipeline_input_from_configuration
# Purpose: Create a named list data structure for storing 
#          and accessing arguments from the configuration file. 
# Input: Arguments file path.
# Output: Arguments to pipeline in a named list. 
extract_pipeline_input_from_configuration <- function(arguments_file){
  
  # Remove non-existant configurations. 
  arg_file <- arguments_file %>% 
    read_tsv(col_types = cols()) %>% 
    filter(!(Argument_Type %>% is.na & Argument %>% is.na)) 
  
  # Extract argument types. 
  arg_type <- arg_file %>% 
    select(Argument_Type) %>% 
    unlist(use.names = F)
  
  # Link argument to argument type in a named list data structure.
  arg_file %>% 
    select(Argument) %>% 
    unlist(use.names = F) %>% 
    as.list() %>% 
    set_names(arg_type)
}




# Function name: decide_alignment_tool
# Purpose: Write directory trees for decided tool.  
# Input: processed configuration file (arguments list).
# Output: Tibble from which decisions on which downstream pipeline
#         decisions are to be made. 
decide_alignment_tool <- function(pipeline_input){
  
  tibble(kallisto = pipeline_input %>% extract2('kallisto'), 
         star = pipeline_input %>% extract2('star')
  ) %>% 
    
    gather() %>% 
    
    # Remove alignment tool which is not selected for implementing. 
    filter(value == 'TRUE') %>% 
    mutate(
      
      directories = map2(key, value, 
                         
                         # Determine which directory tree(s) to create.
                         ~dplyr::case_when(.x == 'kallisto' ~ pipeline_input %>% 
                                      extract2('kallisto_output_dir') %>% 
                                      write_kallisto_directory_tree,
                                    
                                    .x == 'star' ~ pipeline_input %>% 
                                      extract2('star_genomeDir') %>% 
                                      write_STAR_directory_tree
                         )
      )
    )
}




# Function name: write_kallisto_directory_tree
# Purpose: Create directory tree for Kallisto results.
# Input: Alignment pipeline-specific output directory. 
# Output: Kallisto path list. 
write_kallisto_directory_tree <- function(output_directory){
  
  # Attach output directories to the Kallisto directory created.
  Kallisto_paths <- tibble(kallisto = output_directory) %>%
    
    mutate(
      
      kallisto_quantifications = kallisto %>% 
        str_replace(pattern = '$', replacement = 'Kallisto_quantifications/'
        ),
      
      kallisto_log_files = kallisto %>% 
        str_replace(pattern = '$', replacement = 'Kallisto_log_files/')
      
    ) %>% 
    gather() %>% 
    select(value) %>% 
    list()
}




# Function name: write_STAR_directory_tree
# Purpose: Create directory tree for STAR results.
# Input: Alignment pipeline-specific output directory. 
# Output: STAR path list. 
write_STAR_directory_tree <- function(output_directory){
  
  # Save output directory basename to replace with pattern matching.
  STAR_node <- output_directory %>% 
    basename()
  
  STAR_paths <- tibble(STAR = output_directory) %>%
    
    mutate(
      
      star_log_files = STAR %>% 
        str_replace(pattern = STAR_node, replacement = 'STAR_log_files/'
        ),
      
      feature_counts = STAR %>% 
        str_replace(pattern = STAR_node, replacement = 'Feature_counts/')
      
    ) %>% 
    gather() %>% 
    select(value) %>% 
    list()
}




# Function name: assemble_alignment_arguments
# Purpose: Subset arguments particular to an alignment tool 
#          and ensure proper command formatting. 
# Input: Alignment tool name and processed configuration file. 
# Output: Formatted arguments for STAR or Kallisto in a tibble
#         data structure.
assemble_alignment_arguments <- function(alignment, pipeline_input){
  
  # These commands need an equals sign after the command. 
  k_equals <- 'index|bootstrap_samples|seed|fragment_length|sd|threads|output_dir|kmer_size' %>% 
    str_replace_all(pattern = '\\|', replacement = '|kallisto_') %>% 
    str_replace(pattern = '^', 'kallisto_') %>% 
    str_split(pattern = '\\|') %>% 
    unlist()
  
  # Subset arguments based on which tools they belong. 
  arguments <- pipeline_input[grep(pattern = alignment, x = names(pipeline_input))] %>% 
    as_tibble() %>% 
    gather() %>% 
    drop_na() %>% 
    filter(key != alignment) %>% 
    
    # Place equal signs in command where appropriate.
    mutate_at(
      vars(key), 
      list(
        ~(dplyr::if_else(condition = .x %>% is_in(k_equals), 
                  true = .x %>% paste0('='), 
                  false = .x)
        )
      )
    ) %>% 
    mutate(
      # Reformat pipeline arguments to be usable by the alignment tool in the command line. 
      key = key %>% 
        str_replace(pattern = paste0(alignment, '_'), 
                    replacement = '--')
    )
}




# Function name: align_tool_bool
# Purpose: Check for tool existence in the decision data structure.
# Input: Alignment decision tibble, key with alignment tool name,
#        and a string for the alignment name. 
# Output: Boolean value.
align_tool_bool <- function(tib, column_name, str_to_detect){
  tib %>% 
    select(column_name) %>% 
    unlist() %>% 
    str_detect(str_to_detect) %>% 
    any()
}




# Function name: subset_tool_arguments
# Purpose: Subset tool arguments from decision data structure 
#          and pass them to tool-specific command processing functions.    
# Input: Alignment decision data structure, tool name.
# Output: Tibble of arguments from selected tools. 
subset_tool_arguments <- function(alignment_decision, tool_name){
  alignment_decision %>%
    filter(
      str_detect(key, !!tool_name)
    ) %>% 
    select(arguments) %>% 
    unnest(cols = arguments)
}




# Function name: implement_kallisto
# Purpose: Write bash scripts for Kallisto index building and quantification,
#          then execute them. 
# Input: Arguments to Kallisto from configuration file, and all directories
#        created. 
implement_kallisto <- function(Kallisto_arguments, directories, execute_script){
  
  # Separate indexing from quantifying arguments.  
  index_arguments <- 'index|fasta|kmer|unique'
  
  organize_kallisto_arguments <- Kallisto_arguments %>% 
    
    mutate(
      key = map_chr(key, 
                    str_replace_all,
                    '_', 
                    '-'
      ),
      index = map_lgl(key,
                      str_detect,
                      index_arguments
      ),
      quant = map_lgl(key, 
                      str_detect,
                      index_arguments, 
                      negate = TRUE
      )
    ) %>% 
    
    # Mark arguments needed for both indexing and quantifying. 
    mutate_at(
      vars(quant, index), 
      list(
        ~dplyr::if_else(condition = str_detect(key, 'index'), true = TRUE, false = .x)
      )
    ) %>% 
    mutate_at(
      vars(index), 
      list(
        ~dplyr::if_else(condition = str_detect(key, 'path'), true = TRUE, false = .x)
      )
    )
  
  # Extract output directory for quantifications.
  output_dir <- organize_kallisto_arguments %>% 
    filter(
      str_detect(key, 'output')
    ) %>% 
    select(value) %>% 
    unlist(use.names = FALSE)
  
  # Extract output directory for log files.
  log_file_dir <- directories %>% 
    str_subset('Kallisto_log_files')
  
  # Subset indexing-building arguments.
  arguments_to_build_index <- organize_kallisto_arguments %>% 
    filter(index == TRUE) %>% 
    select(key, value)
  
  # Use fasta pattern to move fasta reference to end of command. 
  fasta <- arguments_to_build_index %>% 
    filter(
      str_detect(key, 'fasta-files')
    ) %>% 
    select(value) %>% 
    unlist(use.names = FALSE)
  
  # Identify non-optional arguments.
  arg_core <- arguments_to_build_index %>% 
    filter(
      str_detect(key, pattern = 'path|fasta|index')
    )
  
  # Identify optional arguments using already identified non-optional arguments.
  extraneous_args <- anti_join(arguments_to_build_index, arg_core, by = c('key', 'value')) %>% 
    spread(key, value)
  
  # Assemble optional arguments.
  extraneous_args_type <- extraneous_args %>% 
    colnames()
  
  extraneous_args_arg <- extraneous_args %>% 
    unlist(use.names = FALSE)
  
  extraneous_args <- rbind(extraneous_args_type, extraneous_args_arg) %>% 
    str_c(collapse = ' ') %>% 
    str_remove(' TRUE')
  
  # Assemble non-optional arguments.  
  arg_core <- arg_core %>% 
    spread(key, value)
  
  arg_core_type <- arg_core %>% 
    colnames()
  
  arg_core_arg <- arg_core %>% 
    unlist(use.names = FALSE)
  
  # Assemble index-building shell script. 
  build_index <- rbind(arg_core_type, arg_core_arg) %>% 
    str_c(collapse = ' ') %>% 
    
    # Move path to Kallisto the the start of the command. 
    str_replace(pattern = '^', 
                replacement = str_match(string = ., pattern = '--path.*') %>% paste(' ', sep = '')
    ) %>% 
    
    # Remove unnecessary arguments used as flags before.
    str_remove(pattern = ' --path.*'
    ) %>%
                  
    str_replace('--fasta-files', replacement = str_match(string = ., pattern = '--index.*'
    ) %>% 
      str_replace(pattern = '= ', replacement = '=')
    ) %>% 
    
    str_replace('--index=', 'index --index='
    ) %>% 
    
    str_remove(' --index= .*'
    ) %>% 
    
    paste(extraneous_args, collapse = ''
    ) %>% 
    
    str_replace(pattern = '= ', replacement = '='
    ) %>% 
    
    str_replace(pattern = '$', replacement = str_match(string = ., pattern = fasta)
    ) %>%
    
    str_remove(pattern = str_match(string = ., pattern = paste0(' ' , fasta))
    ) %>% 
    
    str_replace(pattern = fasta, replacement = paste0(' ', fasta)
    ) %>% 
    
    str_replace('--kmer', ' --kmer'
    ) %>% 
    
    paste('2>&1 | tee -a', paste0(log_file_dir, 'indexing.log')
    ) %>% 
    
    str_remove('^--path ')
    
  # Write index building commands to file. 
  indexing_script_location <- getwd() %>% 
    paste0(., '/scripts/Kallisto_build_index.sh')
  
  indexing_script <- indexing_script_location %>% 
    file()
  
  writeLines(build_index, indexing_script)
  
  close(indexing_script)

  FASTQ <- organize_kallisto_arguments %>% 
    filter(
      str_detect(key, 'fastq')
    ) %>% 
    select(value) %>% 
    unlist(use.names = FALSE)
  
  FASTQ_files <- tibble(forward = list.files(FASTQ, full.names = TRUE, pattern = '*_1'), 
                        reverse = list.files(FASTQ, full.names = TRUE, pattern = '*_2')
  ) %>% 
    mutate(
      sample = forward %>% 
        basename %>% 
        str_remove('_.*')
    )
  
  arguments_to_quantify_reads <- organize_kallisto_arguments %>% 
    filter(quant == TRUE) %>% 
    select(key, value)
  
  arg_core <- arguments_to_quantify_reads %>% 
    filter(
      str_detect(key, pattern = 'path|fastq|index|output')
    )
  
  extraneous_args <- anti_join(arguments_to_quantify_reads, arg_core, by = c('key', 'value')) %>% 
    spread(key, value)
  
  extraneous_args_type <- extraneous_args %>% 
    colnames()
  
  extraneous_args_arg <- extraneous_args %>% 
    unlist(use.names = FALSE)
  
  extraneous_args <- rbind(extraneous_args_type, extraneous_args_arg) %>% 
    str_c(collapse = ' ') %>% 
    str_replace_all('= ', '=')
  
  arg_core <- arg_core %>% 
    spread(key, value)
  
  arg_core_type <- arg_core %>% 
    colnames()
  
  arg_core_arg <- arg_core %>% 
    unlist(use.names = FALSE)
  
  quant_paired_end_reads <- rbind(arg_core_type, arg_core_arg) %>% str_c(collapse = ' ') %>% str_replace_all('= ', '=') %>% 
    str_replace('^', replacement = str_match(string = ., pattern = '--path.*') %>% paste(' ', sep = '')) %>% 
    str_remove(' --path.*') %>% 
    str_replace(pattern = '$', replacement = paste(' ', extraneous_args)) %>% 
    str_remove('--path ') %>% 
    str_remove('--fastq-files') %>% 
    str_remove(FASTQ) %>% 
    str_replace('$', paste(' ', FASTQ)) %>% 
    str_replace_all(pattern = '[:space:]{2,}', replacement = ' ') %>% 
    str_replace('--index', 'quant --index') %>% 
    paste('2>&1 | tee ', log_file_dir %>% paste0('quant.log'))
    
  read_tib_to_save <- FASTQ_files %>%
    mutate(
      quantify_pe = quant_paired_end_reads,
      forard_and_reverse = map2(forward,
                                reverse,
                                paste
      ),
      quantify = pmap_chr(
        list(forward, reverse, sample, quantify_pe),
        ~str_replace(string = ..4, pattern = FASTQ, replacement = paste(..1, ..2, sep = ' ')
        ) %>%
          str_replace(string = ., pattern = output_dir, replacement = paste0(output_dir, 'Kallisto_quantifications/', ..3)) %>% 
          str_replace(string = ., pattern = '$', replacement = ' &')
      )
    ) %>%
    select(quantify) %>% 
    # Place she-bang at top of file. 
    add_row(quantify = '#!/bin/bash', .before = .1) %>% 
    unlist(use.names = FALSE)

  quant_script_location <- output_dir %>% 
    str_replace('data.*', 'scripts/') %>% 
    paste0(., 'Kallisto_quantify.sh')
  
  quant_script <- quant_script_location %>% 
    file()
  
  writeLines(read_tib_to_save, quant_script)
  
  close(quant_script)
  
  if (execute_script){
    system(indexing_script_location)
    system(quant_script_location)
  }
  
}





# Function name: implement_STAR
# Purpose: Write bash scripts for STAR index building and mapping,
#          then execute them. 
# Input: Arguments to Kallisto from config file, and all directories
#        created. 
implement_STAR <- function(STAR_arguments, directories, execute_script){
  
  index_arguments <- 'fasta|runThread|genomeDir|genomeFastaFiles|sjdbGTFfile|sjdbOverhang'
  
  organize_STAR_arguments <- STAR_arguments %>% 
    mutate(
      index = map_lgl(key, 
                      str_detect,
                      index_arguments
      ),
      mapping = map_lgl(key, 
                        str_detect,
                        index_arguments, 
                        negate = TRUE
      )
    )
  
  path <- organize_STAR_arguments %>% 
    filter(
      str_detect(key, 'path')
    ) %>% 
    select(value) %>% 
    unlist(use.names = FALSE)
  
  genome_dir <- organize_STAR_arguments %>% 
    filter(
      str_detect(key, 'genomeDir')
    ) %>% 
    select(key, value) %>% 
    unlist(use.names = FALSE) %>% 
    str_c(collapse = ' ')
  
  output_dir <- directories %>% 
    str_subset('STAR_log_files') %>% 
    str_replace(pattern = 'data.*', 'scripts/')
  
  runThread <- organize_STAR_arguments %>% 
    filter(
      str_detect(key, 'runThread')
    ) %>% 
    select(key, value) %>% 
    unlist(use.names = FALSE) %>% 
    str_c(collapse = ' ')
  
  reads <- STAR_arguments %>% 
    filter(
      str_detect(key, 'readFilesIn')
    ) %>% 
    select(value) %>% 
    unlist(use.names = FALSE) %>% 
    list.files(full.names = T) 
  
  # Check if fastq files are zipped.
  reads_zipped <- reads %>% 
    first() %>% 
    str_detect('.gz$')
    
  read_tib <- tibble(
    forward = reads %>% str_subset(pattern = '_1'),
    reverse = reads %>% str_subset(pattern = '_2'),
    sample = map_chr(forward, 
                     ~basename(.x) %>% 
                       str_remove('_.*')
    )
  )
  
  arguments_to_build_index <- organize_STAR_arguments %>% 
    filter(index == TRUE) %>% 
    select(key, value) %>% 
    spread(key, value)
  
  index_argument_type <- arguments_to_build_index %>% 
    colnames()

  index_arguments <- arguments_to_build_index %>% 
    unlist(use.names = FALSE)
  
  build_index_with <- rbind(index_argument_type, index_arguments) %>% 
    str_c(collapse = ' ') %>% 
    paste(path, '--runMode genomeGenerate', ., '2>&1 | tee -a ', directories %>% str_subset('STAR_log_files/') %>% paste0('STAR_index.log')) %>% 
    str_replace(pattern = '$', replacement = ' &') %>% 
    str_replace_all(pattern = '[:space:]{2,}', replacement = ' ')
  
  # Write bash script to build the index to a file in the scripts directory. 
  build_index_with <- paste('#!/bin/bash', build_index_with, sep = '\n')
  
  indexing_script <- paste0(output_dir, 'STAR_build_index.sh') 
  
  indexing_script_file <- indexing_script %>% 
    file()
  
  writeLines(build_index_with, indexing_script_file)
  
  close(indexing_script_file)
  
  FASTQ <- organize_STAR_arguments %>% 
    filter(
      str_detect(key, 'readFilesIn')
    ) %>% 
    select(value) %>% 
    unlist(use.names = FALSE)
  
  FASTQ_files <- FASTQ %>% 
    list.files(full.names = TRUE) %>% 
    str_c(collapse = ' ')
  
  arguments_to_map_reads <- organize_STAR_arguments %>% 
    filter(mapping == TRUE) %>% 
    select(key, value) %>% 
    spread(key, value)

  map_argument_type <- arguments_to_map_reads %>% 
    colnames()
  
  map_arguments <- arguments_to_map_reads %>% 
    unlist(use.names = FALSE)
  
  out_file_prefix <- directories %>% 
    str_subset('genome_Dir')
  
  map_reads_with <- rbind(map_argument_type, map_arguments) %>% 
    str_c(collapse = ' ') %>% 
    str_remove('--path ') %>% 
    paste(genome_dir) %>% 
    paste(runThread) 
   
  # Insert unzipping command if FASTQ files are zipped.
  if(reads_zipped){
    map_reads_with <- map_reads_with %>% 
      paste('--readFilesCommand zcat')
  }
  
  map_reads_with <- map_reads_with %>% 
    paste('--outFileNamePrefix', out_file_prefix) %>% 
    paste(
      paste(' 2>&1 | tee -a ', directories %>% 
              str_subset('STAR_log_files/') %>% 
              paste0('STAR_mapping.log')
      )
    ) %>% 
    str_replace_all(pattern = '[:space:]{2,}', replacement = ' ') %>% 
    str_replace('$', ' &')

  read_tib_to_save <- read_tib %>% 
    mutate(
      # Base command.
      map_reads_with = map_reads_with,
      
      # Combine forward with reverse read.
      forward_and_reverse = map2_chr(forward,
                                     reverse,
                                     paste
      ),
      
      # Write path to sample genome_dir for accessing sample-specific mapped reads. 
      sample_genome_dir = map_chr(sample, 
                                  ~paste0(genome_dir %>% str_remove('^.* '), '/', .x) %>% 
                                    str_replace('//', '/')
      ),
      
      # Create sample genome directories for mapped reads. 
      write_sample_genome_dir = map(sample_genome_dir, 
                                    dir.create,
                                    recursive = TRUE, 
                                    showWarnings = FALSE
      ),
      
      # Concatenate and string replace for writing full commands for mapping with STAR.
      map_reads = pmap_chr(
        list(map_reads_with, forward_and_reverse, sample_genome_dir),
        ~str_replace(string = ..1, pattern = FASTQ, replacement = ..2) %>% 
          str_replace_all(string = ., pattern = genome_dir %>% str_remove('^.* '), replacement = ..3 %>% paste0(., '/')) %>% 
          str_replace(string = ., pattern = ..3, replacement = genome_dir %>% str_remove('^.* ') %>% paste0(., '/')) %>% 
          str_replace_all(string = ., pattern = '//', replacement = '/')
      )
    ) %>% 
    select(map_reads) %>% 
    # Place shebang at top of file. 
    add_row(map_reads = '#!/bin/bash', .before = 1) %>% 
    unlist(use.names = FALSE)
    
  mapping_script_location <- paste0(output_dir %>% str_replace('data.*', 'scripts/'), 'STAR_map_reads.sh')
  
  mapping_script <- mapping_script_location %>% 
    file()
  
  writeLines(read_tib_to_save, mapping_script)
  
  close(mapping_script)
  
  if (execute_script){
    system(indexing_script)
    system(mapping_script_location)
  }
  
  
  
  return(out_file_prefix)
}





# Function name: implement_feature_counts
# Purpose: Write bash script for feature counts to execute.
# Input: path to feature counts output, gtf, genome directory containing sam files,
#        and path to the feature counts program. 
implement_feature_counts <- function(output_dir, annotation, genome_Dir, feature_counts, execute_script){
  
  feature_counts_command <- paste(feature_counts, '-a', annotation, '-o x') %>%  
    str_replace('x$', paste0(output_dir, 'clean_counts.txt'))
  
  # Format feature counts arguments using samples from STAR output.
  feature_counts_commands <- genome_Dir %>% 
    list.files(recursive = T, full.names = TRUE) %>% 
    tibble(sam_output = .) %>% 
    mutate(
      
      sample = map_chr(sam_output,
                       ~dirname(.x) %>% 
                         basename()
      ),

      feature_counts = feature_counts_command,
      
      fc = map2_chr(feature_counts,
                    sample,
                    ~str_replace(string = .x, 
                                 pattern = 'clean_counts.txt', 
                                 replacement = paste(.y, 'clean_counts.tsv', sep = '_')
                    )
      ),
      
      fc_output = map2_chr(fc, 
                           sam_output, 
                           ~paste(.x, .y)
      )
    ) %>% 
    select(fc_output) %>% 
    add_row(fc_output = '#!/bin/bash', .before = 1) %>% 
    unlist(use.names = FALSE)
  
  feature_counts_script_location <- genome_Dir %>% 
    str_replace('data.*', 'scripts/') %>% 
    paste0(., 'feature_counts.sh')
  
  feature_counts_script <- feature_counts_script_location %>% 
    file()
  
  writeLines(feature_counts_commands, feature_counts_script)
  
  close(feature_counts_script)
  
  if (execute_script){
    system(feature_counts_script_location)
  }
}
