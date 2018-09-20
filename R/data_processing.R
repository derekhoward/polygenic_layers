library(tidyverse)
library(data.tree)

# check working directory:
# should be projects/polygenic_layers/R
getwd()
#
setwd("~/projects/polygenic_layers/")

# set up directories/raw data
# currently only setup for a single donor
DATA_FOLDER <- file.path(getwd(), 'data/raw/allen_human_fetal_brain')
donors <- c('lmd_matrix_12566', 'lmd_matrix_12690', 'lmd_matrix_12840', 'lmd_matrix_14751')

# get the fetal ontology as data.tree structure
library('jsonlite')
onto_list <- fromJSON('../molecular_AN/data/fetal21ontology.json', simplifyDataFrame = FALSE)
ontology <- as.Node(onto_list, mode='explicit', check='no-warn')

# get list of structures we care about
transient_structures <- FindNode(ontology, 'transient structures of forebrain')
structures_of_interest <- transient_structures$Get('name')
# there are 652 structure names which we will filter for
length(structures_of_interest)

# create a datafrane with the structures of interest and how they map to specific layer markers (marker's are closer to root of ontology)
get_zone <- function(node) {
  name <- node$name
}
transient_structures$Do(function(node) node$zone <- get_zone(node), filterFun = function(x) x$parent_structure_id == 10506)
zone_markers <- ToDataFrameNetwork(transient_structures, 'name', 'zone', inheritFromAncestors = TRUE) %>% 
  select(name, zone)

# process each donor separately
all_processed_sample_dfs <- list()
for (donor in donors) {
  #DONOR_FOLDER <- file.path(path=DATA_FOLDER, 'lmd_matrix_12566')#donor)
  DONOR_FOLDER <- file.path(path=DATA_FOLDER, donor)
  samples_file <- list.files(path=DONOR_FOLDER, pattern="columns_metadata.csv", recursive = TRUE, full.names = TRUE)
  expression_file <- list.files(path=DONOR_FOLDER, pattern="expression_matrix.csv", recursive = TRUE, full.names = TRUE)
  probes_file <- list.files(path=DONOR_FOLDER, pattern="rows_metadata.csv", recursive = TRUE, full.names = TRUE)
  
  ## load data
  samples <- read_csv(file = samples_file)
  # add a sample_id col to facilitate joining with exp
  samples$sample_id <- c(1:nrow(samples))
  # merge in zone markers, this also restricts samples to those of interest from "transient_structures" ontology
  samples <- left_join(samples, zone_markers, by=c('structure_name' = 'name'))
  exp <- read_csv(file = expression_file, col_names = FALSE)
  names(exp) = c('probe_id', 1:(ncol(exp)-1))
  probes <- read_csv(file = probes_file)
  
  # merge probe info to gene expression
  annotated_sample_exp <- left_join(exp, probes, by=c('probe_id' = 'probeset_id'))
  dim(annotated_sample_exp)
  
  # melt annotated_sample_exp and merge in samples info 
  gene_exp_long_annotated <- annotated_sample_exp %>% 
      select(-c(probeset_name, gene_id, gene_name, entrez_id, chromosome)) %>% 
      # convert = TRUE is required to convert sample_id to numeric data type
      gather(key=sample_id, value=expression, -c(gene_symbol, probe_id), convert=TRUE) %>% 
      left_join(samples, by='sample_id') %>% 
      select(-c(structure_id, well_id, structure_acronym)) %>% 
      filter(!is.na(zone))

  # should probably drop genes/columns that start with A_ here 
  
  # get the mean of gene expression values for each zone/layer
  zone_exp <- gene_exp_long_annotated %>% 
    group_by(zone, gene_symbol) %>% 
    summarise(meanExp = mean(expression))
  
  zone_exp %<>% 
    group_by(zone) %>% 
    mutate(ranking = rank(meanExp))
  
  zone_exp %<>% 
    group_by(gene_symbol) %>% 
    mutate(zscore = (ranking - mean(ranking)) / sd(ranking))
  
  zone_exp$donor_id <- donor
  
  all_processed_sample_dfs[[donor]] <- zone_exp
}

all_donors <- bind_rows(all_processed_sample_dfs)

cortical_zones_expression_matrix <- all_donors %>% 
  select(-meanExp, -ranking) %>% 
  group_by(zone, gene_symbol) %>% 
  summarise(mean_zscore = mean(zscore), n_donors=n_distinct(donor_id))

cortical_zones_expression_matrix %<>% 
            select(-n_donors) %>% 
            spread(key = zone, value = mean_zscore)

result <- cortical_zones_expression_matrix %>% 
  select(-gene_symbol) %>% 
  #map_df(desc) %>% 
  map_df(rank)

result$gene_symbol <- cortical_zones_expression_matrix$gene_symbol

# somehow the following creates twice as many rows as it should?   
#spread(cortical_zones_expression_matrix, key = zone, value = mean_zscore) %>% 
#  select(-n_donors)

#write_csv(result, './data/processed_data_final.csv')
save(result, file='./processed_developmental_zones.Rdata')
