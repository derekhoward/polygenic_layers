library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(stringr)
library(jsonlite)
library(data.tree)

# check working directory:
# should be projects/polygenic_layers
getwd()
setwd("~/projects/polygenic_layers")

DATA_FOLDER <- file.path(getwd(), './data/raw/allen_NHP/')

# get the fetal ontology as data.tree structure
onto_list <- fromJSON('./data/NHP_ontology.json', simplifyDataFrame = FALSE)
ontology <- as.Node(onto_list, mode='explicit', check='no-warn')

neocortical_structures <- FindNode(ontology, 'neocortex')
structures_of_interest <- neocortical_structures$Get('name')
# the hippocampal structures also contain some samples from zones/layers of interest
#hippocampal_structures <- FindNode(ontology, 'hippocampal cortex (hippocampal formation)')

get_zone <- function(node) {
  name <- node$name
}
# assigns a new 'zone' value to each node that has a parent_structure_id == 294021746 (which is the node for neocortex)
neocortical_structures$Do(function(node) node$zone <- get_zone(node), filterFun = function(x) x$parent_structure_id == 294021746)
# pushes the 'zone' value from a parent node to its children and returns a dataframe
zone_markers <- ToDataFrameNetwork(neocortical_structures, 'name', 'zone', inheritFromAncestors = TRUE) %>% 
  select(name, zone) %>% 
  as_tibble()

print(unique(zone_markers$zone))

zones_to_drop <- c("dorsolateral prefrontal cortex (areas 9 and 46)", "medial orbital frontal cortex (area 14)",
                   "caudal orbital frontal cortex (area 13)", "medial frontal cortex (areas 24, 25 and 32)",
                   "rostral cingulate cortex (areas 24 and 32)", "occipital cortex")

zone_markers %<>% filter(!zone %in% zones_to_drop)

## Read in data
samples <- read_csv(file = file.path(DATA_FOLDER, 'lmd_expression_matrix_2014-03-06', 'columns_metadata.csv'))
expression <- read_csv(file = file.path(DATA_FOLDER, 'lmd_expression_matrix_2014-03-06', 'expression_matrix.csv'), col_names = F)
# remove columns which have NA for all expression values
expression %<>% select_if(~sum(!is.na(.)) > 0)
names(expression) = c('probe_id', 1:(ncol(expression)-1))
probes <- read_csv(file = file.path(DATA_FOLDER, 'lmd_expression_matrix_2014-03-06', 'rows_metadata.csv'))
# remove parens() from beginning and end of gene symbols, ex: '(A1BG)' --> 'A1BG' 
probes$gene_symbol <- str_replace(probes$gene_symbol, "(^\\()", "")
probes$gene_symbol <- str_replace(probes$gene_symbol, "(\\)$)", "")
#probes %>% filter(str_detect(gene_symbol, "(^\\()"))

# merge in zone markers, this also restricts samples to those of interest from the ontology
samples <- inner_join(samples, zone_markers, by=c('structure_name' = 'name'))

# merge probe info to gene expression
samples_exp <- left_join(expression, probes, by=c('probe_id' = 'row_num'))
dim(samples_exp)

# melt samples_exp and merge in samples info
# generates a tidy tibble with 32480 expression values per sample (915 samples of interest), along with metadata (sample_id, donor_id, age, structure/zone sampled)
gene_exp_long_annotated <- samples_exp %>% 
  select(-c(gene_id, probe_name, gene_name, entrez_id)) %>% 
  # convert = TRUE is required to convert sample_id to numeric data type
  gather(key=sample_id, value=expression, -c(gene_symbol, probe_id), convert=TRUE) %>% 
  inner_join(samples, by=c('sample_id' = 'column_num')) %>% 
  select(-c(structure_id, well_id, structure_acronym)) %>% 
  filter(!is.na(gene_symbol))
dim(gene_exp_long_annotated)

# get meanExp by gene for each donor
# then, get ranking of genes in zones
# finally, get mean zscore across donors within zones
gene_exp_long_annotated %<>% 
  group_by(donor_id, zone, gene_symbol) %>% 
  summarise(meanExp = mean(expression)) %>% 
  group_by(donor_id, zone) %>% 
  mutate(ranking = rank(meanExp)) %>% 
  #group_by(zone, gene_symbol) %>% 
  group_by(gene_symbol) %>%
  mutate(zscore = (ranking - mean(ranking)) / sd(ranking))
dim(gene_exp_long_annotated)

cortical_zones_expression_matrix <- gene_exp_long_annotated %>% 
  select(-meanExp, -ranking) %>% 
  group_by(zone, gene_symbol) %>% 
  summarise(mean_zscore = mean(zscore), n_donors = n_distinct(donor_id))

cortical_zones_expression_matrix %<>% 
  select(-n_donors) %>% 
  spread(key = zone, value = mean_zscore)

result <- cortical_zones_expression_matrix %>% 
  select(-gene_symbol) %>% 
  map_df(rank)

result$gene_symbol <- cortical_zones_expression_matrix$gene_symbol

NHP_cortical_zones_expression_matrix <- cortical_zones_expression_matrix
NHP_cortical_zones_ranks <- result

write_csv(NHP_cortical_zones_expression_matrix , './data/processed/NHP_cortical_zones_expression_matrix.csv')
write_csv(NHP_cortical_zones_ranks , './data/processed/NHP_cortical_zones_ranks.csv')
save(NHP_cortical_zones_expression_matrix, file='./data/processed/NHP_cortical_zones_expression_matrix.Rdata')
save(NHP_cortical_zones_ranks, file='./data/processed/NHP_cortical_zones_ranks.Rdata')


