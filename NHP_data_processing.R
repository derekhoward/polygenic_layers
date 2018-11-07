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

#zones_to_drop2 <- c("layer I", "layer II", "layer II/III", "layer III", "layer IV", "layer V", "layer VI" )

zone_markers %<>% filter(!zone %in% zones_to_drop)
#zone_markers %<>% filter(!zone %in% zones_to_drop2)


## Read in data
samples <- read_csv(file = file.path(DATA_FOLDER, 'lmd_expression_matrix_2014-03-06', 'columns_metadata.csv'))
expression <- read_csv(file = file.path(DATA_FOLDER, 'lmd_expression_matrix_2014-03-06', 'expression_matrix.csv'), col_names = F)
colnames(expression) = c('probe_id', 1:(ncol(expression)-1))

#only 129 line up
sum(expression$probe_id == 1:(nrow(expression)))
#solution: don't use that first column
expression$probe_id <- 1:(nrow(expression))

probes <- read_csv(file = file.path(DATA_FOLDER, 'lmd_expression_matrix_2014-03-06', 'rows_metadata.csv'))

probes %<>% filter(!is.na(gene_symbol))
# remove parens() from beginning and end of gene symbols, ex: '(A1BG)' --> 'A1BG' 
#probes$gene_symbol <- str_replace(probes$gene_symbol, "(^\\()", "")
#probes$gene_symbol <- str_replace(probes$gene_symbol, "(\\)$)", "")
#remove paren probes
probes %<>% filter(!str_detect(gene_symbol, "(^\\()"))

# merge in zone markers, this also restricts samples to those of interest from the ontology
samples <- inner_join(samples, zone_markers, by=c('structure_name' = 'name'))

#only use first to agepoints
samples %<>% filter(age %in% c("E40", "E50"))


# merge probe info to gene expression
samples_exp <- left_join(probes, expression, by=c('row_num' = 'probe_id')) %>% rename(probe_id = row_num)
dim(samples_exp)

# remove columns which have NA for all expression values
#samples_exp %>% select_if(~sum(is.na(.)) > 0)
samples_exp %<>% select_if(~sum(!is.na(.)) > 0)
#as.data.frame(samples_exp[1,])

# melt samples_exp and merge in samples info
# generates a tidy tibble with 32480 expression values per sample (915 samples of interest), along with metadata (sample_id, donor_id, age, structure/zone sampled)
gene_exp_long_annotated <- samples_exp %>% 
  select(-c(gene_id, probe_name, gene_name, entrez_id)) %>% 
  # convert = TRUE is required to convert sample_id to numeric data type
  gather(key=sample_id, value=expression, -gene_symbol, -probe_id, convert=TRUE) %>% 
  inner_join(samples, by=c('sample_id' = 'column_num')) %>% 
  select(-c(structure_id, well_id, structure_acronym)) 

dim(gene_exp_long_annotated)
gene_exp_long_annotated %>% filter(gene_symbol == "DAAM1", structure_name == "ventricular zone of V1") #spot check
# get meanExp by gene for each donor
# then, get ranking of genes in zones
# finally, get mean zscore across donors 
gene_exp_long_annotated %<>% 
  group_by(donor_id, zone, gene_symbol) %>% 
  summarise(meanExp = mean(expression, na.rm = T)) %>% 
  group_by(donor_id, zone) %>% 
  mutate(ranking = rank(meanExp)) %>% 
  group_by(gene_symbol) %>%
  mutate(zscore = (ranking - mean(ranking, na.rm = T)) / sd(ranking, na.rm = T))
dim(gene_exp_long_annotated)

#gene_exp_long_annotated %<>% filter(!is.na(meanExp))

cortical_zones_expression_matrix <- gene_exp_long_annotated %>% 
  select(-meanExp, -ranking) %>% 
  group_by(zone, gene_symbol) %>% 
  summarise(mean_zscore = mean(zscore, na.rm = T), n_donors = n_distinct(donor_id))

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

