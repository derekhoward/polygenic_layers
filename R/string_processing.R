library(stringr)
library(homologene)

process_input_genes <- function(input_genes) {
  processed_genes <- unlist(str_split(trimws(input_genes), "[, \n]+") )
  return(processed_genes)
}

convert_genes <- function(input_genes) {
  human_genes <- mouse2human(input_genes)
  return(unique(human_genes$humanGene))
}


#genes <- c('CADM2', 'ZNF704', 'NCAM1', 'RABEP2', 'ATP2A1')

#genes <- c('CYP3A5','CYP3A4','CYP3A7','CYP39A1','CYP3A43')

#convert_genes(genes)
