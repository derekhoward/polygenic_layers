library(stringr)
library(homologene)

process_input_genes <- function(input_genes) {
  processed_genes <- unlist(str_split(trimws(input_genes), "[, \n]+") )
  return(processed_genes)
}

convert2human <- function(input_genes, in_species) {
  if (in_species == 'Mouse') {
    human_genes <- mouse2human(input_genes)
    return(unique(human_genes$humanGene))
  } else if (in_species == 'Rhesus Macaque') {
    human_genes <- homologene(input_genes, inTax = 9544, outTax = 9606)
    return(unique(human_genes$humanGene))
  } else {
    return(unique(input_genes))
  }
  #return(unique(human_genes$humanGene))
}

convert2nhp <- function(input_genes, in_species) {
  if (in_species == 'Mouse') {
    nhp_genes <- homologene(input_genes, inTax = 10090, outTax = 9544)
    return(unique(nhp_genes$`9544`))
  } else if (in_species == 'Human') {
    nhp_genes <- homologene(input_genes, inTax = 9606, outTax = 9544)
    return(unique(nhp_genes$`9544`))
  } else {
    return(unique(input_genes))
  }
  #return(unique(nhp_genes$`9544`))
}

#genes <- c('CADM2', 'ZNF704', 'NCAM1', 'RABEP2', 'ATP2A1')

#genes <- c('CYP3A5','CYP3A4','CYP3A7','CYP39A1','CYP3A43')

#convert_genes(genes)
