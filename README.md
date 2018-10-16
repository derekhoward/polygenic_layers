# polygenic_layers

## to download the Developing Human Brain Atlas ontology:
curl -o ./data/developing_onto.json http://api.brain-map.org/api/v2/structure_graph_download/16.json 

## to download the Non-Human Primate Brain Atlas ontology:
curl -o ./data/NHP_ontology.json http://api.brain-map.org/api/v2/structure_graph_download/8.json

NOTE:
requires minimal processing of .json files after: extract 'msg' field for the whole ontology

## Test genelists
Genelist file taken from:
https://medicine.yale.edu/lab/rakic/transcriptome/

from supplemental of the paper
Ayoub, Albert E., et al. "Transcriptional programs in transient embryonic zones of the cerebral cortex defined by high-resolution mRNA sequencing." Proceedings of the National Academy of Sciences (2011): 201112213.

### Gene list processing
Run
> python process_Rakic_genes.py
to get 6 gene lists from the .xls file
