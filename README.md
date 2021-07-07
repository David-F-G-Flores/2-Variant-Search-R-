# Variant-Search-R
This code follows on from the MGI-search-R- repository. 
The final ``` mousegenes ``` object contains genes and supporting evidence.

Load bioMart and pull mouse dataset
```R
library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
```
