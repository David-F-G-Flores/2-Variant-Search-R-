# Variant-Search-R
This code follows on from the MGI-search-R- repository. 
The ```R mousegenes ```


Load bioMart and pull mouse dataset
```R
library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
```
