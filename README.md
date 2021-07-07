# Variant-Search-R
This code follows on from the MGI-search-R- repository. 
The final ``` mousegenes ``` object contains genes and supporting evidence.

The ``` mousegenes ``` object held duplicated genes because of multiple studies. The object was sorted by publication and duplicate genes removed.
```R
library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
```
