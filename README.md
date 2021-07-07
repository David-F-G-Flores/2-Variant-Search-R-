# Variant-Search-R
This code follows on from the MGI-search-R- repository. 
The final ``` mousegenes ``` object contains genes and supporting evidence.

Because multiple studies (supporting evidence) on the same gene, the ``` mousegenes ``` object held duplicated genes. The object was sorted by publication and duplicate genes removed.

Open bioMart library.
```R
library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
```
