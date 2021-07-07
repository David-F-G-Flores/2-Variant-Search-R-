# Variant-Search-R
This code follows on from the MGI-search-R- repository. 
The final ``` mousegenes ``` object contains genes and supporting evidence.

Because multiple studies (supporting evidence) on the same gene, the ``` mousegenes ``` object held duplicated genes. The object was sorted by publication and duplicate genes removed, using R.

Open bioMart library.
```R
library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
```
Extract genes from object, Mouse Genome Informatics (MGI) nomenclature.
```R
mousegenes<-mousegenes$OntologyAnnotation.subject.primaryIdentifier
```

A useful bioMart piece of code below. 
```R
mousegenes<-mousegenes$OntologyAnnotation.subject.primaryIdentifier
```
