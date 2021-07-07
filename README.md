# Variant-Search-R
Below was used to search for Bos taurus homologs of the previously identified mouse genes associated with perinatal lethality.
This code follows on from the MGI-search-R- repository, and uses the final ``` mousegenes ``` object which contains mouse genes and the supporting evidence.

Because multiple studies (supporting evidence) of the same gene were identified, the ``` mousegenes ``` object held duplicated genes. The object was sorted by publication and duplicate genes removed, using R.

Open bioMart library.
```R
library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
```
Extract genes from ``` mousegenes ``` object, Mouse Genome Informatics (MGI) nomenclature.
```R
mousegenes2<-mousegenes$OntologyAnnotation.subject.primaryIdentifier
```

A useful bioMart piece of code below to convert Nomenclature.
Initiate Bos taurus list, query bioMart for ensembl Bos taurus homolog using MGI ID.  
Populating this list may take some time.
```R
btauLr<-list()

for (genes in mousegenes2){
  data = getBM( attributes=c("external_gene_name","ensembl_gene_id","btaurus_homolog_ensembl_gene",
                             "btaurus_homolog_associated_gene_name","description"), filters= ("mgi_id"), values =genes,mart=mouse)
  btauLr[[genes]] <- data
}
btauLr = do.call('rbind',l)
btauLr$OntologyAnnotation.subject.primaryIdentifier <- row.names(btauLr) ## rownames to column for a later merge
```
Some genes from above did not return homologs using MGI ID however various other gene names are identified and may be rerun using the above code and changing the attributes and filters.
*Note* UMD3.1 Ensembl release 94 was updated to ARS-UCD1.2 during this project. Next command uses archived snp database.
Using correct bioMar assemblies, biomaRt btau and btauSNP marts loaded. Gene start/end positions will be searched, and used to search for SNP within.

To load btau biomaRt.
```R
cowinfo <- useEnsembl(biomart = "ensembl", 
                   dataset = "btaurus_gene_ensembl", 
                   version = "94")
```
To load btau SNP biomaRt.
```R
snpmart <- useEnsembl(biomart = "ENSEMBL_MART_SNP", 
                   dataset = "btaurus_snp", 
                   version = "94")
```
Loops over btau ensembl IDs, gets position of gene. Queries for snp within that location. Filters snp on consequence into ```TopSNPs```. This can take some time, consider writing table and saving ```TopSNPs``` into external file.
```R
variants<-list()
for (gene in btauLr$btaurus_homolog_ensembl_gene) {
  print(gene)
  data=getBM(attributes=c("start_position","end_position","chromosome_name"),filters="ensembl_gene_id",values=gene, mart=cowinfo)
  chr=c(data[,3])
  start=data[,1]
  end=data[,2]
  print('Got position')
  SNPS<-getBM(c("ensembl_gene_stable_id","refsnp_id","chr_name","chrom_strand","allele","chrom_start","ensembl_type",
                "consequence_type_tv","sift_prediction","sift_score","distance_to_transcript"), 
              filters=c("start", "end","chr_name"),values=list(start,end,chr), mart=snpmart)
  print('Got SNP')
  TopSNPs=SNPS[SNPS$consequence_type_tv %in% c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", 
                                               "stop_gained", "frameshift_variant", "stop_lost", "start_lost", 
                                               "transcript_amplification", "inframe_insertion", "inframe_deletion",
                                               "missense_variant", "protein_altering_variant"),]
  print('Filtered on consequence')
}
variants<-do.call('rbind',variants)
```
Merge back on the suppporting evidence. Merge supporting evidence with variants.
```R
variants2 <- merge(btauLr,mousegenes, by="OntologyAnnotation.subject.primaryIdentifier")
variantsMaster<- merge(variants2,variants,by.x="btaurus_homolog_ensembl_gene",by.y ="ensembl_gene_stable_id")
```
