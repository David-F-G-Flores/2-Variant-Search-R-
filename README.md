# Variant-Search-R
The code below was used to search for the *Bos taurus* homologs, of mouse genes associated with prenatal lethality.
This code follows on from the repository "MGI-search-R-", and used the final ``` mousegenes ``` object which contained annotated mouse genes.

Because multiple studies (supporting evidence) of the same gene were identified, the ``` mousegenes ``` object held duplicate genes records. The object was sorted by the year of publication to retain the most recent, and duplicate genes removed. R code not shown for this sorting operation.

This method used the biomaRt library to find homologs, the chromosomal positions of the genes, variants within and associated SIFT score information. Specific genome assemblys were requested in these biomaRt requests.

Open bioMart library.
```R
library(biomaRt)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
```
Extracted genes from ``` mousegenes ``` object, Mouse Genome Informatics (MGI) nomenclature.
```R
mousegenes2<-mousegenes$OntologyAnnotation.subject.primaryIdentifier
```

Useful bioMart code below to convert Nomenclature from different sources.
This code initiates a Bos taurus list, and queries bioMart for ensembl *Bos taurus* homolog, using MGI ID.  
Populating this list took some time and was saved into a csv as a security measure because of a lost connection, code not shown here.
```R
btauLr<-list()

for (genes in mousegenes2){
  print(genes)
  data = getBM( attributes=c("external_gene_name","ensembl_gene_id","btaurus_homolog_ensembl_gene",
                             "btaurus_homolog_associated_gene_name","description","btaurus_homolog_orthology_confidence"), filters= ("mgi_id"), values =genes,mart=mouse)
  btauLr[[genes]] <- data
}
btauLr = do.call('rbind',btauLr)
btauLr$OntologyAnnotation.subject.primaryIdentifier <- row.names(btauLr) ## rownames to column for a later merge 
```
Some mice genes in the query above did not return homologs using MGI ID, however, various other gene names were identified e.g. "MGI:104982" returned *CEBPG* and "ENSMUSG00000056216" which is mouse ensembl nomenclature, but no *Bos taurus* ensembl number. Furthermore, numerous ensembl numbers were retrieved e.g. *ANK2* returned two *Bos taurus* ensembl numbers, "ENSBTAG00000054394" and "ENSBTAG00000002392" however validating these through the ensembl browser, this issue appeared to be related to changing of genome versions. All homologs were extracted and used to query the desired *Bos taurus* genome. If working with a limited number of genes, these alternative names may be helpful to rerun the analysis with some modifications to the attributes and filters. If in doubt the ensembl website https://www.ensembl.org/index.html should be referred to. Also note the attribute ```"btaurus_homolog_orthology_confidence"``` which may help with filtering and identifying best genes if working with a large number.



*Note* *Bost taurus* UMD3.1 Ensembl release 94 was updated to ARS-UCD1.2 during this project. Next command uses archived database.
Using the correct bioMart assemblies, biomaRt btau and btauSNP marts were loaded. Gene start/end positions was searched in btau, and used to search for SNP within using btauSNP.

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
The code below looped over btau ensembl IDs, and returned the position of the gene. This was followed by a query for snp within that location. These variants are further filtered on sift score consequence, and stored into ```TopSNPs```. This can take some time, consider saving ```TopSNPs``` into an external file on each iteration. Furthermore, the consequence type may be filtered on the actual sift score instead of the consequence as is carried out below but should be aware not all variants are associated with a SIFT score and is blank.
```R
variants<-list()
for (gene in btauLr$btaurus_homolog_ensembl_gene) {
  print(gene)
  data=getBM(attributes=c("start_position","end_position","chromosome_name"),filters="ensembl_gene_id",values=gene, mart=cowinfo)
  chr=c(data[,3])
  start=data[,1]
  end=data[,2]
  if(length(start)!=0){
    print('Got position')
    SNPS<-getBM(c("ensembl_gene_stable_id","refsnp_id","chr_name","chrom_strand","allele","chrom_start","ensembl_type",
                  "consequence_type_tv","sift_prediction","sift_score","distance_to_transcript"), 
                filters=c("start", "end","chr_name"),values=list(start,end,chr), mart=snpmart)
    print('Got SNP')
    TopSNPs=SNPS[SNPS$consequence_type_tv %in% c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", 
                                                 "stop_gained", "frameshift_variant", "stop_lost", "start_lost", 
                                                 "transcript_amplification", "inframe_insertion", "inframe_deletion",
                                                 "missense_variant", "protein_altering_variant"),]
    variants[[gene]]<- TopSNPs
    print('Filtered on consequence')
    }
}
variants<-do.call('rbind',variants)
```
Merge back on the suppporting evidence. Second merge of supporting evidence with variants.
```R
variants2 <- merge(btauLr,mousegenes, by="OntologyAnnotation.subject.primaryIdentifier")
variantsMaster<- merge(variants2,variants,by.x="btaurus_homolog_ensembl_gene",by.y ="ensembl_gene_stable_id")
```
