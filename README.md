# DEG-enrichment
Pipeline for DEG enrichment analysis of non-model animals

## Package preparation

Packages required in this project are mostly from Bioconductor packages.

```{r, eval = FALSE}

library(DESeq2)
library(affy)
library(GOstats)
library(KEGG.db)
library(GSEABase)
library(DOSE)
library(RColorBrewer)
library(clusterProfiler)

```

## Input files
In the current work directory, you need to make a new directory named "Input_files", which contains transcript quantification data from RSEM formatted as follows. Name the quantificaiton file using the following format: "commonFeatures_sample1_sample2.txt". Then we use a list to store these information. This is the only thing you need to change if you have different samples. 
Also make sure that annotation files of OGS from pannzer ("GO.txt") and KAAS ("koID.txt" - gene_id to ko_id, "kegg_desc.txt" - ko_id to kegg description) are formated as follows. 

### adult_female_male.txt
```
gene_id adult_female1 adult_female2 adult_female3 adult_male1 adult_male2 adult_male3
PSOL00001-PA  75.50 108.50  70.00 43.50 37.00 37.00
PSOL00002-PA  9.12  85.67 10.36 15.38 15.95 49.45
```

### GO.txt
```
qpid  ontology  goid  desc
PSOL00001-PA  MF  GO:0004831  tyrosine-tRNA ligase activity
PSOL00001-PA  BP  GO:0006437  tyrosyl-tRNA aminoacylation
```

### koID.txt
```
PSOL00001-PA  K01866
PSOL00002-PA  K04950
```

### kegg_desc.txt
```
K00844   hexokinase
K12407   glucokinase
```


