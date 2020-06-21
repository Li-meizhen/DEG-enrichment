# DEG-enrichment
Pipeline for DEG enrichment analysis of non-model animals

### Package preparation

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

### Input files
In the current work directory, you need to make a new directory named "Input_files", which contains transcript quantification data from RSEM, annotation of OGS from pannzer ("GO.txt") and KAAS ("koID.txt" - gene_id to ko_id, "kegg_desc.txt" - ko_id to kegg description). 

Name the quantificaiton file using the following format: "commonFeatures_sample1_sample2.txt". Then we use a list to store these information. This is the only thing you need to change if you have different samples. 

