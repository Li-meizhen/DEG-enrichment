#BiocManager::install("DESeq2")
#BiocManager::install("affy")
#BiocManager::install("GOstats")
#BiocManager::install("KEGG.db")
#BiocManager::install("GSEABase")
#BiocManager::install("DOSE")
#install.packages("RColorBrewer")

# adult_female_male.txt
# gene_id adult_female1 adult_female2 adult_female3 adult_male1 adult_male2 adult_male3
# PSOL00001-PA  75.50 108.50  70.00 43.50 37.00 37.00
# PSOL00002-PA  9.12  85.67 10.36 15.38 15.95 49.45

# GO.txt
# qpid  ontology  goid  desc
# PSOL00001-PA  MF  GO:0004831  tyrosine-tRNA ligase activity
# PSOL00001-PA  BP  GO:0006437  tyrosyl-tRNA aminoacylation

# koID.txt
# PSOL00001-PA  K01866
# PSOL00002-PA  K04950

# kegg_desc.txt
# K00844   hexokinase
# K12407   glucokinase

library(DESeq2)
library(affy)
library(GOstats)
library(KEGG.db)
library(GSEABase)
library(DOSE)
library(RColorBrewer)
library(clusterProfiler)
mycolor <- brewer.pal(12,"Paired")

###########################################
# give input files, output files and list variables to store all analysis results 
###########################################

inputfile <- c("adult_female_male.txt", "larva_female_male.txt", "female_larva_adult.txt", "male_larva_adult.txt")

sample <- list(length = length(inputfile))  # store sample information in list "sample"
for (i in 1:length(inputfile) ){
    title <- gsub("(\\w+)_(\\w+)_(\\w+).txt",replacement = "\\1", inputfile[i]) 
    sample1 <- gsub("(\\w+)_(\\w+)_(\\w+).txt",replacement = "\\2", inputfile[i])   
    sample2 <- gsub("(\\w+)_(\\w+)_(\\w+).txt",replacement = "\\3", inputfile[i])   
    sample[[i]] <- c(title, sample1, sample2)   # store sample names for future references
}

# some output files
resdata <- list(length = length(inputfile))   # DESeq2 counts
outputfile <- gsub(pattern = ".txt", replacement = "_exp.csv", inputfile)   # DESeq2 output files
volcanoPlot <- gsub(pattern = ".txt", replacement = "_volcano.png", inputfile)
MAplotBefore <- gsub(pattern = ".txt", replacement = "_MA_before_normalization.png", inputfile)
MAplotAfter <- gsub(pattern = ".txt", replacement = "_MA_after_normalization.png", inputfile)

# Result directories
dir.create("Results")
dir.create("Results/exp")
dir.create("Results/MA_plots")
dir.create("Results/Volcano_plots")
dir.create("Results/GOstats_KEGG")
dir.create("Results/GOstats_GO")
dir.create("Results/ClusterProfiler_KEGG")
dir.create("Results/ClusterProfiler_GO")

##########################################
# customized MA plot function
##########################################
myMAplot <- function(dataframe, filename, ...) {
    meanFrame <- as.data.frame( cbind(rowMeans(dataframe[, 1:3]), rowMeans(dataframe[, 4:6])) )
    meanFrame <- meanFrame + 1
    png(file = filename, height = 1200, width = 1600)
    par(mar = c(8,8,10,3))
    ma.plot( rowMeans(log2(meanFrame)), log2(meanFrame[,1])-log2(meanFrame[,2]), cex=1.5, xlim = c(1,20), ...)
    dev.off()
    return(meanFrame)
}


###########################################
# DESeq2 analysis
###########################################

for  (i in 1:length(inputfile) ) {
    # reads count data
    countData <- read.table(paste("Input_files/",inputfile[i], sep = ""), header = T, row.names = 1)
    
    ########################################
    # MA plot before normalization
    meanframeB <- myMAplot(countData, paste("Results/MA_plots/", MAplotBefore[i], sep = ""), 
                           main = paste("MA plot of gene expression in", sample[[i]][2], "and", sample[[i]][3], sample[[i]][1], "mealybugs" ),
                           mgp = c(5,1.5,0), cex.main = 4, cex.lab = 3, cex.axis = 2)
    
    #########################################
    # DESeq2 analysis
    # conditions
    condition <- c(rep("C1",3), rep("C2",3))
    colData <- data.frame(row.names = colnames(countData), condition)

    # interger mode required
    countData <- round(countData)

    # normalization and filtered out low-count reads
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design =~ condition )
    dds2 <- DESeq(dds)
    dds2 <- dds2[ rowSums(counts(dds2))>1, ]
    
    # merge expression analysis data with normalized count
    res <- results(dds2) 
    resdata[[i]] <- merge(as.data.frame(res), as.data.frame(counts(dds2, normalize=TRUE)), by="row.names", sort=FALSE)  
    resdata[[i]] <- na.omit(resdata[[i]])
    
    # write dateframe resdata to file
    write.csv(resdata[[i]], file = paste("Results/exp/", outputfile[i], sep = ""), row.names = F)
    
    #############################################
    # MA plot after DESeq2 normalization
    deseq2frame <- resdata[[i]][, 8:13]
    meanframeB <- myMAplot(deseq2frame, paste("Results/MA_plots/", MAplotAfter[i], sep = ""), 
                           main = paste("MA plot of gene expression in", sample[[i]][2], "and", sample[[i]][3], sample[[i]][1], "mealybugs" ),
                           mgp = c(5,1.5,0), cex.main = 4, cex.lab = 3, cex.axis = 2)
    
    ###############################################
    # volcano plot
    # extract data for volcano plot
    expdata <- resdata[[i]]    # just for convenience
    
    # log2 fold change > 1 and corresponding pvalues
    log2foldChange2 <- expdata$log2FoldChange[expdata$log2FoldChange > 1]
    pvalue2 <- -log10(expdata$padj[expdata$log2FoldChange > 1])
    
    # log2 fold change < -1 and corresponding pvalues
    log2foldChange1 <- expdata$log2FoldChange[expdata$log2FoldChange < -1]
    pvalue1 <- -log10(expdata$padj[expdata$log2FoldChange < -1])
    
    # -1 <= log2 fold change <= 1 and corresponding pvalues
    log2foldChange0 <- expdata$log2FoldChange[expdata$log2FoldChange <= 1 & expdata$log2FoldChange > -1]
    pvalue0 <- -log10(expdata$padj[expdata$log2FoldChange <= 1 & expdata$log2FoldChange > -1])
   
    # plot
    png(file = paste("Results/Volcano_plots/", volcanoPlot[i], sep = ""), width = 1200, height = 1000)
    
    par(mar = c(8,8,10,3))
    plot(x = log2foldChange0, y = pvalue0, type = "p", pch = 16, cex = 1, col = "gray",
         main = paste("Gene expression in", sample[[i]][1], "mealybugs"), ylab = "-log10 pvalue",
         xlab = paste("log2 fold change (", sample[[i]][3], "/", sample[[i]][2], ")"), mgp = c(5,1.5,0),
         cex.main = 4, cex.lab = 3, cex.axis = 2, ylim = c(0, 80), xlim = c(-20,20) )
    
    lines(x = log2foldChange1, y = pvalue1, type = "p", pch = 16, cex = 1, col = mycolor[8])
    lines(x = log2foldChange2, y = pvalue2, type = "p", pch = 16, cex = 1, col = mycolor[2])
    
    # add horizontal line of pvalue = 0.01
    abline(h = 2, col = "red", lty = 5)    
    dev.off()
    
}

#######################################
# GOstats for enrichment analysis
#######################################

#########################################
# GO ontology
# universal settings

GOprediction <- read.table("Input_files/GO.txt", header = T, sep = "\t", quote = "")
types <- c("BP", "MF", "CC")

# do enrichment analysis according to molecular function, biological process, cellular component respectively
for (type in types) {
    oriGeneSet <- GOprediction[GOprediction$ontology == type, ]
    oriGeneSet$frame.Evidence <- rep("ISS", nrow(oriGeneSet))         # add gene evidence code column
    geneSet <- data.frame(oriGeneSet$goid, oriGeneSet$frame.Evidence, oriGeneSet$qpid)   # change column order
    colnames(geneSet) <- c("frame.go_id", "frame.Evidence", "frame.gene_id")
    goFrame <- GOFrame(geneSet, organism = "psol")
    goAllFrame <- GOAllFrame(goFrame)
    gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
    universe <- as.character( unique(geneSet$frame.gene_id) )
    
    # enrichment analysis on each sample couple
    for (i in 1:length(inputfile)) {
        genelist1 <- resdata[[i]]$Row.names[ resdata[[i]]$log2FoldChange < -1 & resdata[[i]]$padj < 0.05 ]
        genelist2 <- resdata[[i]]$Row.names[ resdata[[i]]$log2FoldChange > 1 & resdata[[i]]$padj < 0.05 ]
        gene1 <- intersect(genelist1, universe)    # some genes didn't have annotations
        gene2 <- intersect(genelist2, universe)
        
        # set parameters for sample1
        params1 <- GSEAGOHyperGParams(name = "my params 1",
                             geneSetCollection = gsc,
                             geneIds = gene1,
                             universeGeneIds = universe,
                             ontology = type,
                             pvalueCutoff = 0.05,
                             conditional = FALSE,
                             testDirection = "over")

        Over1 <- hyperGTest(params1)
        htmlReport(Over1, file = paste("Results/GOstats_GO/", sample[[i]][1], "_", sample[[i]][2], "_", type, "_enrichment.html", sep = ""))
    
        # set parameters for sample1
        params2 <- GSEAGOHyperGParams(name = "my params 2",
                             geneSetCollection = gsc,
                             geneIds = gene2,
                             universeGeneIds = universe,
                             ontology = type,
                             pvalueCutoff = 0.05,
                             conditional = FALSE,
                             testDirection = "over")

        Over2 <- hyperGTest(params2)
        htmlReport(Over2, file = paste("Results/GOstats_GO/", sample[[i]][1], "_", sample[[i]][3], "_", type, "_enrichment.html", sep = ""))
    }

}


######################################
# KEGG pathway
# universal settings
kegg <- read.table("Input_files/koID.txt", header = F)
keggframeData <- data.frame(kegg$V2, kegg$V1)
colnames(keggframeData) <- c("frame.path_id", "frame.gene_id")
keggframeData$frame.path_id <- gsub("K", replacement="", keggframeData$frame.path_id)  # remove "K" from koID

keggFrame <- KEGGFrame(keggframeData, organism = "psol")
gscK <- GeneSetCollection(keggFrame, setType = KEGGCollection())
universeK <- as.character( unique(keggframeData$frame.gene_id) )

for (i in 1:length(inputfile)) {
    genelist1 <- resdata[[i]]$Row.names[ resdata[[i]]$log2FoldChange < -1 & resdata[[i]]$padj < 0.05 ]
    genelist2 <- resdata[[i]]$Row.names[ resdata[[i]]$log2FoldChange > 1 & resdata[[i]]$padj < 0.05 ]
    gene1 <- intersect(genelist1, universeK)    # some genes didn't have annotations
    gene2 <- intersect(genelist2, universeK)
  
    # set parameters for sample1
    kparams1 <- GSEAKEGGHyperGParams(name = "my paramsKEGG 1",
                                     geneSetCollection = gscK,
                                     geneIds = gene1,
                                     universeGeneIds = universeK,
                                     pvalueCutoff = 0.05,
                                     testDirection = "over")
  
    kOver1 <- hyperGTest(kparams1)
    htmlReport(kOver1, file = paste("Results/GOstats_KEGG/", sample[[i]][1], "_", sample[[i]][2], "_", "KEGG_enrichment.html", sep = ""))
  
    # set parameters for sample2
    kparams2 <- GSEAKEGGHyperGParams(name = "my paramsKEGG 2",
                                     geneSetCollection = gscK,
                                     geneIds = gene2,
                                     universeGeneIds = universeK,
                                     pvalueCutoff = 0.05,
                                     testDirection = "over")
  
    kOver2 <- hyperGTest(kparams2)
    htmlReport(kOver2, file = paste("Results/GOstats_KEGG/", sample[[i]][1], "_", sample[[i]][3], "_", "KEGG_enrichment.html", sep = ""))
  
}


#########################################################
# use ClusterProfiler to generate barplots and dotplots
#########################################################

######################################
# GO enrichment
for (type in types) {
    
    oriGeneSet <- GOprediction[GOprediction$ontology == type, ]
    universe <- as.character( unique(oriGeneSet$qpid) )
    term2gene <- data.frame(oriGeneSet$goid, oriGeneSet$qpid)
    term2name <- data.frame(oriGeneSet$goid, oriGeneSet$desc)
    
    for (i in 1:length(inputfile)) {
        
        genelist1 <- resdata[[i]]$Row.names[ resdata[[i]]$log2FoldChange < -1 & resdata[[i]]$padj < 0.05 ]
        genelist2 <- resdata[[i]]$Row.names[ resdata[[i]]$log2FoldChange > 1 & resdata[[i]]$padj < 0.05 ]
        gene1 <- intersect(genelist1, universe)    # some genes didn't have annotations
        gene2 <- intersect(genelist2, universe)

        x1 <- enricher(gene1, TERM2GENE=term2gene, TERM2NAME=term2name)
        x2 <- enricher(gene2, TERM2GENE=term2gene, TERM2NAME=term2name)
        
        
        # barplot sample1
        png(file = paste("Results/ClusterProfiler_GO/", sample[[i]][1], "_", sample[[i]][2], "_", type, "_barplot.png", sep = ""), 
            width = 1200, height = 800)
        print( barplot(x1, drop=TRUE, showCategory=30) )
        dev.off()
        
        # barplot sample2
        png(file = paste("Results/ClusterProfiler_GO/", sample[[i]][1], "_", sample[[i]][3], "_", type, "_barplot.png", sep = ""), 
            width = 1200, height = 800)
        print( barplot(x2, drop=TRUE, showCategory=30) )
        dev.off()
        
        # dotplot sample1
        png(file = paste("Results/ClusterProfiler_GO/", sample[[i]][1], "_", sample[[i]][2], "_", type, "_dotplot.png", sep = ""), 
            width = 1200, height = 800)
        print( dotplot(x1, showCategory=30) )
        dev.off()
        
        # dotplot sample2
        png(file = paste("Results/ClusterProfiler_GO/", sample[[i]][1], "_", sample[[i]][3], "_", type, "_dotplot.png", sep = ""), 
            width = 1200, height = 800)
        print( dotplot(x2, showCategory=30) )
        dev.off()
    }   

}

####################################
# KEGG enrichment

term2gene <- keggframeData
term2name <- read.table("Input_files/kegg_desc.txt", header = F, sep = "\t", quote = "")
term2name$V1 <- gsub("K", replacement = "", term2name$V1)
universeK <- as.character(term2gene$frame.gene_id)

for (i in 1:length(inputfile)) {
    i <- 1
    genelist1 <- resdata[[i]]$Row.names[ resdata[[i]]$log2FoldChange < -1 & resdata[[i]]$padj < 0.05 ]
    genelist2 <- resdata[[i]]$Row.names[ resdata[[i]]$log2FoldChange > 1 & resdata[[i]]$padj < 0.05 ]
    gene1 <- intersect(genelist1, universeK)    # some genes didn't have annotations
    gene2 <- intersect(genelist2, universeK)
  
    x1 <- enricher(gene1, TERM2GENE=term2gene, TERM2NAME=term2name)
    x2 <- enricher(gene2, TERM2GENE=term2gene, TERM2NAME=term2name)
  
  # barplot sample1
    png(file = paste("Results/ClusterProfiler_KEGG/", sample[[i]][1], "_", sample[[i]][2], "_barplot.png", sep = ""), 
        width = 1200, height = 800)
    print( barplot(x1, drop=TRUE, showCategory=30) )
    dev.off()
  
    # barplot sample2
    png(file = paste("Results/ClusterProfiler_KEGG/", sample[[i]][1], "_", sample[[i]][3], "_barplot.png", sep = ""), 
        width = 1200, height = 800)
    print( barplot(x2, drop=TRUE, showCategory=30) )
    dev.off()
  
    # dotplot sample1
    png(file = paste("Results/ClusterProfiler_KEGG/", sample[[i]][1], "_", sample[[i]][2], "_dotplot.png", sep = ""), 
        width = 1200, height = 800)
    print( dotplot(x1, showCategory=30) )
    dev.off()
  
    # dotplot sample2
    png(file = paste("Results/ClusterProfiler_KEGG/", sample[[i]][1], "_", sample[[i]][3], "_dotplot.png", sep = ""), 
        width = 1200, height = 800)
    print( dotplot(x2, showCategory=30) )
    dev.off()
}   

