#Author: Colin Tang
#
#Scripts to perform differential expression analysis with htseq count data
#wd must be folder with files ending with _counts.txt 
#
#based on htseq files containing gene names and gene IDs. will delete gene IDs. Will compare 
#columns/sample names containing "_untreated" against "_treated"
#will output a differential expression file and a GSEA prerank .rnk file for GSEA.

library("DESeq2")
library("apeglm")

diffexpress<-function(countsfile) {
  rnaData <- read.table(countsfile, sep = "\t", header = T)
  geneKey <- rnaData[,1:2]
  #reorganizes the data frame. takes gene names as column1 and data only
  untreatCols <- grep("_untreated", colnames(rnaData))
  treatCols <- grep("_treated", colnames(rnaData))
  rnaData <- subset(rnaData, select = c(1, untreatCols, treatCols))
  rnaData <- as.matrix(as.data.frame(rnaData[,-1], row.names = as.character(rnaData[,1])))
  conditions <- data.frame(sample = colnames(rnaData), 
                           condition = c(rep("untreated", length(untreatCols)), rep("treated", length(treatCols))))
  
  dds <- DESeqDataSetFromMatrix(countData = rnaData, colData = conditions, design = ~ condition)
  
  #prefiltering removes low count features
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  #factor leveling
  dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
  #differential expression calculation
  dds <- DESeq(dds)
  res <- results(dds)
  res <- results(dds, name="condition_treated_vs_untreated")
  res <- results(dds, contrast=c("condition","treated","untreated"))
  #LFC shrinkage with apeglm
  resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
  
  resOrdered <- res[order(res$pvalue),]
  summary(res)
  
  resOrdered <- data.frame(GeneID = rownames(resOrdered), GeneSymbol = geneKey[match(rownames(resOrdered), geneKey[,1]),2], resOrdered)
  write.table(resOrdered, "LFC.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(subset(resOrdered, select = c(2,4)), "LFC_GSEA.rnk", sep = "\t", col.names = F, row.names = F, quote = F)
  
  
  #plots the shrunken RNA seq data with the apeglm method published here (removes noise from low count transcripts leading to high LFC):
  #Nature Methods, 13:7. 10.1038/nmeth.3885
  pdf("RNAseqMAplot.pdf")
  plotMA(resLFC, ylim=c(-3,3))
  dev.off()
}

#used for SF268 GCN2KO rna-seq data

diffexpress.2<-function(countsfile) {
  rnaData <- read.table(countsfile, sep = "\t", header = T)
  geneKey <- rnaData[,1:2]
  #reorganizes the data frame. takes gene names as column1 and data only
  untreatedCols <- grep("N675_", colnames(rnaData))
  treatedCols <- grep("Tm_", colnames(rnaData))
  rownames(rnaData) <- rnaData[,1]
  rnaData <- subset(rnaData, select = c(2, untreatedCols, treatedCols))
  rnaData <- as.matrix(as.data.frame(rnaData[,-1]))
  rnaData <- round(rnaData)
  conditions <- data.frame(sample = colnames(rnaData), 
                           condition = c(rep("untreated", length(untreatedCols)), rep("treated", length(treatedCols))))
  
  dds <- DESeqDataSetFromMatrix(countData = rnaData, colData = conditions, design = ~ condition)
  
  #prefiltering removes low count features
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  #factor leveling
  dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
  #differential expression calculation
  dds <- DESeq(dds)
  res <- results(dds)
  res <- results(dds, name="condition_treated_vs_untreated")
  res <- results(dds, contrast=c("condition","treated","untreated"))
  #LFC shrinkage with apeglm
  resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
  
  resOrdered <- res[order(res$pvalue),]
  summary(res)
  
  resOrdered <- data.frame(GeneID = rownames(resOrdered), GeneSymbol = geneKey[match(rownames(resOrdered), geneKey[,1]),2], resOrdered)
  write.table(resOrdered, "LFC.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(subset(resOrdered, select = c(2,4)), "LFC_GSEA.rnk", sep = "\t", col.names = F, row.names = F, quote = F)
  
  #plots the shrunken RNA seq data with the apeglm method published here (removes noise from low count transcripts leading to high LFC):
  #Nature Methods, 13:7. 10.1038/nmeth.3885
  pdf("RNAseqMAplot.pdf")
  plotMA(resLFC, ylim=c(-3,3))
  dev.off()
}




#uses a tab delimited conditions file to make all the LFC files. (Sample, Condition, OutputFileName)
#requires the salmon counts file to contain all the samples in the conditions file.
#sample names are the same as the filenames without the .sf extension.
#LFC analyses are separated by unique outputfilenames.

diffexpress.3<-function(countsfile, conditionsfile) {
  rnaData <- read.table(countsfile, sep = "\t", header = T)
  #Checks for duplicated genes in the counts file and deletes the second entry if found.
  duplicatedGenes <- which(duplicated(rnaData[,1]) == T)
  if (length(duplicatedGenes) > 0) {
    print("Duplicate ENSGs. Deleted entries:")
    print(rnaData[duplicatedGenes,1:2])
    rnaData <- rnaData[-duplicatedGenes,]
  }
  rownames(rnaData) <- rnaData[,1]   #makes rownames the unique ensembl gene ID
  geneKey <- rnaData[,1:2]
  conditionsdf <- read.table(conditionsfile, sep = "\t", header = T)
  
  #determines number of analyses for deseq loop
  listAnalysis <- levels(conditionsdf$OutputFileName)
  
  for (i in listAnalysis) {
    #gets RNAdata columns for samples in the ith analysis
    treatedSamples <- conditionsdf$Sample[which(conditionsdf$OutputFileName == i & conditionsdf$Condition == "treated")]
    treatedCols <- match(treatedSamples, colnames(rnaData))
    untreatedSamples <- conditionsdf$Sample[which(conditionsdf$OutputFileName == i & conditionsdf$Condition == "untreated")]
    untreatedCols <- match(untreatedSamples, colnames(rnaData))
    #creates a RNA data subset df for each analysis
    tempRNAdata <- as.matrix(subset(rnaData, select = c(untreatedCols, treatedCols)))
    conditions <- conditionsdf[which(conditionsdf$OutputFileName == i),1:2]
    colnames(conditions) <- c("sample", "condition")
    
    dds <- DESeqDataSetFromMatrix(countData = tempRNAdata, colData = conditions, design = ~ condition)
    
    #prefiltering removes low count features of < 10 reads for a row of a feature
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    #factor leveling
    dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
    #differential expression calculation
    dds <- DESeq(dds)
    res <- results(dds)
    res <- results(dds, contrast=c("condition","treated","untreated"))
    resOrdered <- res[order(res$pvalue),]
    summary(res)
    resOrdered <- data.frame(GeneID = rownames(resOrdered), GeneSymbol = geneKey[match(rownames(resOrdered), geneKey[,1]),2], resOrdered)
    
    #LFC shrinkage with apeglm. apeglm method published here (removes noise from low count transcripts leading to high LFC):
    #Nature Methods, 13:7. 10.1038/nmeth.3885
    #shrinkLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
    
    #writes LFC and gsea files
    write.table(resOrdered, paste(i, ".txt", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(subset(resOrdered, select = c(2,4)), paste(i, "_GSEA.rnk", sep = ""), sep = "\t", col.names = F, row.names = F, quote = F)
    
    #makes MA plot
    #pdf("RNAseqMAplot.pdf")
    #plotMA(resLFC, ylim=c(-3,3))
    #dev.off()
  }   #closes loop of deseq calls
}

