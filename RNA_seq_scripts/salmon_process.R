#script to process the salmon quant.sf transcript files into gene abundance working directory should
#be directory with all .sf files

library("tximport")
library("biomaRt")

salmon.process<-function() {
  dir <- system.file("extdata", package = "tximportData")
  tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v27.csv"))
  
  quantFiles <- list.files(pattern = ".sf")
  txi <- tximport(quantFiles, type = "salmon", tx2gene = tx2gene)

  sampleNames <- gsub(".sf", "", quantFiles)
  abundancedf <- data.frame(rownames(txi$abundance), txi$abundance)
  countsdf <- data.frame(rownames(txi$counts), txi$counts)
  colnames(abundancedf) <- c("gene_id", sampleNames)
  colnames(countsdf) <- c("gene_id", sampleNames)
  
  abundancedf$ensembl_id <- gsub('\\..+$', '', abundancedf$gene_id)
  countsdf$ensembl_id <- gsub('\\..+$', '', countsdf$gene_id)
  #gets symbols from biomart
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  symbol <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), filters = "ensembl_gene_id", values = abundancedf$ensembl_id, mart = mart)
  
  abundancedf <- merge(x = symbol, y = abundancedf, by.x="ensembl_gene_id", by.y="ensembl_id")
  countsdf <- merge(x = symbol, y = countsdf, by.x="ensembl_gene_id", by.y="ensembl_id")
  #removes duplicates
  #abundancedf <- abundancedf[-which(duplicated(abundancedf[,3]) == TRUE),]
  
  #removes empty gene names, and ensembl id column
  abundancedf <- abundancedf[-which(abundancedf[,2] == ""),]
  countsdf <- countsdf[-which(countsdf[,2] == ""),]
  abundancedf$ensembl_gene_id <- NULL
  countsdf$ensembl_gene_id <- NULL
  abundancedf <- subset(abundancedf, select = c(2,1,3:length(abundancedf)))
  countsdf <- subset(countsdf, select = c(2,1,3:length(countsdf)))
  countsdf[,3:length(countsdf)] <- round(countsdf[,3:length(countsdf)])   #rounds counts to nearest integer
  
  colnames(abundancedf)[1:2] <- c("GeneID", "GeneSymbol")
  colnames(countsdf)[1:2] <- c("GeneID", "GeneSymbol")
  write.table(abundancedf, "salmongene_abundance.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  write.table(countsdf, "salmongene_counts.txt", sep = "\t", col.names = T, row.names = F, quote = F)
}
