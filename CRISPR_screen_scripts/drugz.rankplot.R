library("ggplot2")
library("ggrepel")

# Script to make labeled rank plot from the drugZ output file

drugz.rankplot<-function(dataFile, genelabels = NULL, width = 6, height = 6, outputName = "drugZ_rankplot.pdf") {
  datatable <- read.table(dataFile, sep = "\t", header = T, stringsAsFactors = F)
  datatable <- subset(datatable, select = c(1,4,6))
  datatable <- datatable[order(datatable$normZ),]
  datatable <- na.omit(datatable)
  labeltext <- rep("", length(datatable[,1]))
  if (is.null(genelabels) == FALSE) {
    generows <- data.frame(gene = genelabels, index = match(genelabels, datatable$GENE))
    generows <- na.omit(generows)
    labeltext[generows[,2]] <- as.character(generows[,1])
  }
  ggdata <- ggplot(datatable, aes(x = rank_synth, y = normZ, color=normZ))
  graph <- ggdata + geom_hline(yintercept = 0, color="black", linetype="dashed", size=0.5) +
    geom_point(size = 2) + 
    scale_color_gradientn(colors = c("steelblue","grey","red"), 
                          values = scales::rescale(c(min(datatable$normZ),-2.5, 0, 2.5,max(datatable$normZ)))) + 
    ylab("Normalized Z score") + xlab("sgRNA Rank") + theme_bw() + 
    theme(axis.text = element_text(size=18), axis.title = element_text(size=18, face="bold"), legend.position = "none") +
    #coord_cartesian(xlim = c(-4,2.5), ylim = c(-1.25,7.25))   #for manual lims
    ylim(c(min(datatable$normZ), max(datatable$normZ))) + #for automatic lims
    scale_y_continuous(breaks = seq(-5,20,5), labels = seq(-5,20,5)) + 
    geom_text_repel(label = labeltext, box.padding = 1, min.segment.length = 0.1, size = 7)
  ggsave(outputName, width = width, height = height)
}


# Debug at line 2 and run for manual gene label changes
#datatable[which(datatable[,1] == "PPP1R15A"),1] <- "GADD34"
#datatable[which(datatable[,1] == "EIF2AK4"),1] <- "GCN2"


drugz.rankplot2<-function(dataFile, hits = 3, width = 6, height = 6) {
  outputName <- paste(dataFile, "_rankplot.pdf")
  datatable <- read.table(dataFile, sep = "\t", header = T, stringsAsFactors = F)
  datatable <- subset(datatable, select = c(1,4,6))
  datatable <- datatable[order(datatable$normZ),]
  datatable <- na.omit(datatable)
  genenum <- nrow(datatable)
  labeltext <- rep("", genenum)
  if (hits > 0) {
    rowIndex <- c(seq(1,hits), seq(genenum - hits + 1, genenum))
    generows <- data.frame(gene = datatable$GENE[rowIndex], index = rowIndex)
    labeltext[generows[,2]] <- as.character(generows[,1])
  }
  ggdata <- ggplot(datatable, aes(x = rank_synth, y = normZ, color=normZ))
  graph <- ggdata + geom_hline(yintercept = 0, color="black", linetype="dashed", size=0.5) +
    geom_point(size = 2) + 
    scale_color_gradientn(colors = c("steelblue","grey","red"), 
                          values = scales::rescale(c(min(datatable$normZ),-2.5, 0, 2.5,max(datatable$normZ)))) + 
    ylab("Normalized Z score") + xlab("sgRNA Rank") + theme_bw() + 
    theme(axis.text = element_text(size=18), axis.title = element_text(size=18, face="bold"), legend.position = "none") +
    ggtitle(gsub(".txt", "",dataFile)) + 
    #coord_cartesian(xlim = c(-4,2.5), ylim = c(-1.25,7.25))   #for manual lims
    #ylim(c(min(datatable$normZ), max(datatable$normZ))) + #for automatic lims
    scale_y_continuous(breaks = seq(-5,5,2.5), labels = seq(-5,5,2.5)) + 
    geom_text_repel(label = labeltext, box.padding = 0.4, min.segment.length = 0.25, size = 5.5, max.overlaps = 20)
  ggsave(outputName, width = width, height = height)
}
