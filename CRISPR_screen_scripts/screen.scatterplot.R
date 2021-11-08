#script to make scatter plot of a data frame with two samples
#will scatter plot against them and add labels to the top 10
#sorted by column three
#df is formatted as gene, sample1, sample2

library("ggplot2")
library("wesanderson")
library("ggrepel")

screen.scatterplot<-function(datatable, ngenes = 10, top = TRUE) {
  datatable <- datatable[order(datatable[,3], decreasing = top),]
  topn <- datatable[ngenes,3]
  ggdata <- ggplot(datatable, aes(x = datatable[,2], y = datatable[,3]))
  if (top == TRUE && ngenes > 0) {
    labeltext <- geom_text(aes(label=ifelse(datatable[,3]>=topn,as.character(gene),'')),hjust=0,vjust=0)
  } else if (top == FALSE && ngenes > 0) {
    labeltext <- geom_text(aes(label=ifelse(datatable[,3]<=topn,as.character(gene),'')),hjust=0,vjust=0)
  }
  graph <- ggdata + geom_point(size = 1, color= "grey") + 
    labeltext +
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.5) +
    ylab(paste(colnames(datatable)[3], "sgRNA Frequency", sep = "\n")) + xlab(paste(colnames(datatable)[2], "sgRNA Frequency", sep = "\n")) + 
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.7, linetype = "solid")) +
    ylim(c(min(datatable[,3]), max(datatable[,3])))
  print(graph)
  ggsave("sgRNA_scatterplot.pdf", width = 10, height = 7)
}


# updated version uses ggrepel and uses ngenes is used to define number of top genes if top is TRUE, 
# otherwise it can accept a vector containing a list of specific gene names desired to label. 
# Black labels is a new argument for controlling the label color.
screen.scatterplot2<-function(datatable, ngenes = 10, top = TRUE, blackLabels = FALSE) {
  datatable <- datatable[order(datatable[,3], decreasing = top),]
  datatable$Effect <- datatable[,3] - datatable[,2]
  datatable <- na.omit(datatable)
  if (top == TRUE && ngenes > 0) {
    labeltext <- rep("", length(datatable[,1]))
    labeltext[1:ngenes] <- datatable[1:ngenes,1]
  } else if (top == FALSE && ngenes > 0) {   #enters if wanting a custom label of genes
    labeltext <- rep("", length(datatable[,1]))
    genesofInterest <- match(ngenes, datatable[,1])
    labeltext[genesofInterest] <- datatable[genesofInterest,1]
  }
  if (blackLabels == FALSE && ngenes > 0) {
    ggdata <- ggplot(datatable, aes(x = datatable[,2], y = datatable[,3], color=Effect, label = labeltext)) + 
      geom_point(size = 1.5) + geom_text_repel()
  } else if (blackLabels == TRUE && ngenes > 0) {
    ggdata <- ggplot(datatable, aes(x = datatable[,2], y = datatable[,3], color=Effect, label = labeltext)) + 
      geom_point(size = 1.5) + geom_text_repel(color = "black", box.padding = 0.3, force = 0.1, min.segment.length = 0.1)
  } else {
    ggdata <- ggplot(datatable, aes(x = datatable[,2], y = datatable[,3], color=Effect)) + geom_point(size = 1.5)
  }
  graph <- ggdata + 
    scale_color_gradientn(colors = c("steelblue","grey","red"), 
                          values = scales::rescale(c( min(datatable$Effect), 0, max(datatable$Effect)))) + 
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.5) +
    ylab(paste(colnames(datatable)[3], "sgRNA Frequency", sep = " ")) + xlab(paste(colnames(datatable)[2], "sgRNA Frequency", sep = " ")) + 
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.7, linetype = "solid"), 
          axis.text = element_text(size=14), axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold")) +
    #coord_cartesian(xlim = c(-4,2.5), ylim = c(-1.25,7.25))   #for manual lims
    ylim(c(min(datatable[,3]), max(datatable[,3])))   #for automatic lims
  print(graph)
  ggsave("sgRNA_scatterplot2.pdf", width = 8, height = 5)
}
