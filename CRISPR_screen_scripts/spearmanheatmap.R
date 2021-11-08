#script to make correlation heatmap
#input is in format gene then samples
#performs spearman rank correlation
#

library("ggplot2")
  
screen.corheatmap<-function(datatable) {
  datatable$gene <- NULL
  numsamples <- length(datatable)
  cortable <- data.frame(colnames(datatable))
  #makes correlation table
  for (i in 1:numsamples) {
    for(j in 1:numsamples) {
      cortable[(j),(i+1)] <- cor.test(datatable[,i],datatable[,j], method = "spearman")$estimate
    }
  }
  colnames(cortable) <- c("Sample", colnames(datatable))
  cortable.reordered <- melt(cortable)
  colnames(cortable.reordered) <- c("SampleX", "SampleY", "Rank_correlation")
  ggdata <- ggplot(cortable.reordered, aes(x = SampleX, y = SampleY))
  graph <- ggdata + geom_tile(aes(fill = Rank_correlation), color ="black") + 
    geom_text(aes(label=round(Rank_correlation, digits = 2))) + 
    scale_fill_gradientn(colors = c("steelblue","white") ,
                         values = scales::rescale(c(0.5,1))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  ggsave("corrleation_heatmap.pdf", width = 7, height = 5)
}