## DegPlotWide == modified function
library(dplyr)

myDegPlotWide <- function (counts, genes, group, metadata = NULL, batch = NULL) {
  if (is.null(metadata)) 
    metadata = data.frame(colData(counts))
  metadata = data.frame(metadata)
  if (class(counts)[1] == "DESeqDataSet") {
    dd = bind_rows(lapply(genes, function(gene) {
      plotCounts(counts, gene, intgroup = group, returnData = TRUE) %>% 
        mutate(count = log2(count + 1)) %>% mutate(gene = gene, 
                                                   sample = row.names(metadata))
    }))
    dd$group = dd[[group]]
  }
  else if (class(counts)[1] == "matrix") {
    dd = melt(counts[genes, ])
    colnames(dd) = c("gene", "sample", "count")
    dd$group = as.factor(metadata[as.character(dd$sample), 
                                  group])
  }
  else {
    stop("No supported for class", class(counts)[1])
  }
  if (is.null(group)) {
    dd$treatment = "one_group"
  }
  else {
    dd$treatment = dd[, "group"]
  }
  p = ggplot(dd, aes_string(x = "treatment", y = "count", color = "treatment"))
  if (!is.null(batch)) {
    dd$batch = as.factor(metadata[dd$sample, batch])
    p = ggplot(dd, aes(x = "treatment", y = "count", shape = "batch"))
  }
  p = p + geom_point(position = position_jitterdodge(dodge.width = 0.9), alpha=.5) + 
    stat_summary(fun = "median", geom = "crossbar") + 
    facet_wrap(~gene, scales="fixed", nrow=1) + 
    xlab("") + ylab("Normalized Counts") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45L, hjust = 1L))
  p
}

