#!/usr/bin/Rscript

source('funcs.r')
source('utils.r')

library(ggplot2)
library(gridExtra)

preload <- T

if (preload) {

    load('data/all_samples.R')

} else {
    
    cpths <- c('COUNT_DIRS')
    mpths <- c('META_DIRS')

    cpths <- sort(cpths)
    mpths <- sort(mpths)

    dta <- load_multiple(cpths,mpths)

} 

gsets <- read.table("/tmp/autocorr/20190726064916.autocorr.clustering.tsv",
                    sep = '\t',
                    row.names = 1,
                    stringsAsFactors = F
                    )
print(head(gsets))
unisection <- unique(dta$meta_data$unique_id)
uniclust <- unique(gsets$cluster)
nsections <- length(unisection)

pdf('/tmp/allviz.pdf',
    width = 13 * 7,
    height = 3 * 7)

for (cl in uniclust) {
    gset <- rownames(gsets)[gsets$cluster == cl]
    par(mfrow=c(13,3))
    gl <- list()
    for (section in 1:nsections) {
    
        idx <- which(dta$meta_data$unique_id == unisection[section])
        escore <- getEnrichementScore(dta$count_data[idx,],
                                      gset)
    
        plotdf <- data.frame(x = as.numeric(dta$meta_data[idx,c('xcoord')]),
                             y = as.numeric(dta$meta_data[idx,c('ycoord')]),
                             annot = dta$meta_data[idx,c('tumor')],
                             val = escore
                             )
    
        gl[[section]] <- ggplot(plotdf,
                                 aes(x = x,
                                     y = y)
                                 ) + 
                          geom_point(pch = 21,
                                     size = 5,
                                     aes(fill = val,
                                         color = annot
                                         )
                                     ) + 

                          scale_fill_gradient(low = "#FFFFFF",
                                            high = "#56B1F7",
                                            space = "Lab",
                                            na.value = "grey50",
                                            guide = "colourbar",
                                            aesthetics = "fill",
                                            oob = scales::squish) +
             
                          scale_colour_manual(values=alpha(c("green",
                                                             "red"),
                                                           0.3)) 

    }

    print(sprintf('Visualized cluster %d',cl))
    grid.arrange(grobs = gl,ncol = 13)
    grid::grid.newpage()
}

dev.off()

