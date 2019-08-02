#!/usr/bin/Rscript
# Small script for easy visualization
#   of gene expression across a section

library(ggplot2)
library(argparse)
source("spatial_funcs.r")
source("utils.r")

# Parser -----------------------------
parser <- ArgumentParser()

parser$add_argument("-c","--count_file",
                    type = "character",
                    default = NULL
                    )

parser$add_argument("-m","--meta_file",
                    type = "character",
                    default = NULL
                    )

parser$add_argument("-o","--out_dir",
                    type = "character",
                    default =  "/tmp/geneExpression",
                    nargs = "+")

parser$add_argument("-gs","--gene_set",
                    type = "character",
                    default = NULL,
                    nargs = "+")

parser$add_argument("-an","--show_annot",
                    default = F,
                    action = 'store_true')




args <- parser$parse_args()

# Main -----------------------------

# use example data if no count or meta
# data is provided
if (is.null(args$count_file)){
    load("data/example_data.Rdata")
    ct <- ex_data$count_data
    ct <- sweep(ct, 1, rowSums(ct),'/') * 100
    mt <- ex_data$meta_data
    bname <- "ex_data"
} else {
    data <- readSTsection(args$count_file,
                          args$meta_file)
    ct <- data$count_data
    mt <- data$meta_data
    bname <- gsub("count_data\\-|\\.tsv",
                      "",
                      basename(args$count_file))
 
}

   

# Read set of genes either from file
# or spae separated strings

if (is.null(args$gene_set[1])) {
    genelist <- c("ERBB2")
} else if(file.exists(args$gene_set[1])){
    genelist <- read.table(args$gene_set,
                           sep = '\n',
                           header = 1,
                           row.names = NA)
} else {
    genelist <- args$gene_set
}

# create ouptu directory if it does not exists
checkMakeDir(args$out_dir)
# get coordinates
crd <- getCoordinates(rownames(ct))
# create long-formatted data frame to be plotted
plotdf <- data.frame(x = crd[,1],
                     y = crd[,2])
    
if (args$show_annot) {
    plotdf$annot <- mt$tumor
    clrm <- scale_colour_manual(values=alpha(c("green",
                                         "red"),0.7))
} else {
    plotdf$annot <- "noanno"
    clrm <- scale_colour_manual(values="black")
}

# plot distribution of each gene
for (gene in genelist) {
    if (!(gene %in% colnames(ct))) {
            print(sprintf(" %s is not present in the data",gene))
            next
    }

    val <- ct[,gene] 
    plotdf$val <- ct[,gene]

    g <- ggplot(plotdf, aes(x = x, y = y)) + 
            geom_point(size = 5, pch = 21, aes(fill = val,color = annot)) +
            scale_fill_gradient(low = "#FFFFFF",
                                high = "#56B1F7",
                                space = "Lab",
                                na.value = "grey50",
                                guide = "colourbar",
                                aesthetics = "fill",
                                oob = scales::squish) +
            clrm
    # set output name 
    oname <- paste(args$out_dir,
                   paste(gene,bname,"viz.png",sep = '.'),
                   sep = '/')

    png(oname)
    print(g)
    dev.off()
}
