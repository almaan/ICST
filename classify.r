#!/usr/bin/Rscript

library("argparse")
library("iC10")
library("ggplot2")


# Functions ---------------

classify <- function(x) {
  features <- matchFeatures(Exp=ct, Exp.by.feat="gene")
  features <- normalizeFeatures(features, "scale")
  res <- iC10(features)
  
  return(res)
} 

# Parser ----------------------

parser <- ArgumentParser()
parser$add_argument("-c","--count_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-m","--meta_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-o","--output_dir",
                    type = "character",
                    default = "/tmp/ic10")

args <- parser$parse_args()

# Main ----------------------


cpth <- args$count_file
mpth <- args$meta_file
odir <- args$output_dir

# sort to match meta and count data
cpth <- sort(cpth)
mpth <- sort(mpth)

# to store loaded data
clist <- list()
mlist <- list()


for (num in 1:length(cpth)) {
    print(sprintf("Analyzing sample %s | %d/%d",cpth[num],num,length(cpth)))
    # load count data
    ct <- read.table(cpth[num],
                     sep = '\t',
                     header = 1,
                     row.names = 1)
    #load meta data
    mt <- read.table(mpth[num],
                     sep = '\t',
                     header = 1,
                     row.names = 1)
    
    # select only spots found in meta and count data
    inter <- intersect(rownames(ct),rownames(mt))
    ct <- t(ct[inter,])
    mt <- mt[inter,]
    
    # Classify spots
    print('Initiating iC10 classification... ')
    res <- tryCatch(classify(ct), error = function(e) F)
    
    # if classification is successfull save results    
    if (!is.logical(res)) {
      print("Classification Successfull")
      
      posterior <- res$posterior
      class <- res$class
      # add class and posterior information to metadata
      mt <- cbind(mt,class)
      mt <- cbind(mt,posterior)
      # save new meta data
      write.csv(mt, file = paste(c(odir,basename(mpth[num])),collapse = "/"))
      
      # append to lists (for additonal analysis)
      clist[[num]] <- ct
      mlist[[num]] <- mt
  
      # open image file and write
      png(file = paste(c(odir,gsub("\\.tsv","\\.png",basename(cpth[num]))),collapse = "/"),
          width = 720,
          height = 720)
      
      # visualize classification
      g <- ggplot(mt, aes(x = xcoord, y = ycoord)) +
            geom_point(shape = 21, aes(fill = class, color = tumor), size = 10) +
            scale_colour_manual(values=c("tumor"="black","non"="white"))
      print(g)
      #close image file
      dev.off()
      
    }  else {
      print("Unsuccessful classification")
  }
}
