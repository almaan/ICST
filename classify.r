#!/usr/bin/Rscript

library("argparse")
library("iC10")
library("ggplot2")
library("futile.logger")
library("umap")


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

print(args)

# Main ----------------------

MNS <- 4
MXS <- 10

cpth <- args$count_file
mpth <- args$meta_file
odir <- args$output_dir

# sort to match meta and count data
cpth <- sort(cpth)
mpth <- sort(mpth)

# to store loaded data
clist <- list()
mlist <- list()
plist <- list()
slist <- list()
tlist <- list()

index <- c(0)

for (num in 1:length(cpth)) {
    flog.info(sprintf("Analyzing sample %s | %d/%d",cpth[num],num,length(cpth)))
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
    flog.info('Initiating iC10 classification... ')
    res <- tryCatch(classify(ct), error = function(e) F)
    
    # if classification is successfull save results    
    if (!is.logical(res)) {

      flog.info("Classification Successfull")
      
      # grep for sample id
      grp <- regexpr("[0-9]{5}_[A-Z][0-9]",basename(cpth[num]))
      idx <- substr(basename(cpth[1]), grp[1], grp[1] + attr(grp,'match.length'))
      
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
      plist[[num]] <- posterior
      tlist[[num]] <- idx
      index <- c(index,dim(mt)[1])

      # open image file and write
      png(file = paste(c(odir,gsub("\\.tsv","\\.png",basename(cpth[num]))),collapse = "/"),
          width = 720,
          height = 720)
      
      # scale according to spot library size
      libsize <- apply(ct,2,sum)
      mxlib <- max(libsize)
      scale <- MNS + (libsize / mxlib) * (MXS - MNS)
      
      # store scaling factors
      slist[[num]] <-scale
      
      # visualize classification
      g <- ggplot(mt, aes(x = xcoord, y = ycoord)) +
        geom_point(shape = 21, aes(fill = class, color = tumor), size = scale) +
        scale_colour_manual(values=c("tumor"="black","non"="white")) +
        ggtitle(paste(c(idx,'highest posterior iC10'),collapse = ' '))
      print(g)
      
      #close image file
      dev.off()
      
    }  else {
      flog.error("Unsuccessful classification")
  }
}

# compute indices for each replicate in joint matrix
index <- cumsum(index)

# create joint matrix
joint_posterior <- data.frame(matrix(0,index[length(index)],dim(posterior)[2]))
for (pos in 1:(length(index)-1)) {
  joint_posterior[(index[pos]+1):(index[pos+1]),] <- plist[[pos]] 
}

# compress with umap
cnfg <- umap.defaults
cnfg$n_components <- 3

clr <- umap(joint_posterior, config = cnfg)
clr <- clr$layout
# normalize 
mx <- apply(clr, 2, max)
mn <- apply(clr,2,min)
clr <- sweep(clr,2,mn,'-')
clr <- sweep(clr,2,(mx-mn),'/')

# split joint matrix and plot compressed visualization
for (pos in 1:(length(index)-1)) {
    srgb <- clr[(index[pos]+1):index[pos+1],]
    srgb = rgb(r = srgb[,1],g = srgb[,2],b = srgb[,3])
    
    
    png(file = paste(c(odir,gsub("\\.tsv","\\.umap\\.png",basename(mpth[pos]))),collapse = "/"),
        width = 720,
        height = 720)
    
    g <- ggplot(mlist[[pos]], aes(x = xcoord, y = ycoord)) +
      geom_point(shape = 21, fill = srgb, aes(color = tumor), size = slist[[pos]]) +
      scale_colour_manual(values=c("tumor"="black","non"="white")) + 
      ggtitle(paste(c(tlist[[pos]],'UMAP compression of posterior'),collapse = ' '))
    print(g)
    dev.off()
}