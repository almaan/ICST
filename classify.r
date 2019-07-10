#!/usr/bin/Rscript

sh <- suppressPackageStartupMessages

sh(library("argparse"))
sh(library("iC10"))
sh(library("ggplot2"))
sh(library("futile.logger"))
sh(library("umap"))
sh(library("RANN"))
sh(library("igraph"))
sh(library("base"))
sh(library("RColorBrewer"))

# Functions ---------------

l2_dist <- function(x) {
  # Similarity measure based on l2 norm
  # 
  # Args:
  #   x: 2xnfeatures matrix containing 
  #       where rows represent the vectors
  #       for which similarity is to be
  #       determined
  # Returns:
  #   similarity l2-based similariy of vectors
  return(1.0 / dist(x) ^ 2)
}

cos_dist <- function(x) {
  # Similarity measure based on angle between vectors
  # 
  # Args:
  #   x: 2xnfeatures matrix containing 
  #       where rows represent the vectors
  #       for which similarity is to be
  #       determined
  # Returns:
  #   cosine based similariy of vectors
  
  nx <- x[1,] / sqrt(x[1,] %*% x[1,])
  ny <- x[2,] / sqrt(x[2,] %*% x[2,])
  return( 1.0 + (nx %*% ny))
}


pool_by_membership <- function(cnt, labels) {
  # Pool spots with same labels
  #
  # Args:
  #   cnt: n_spots x n_genes - count matrix
  #   labels: n_spots - vector of labels
  # Returns:
  #   n_labels x n_genes - pooled count matrix
  
  n_labels <- length(unique(labels))
  n_genes <- dim(cnt)[2]
  pooled <- as.data.frame(matrix(0,nrow = n_labels,ncol = n_genes))
  colnames(pooled) <- colnames(cnt)
  for (label in unique(labels)) {
    idx <- which(labels == label)
    if (length(idx) > 1) {
      pooled[label,] <- colMeans(cnt[idx,])
    } else {
      pooled[label,] <- cnt[idx,]
    }
  }
  
  return(pooled)
}


classify <- function(x,
                     feat = 'gene',
                     method = 'scale') {
  # classify samples using iC10
  # 
  # Args:
  #   x: n_genes x n_samples - count matrix
  #   feat: string - type of gene identifiers (gene for HGNC)
  #   method: string - type of normalizing method.
  #                     see normalizeFeatures for more
  # Returns:
  #   iC10 results object
  
  features <- matchFeatures(Exp=x, Exp.by.feat=feat)
  features <- normalizeFeatures(features, method)
  res <- iC10(features)
  
  return(res)
} 


louvain_clustering <- function(xx,yy,datum,nn = 8) {
  # Assemble graph and perform louvain clustering
  # based on spatial and gene-expression information
  # Args:
  #   xx: n_samples x 1 - vector of x-coordinates
  #   yy: n_samples x 1 - vector of y-coordinates
  #   datum: n_samples x n_genes - matrix of expression data 
  #   nn: integer - number of neighbours to use
  #
  # Returns:
  #   list with graph, layout, labels, number of communities
  
  eps <- 0.01
  radius <- c('4' = (1 + eps),
              '8' = (sqrt(2) + eps))
  
  crd <- data.frame(x= xx, y = yy)
  nbr <- nn2(crd,crd,
             k = nn + 1,
             searchtype = 'radius',
             radius = radius[as.character(nn)])[[1]]
  
  tov <- c()
  frv <- c()
  weights <- c()
  
  for ( from in 1:length(xx)) {
    for (to in 1:length(xx)) {
      if ((to %in% nbr[from,]) &  !(to == from)) {
        verts <- sort(c(to,from))
        tov <- c(tov,verts[1])
        frv <- c(frv,verts[2])
        delta <- dist_fun(datum[c(verts[1],verts[2]),])
        weights <- c(weights,delta)
      }
    }
  }
  
  tedges <- data.frame(from = tov, to = frv)
  rownames(tedges) <- paste('x',c(1:dim(tedges)[1]))
  edges <- unique(tedges)
  weights <- weights[which(rownames(edges) %in% rownames(tedges))] 
  
  graph <- graph_from_data_frame(edges,
                                 directed=FALSE,
                                 vertices=rownames(crd))
  
  lov <- cluster_louvain(graph, weights = weights)
  lo <- norm_coords(as.matrix(crd))
  ncomms <- length(unique(lov$membership))
  
  result <- list(graph = graph,
                 layout = lo,
                 membership = lov$membership,
                 ncomms = ncomms )
  
  
  return(result)
}



# Parser ----------------------

parser <- ArgumentParser()
parser$add_argument("-c","--count_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-m","--meta_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-dm","--distance_metric",
                    type = "character",
                    default = "l2")

parser$add_argument("-nn","--n_neighbours",
                    type = "integer",
                    default = 8)

parser$add_argument("-o","--output_dir",
                    type = "character",
                    default = "/tmp/ic10")

args <- parser$parse_args()

print(args)

# Main ----------------------

# Generate specidic colormap for consistency
NCL <- 10
brwr <- brewer.pal(NCL,'Spectral')
CMAP <- list()
for (ii in 1:NCL) {
  CMAP[[as.character(ii)]] <- brwr[ii]
}

# Set variables for visualization
MNS <- 4
MXS <- 10
FIGHEIGHT <- 720
FIGWIDTH <- 720

cpth <- args$count_file
mpth <- args$meta_file
odir <- args$output_dir

# set logger details
tag <- format(Sys.time(), "%Y%m%d%H%m%s")
logfile <- paste(odir,
                 paste('STiC10-analysis',
                       tag,
                       'log',
                       sep = '.'),
                 sep='/')

flog.appender(appender.tee(logfile),
              name='ROOT')

args$n_neighbours <- ifelse(test = args$n_neighbours %in% c(4,8),
                            yes = args$n_neighbours,
                            no = 4)

flog.info(sprintf("Using %d member neighbourhood", args$n_neighbours))

# set desired distance function
if (args$distance_metric == 'cos') {
  dist_fun <- cos_dist
  flog.info("Using Cosine distance as similarity measure")
} else {
  dist_fun <- l2_dist
  flog.info("Using Euclidean distance as similarity measure")
}

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
    
    # grep for sample id
    grp <- regexpr("[0-9]{5}_[A-Z][0-9]",basename(cpth[num]))
    idx <- substr(basename(cpth[num]), grp[1], grp[1] + attr(grp,'match.length'))
  
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
    ct <- ct[inter,]
    mt <- mt[inter,]
    
    
    # Group spots by Louivain clustering
    flog.info("Perform Louvain Clustering")
    nct <- sweep(ct,1,apply(ct,1,sum),'/')
    nct[is.na(nct)] <- 0.0
    nct <- as.matrix(nct)
    xcrd <- as.numeric(unlist(mt['xcoord']))
    ycrd <- as.numeric(unlist(mt['ycoord']))
    lov_res <- louvain_clustering(xx = xcrd,
                                  yy = ycrd,
                                  datum = nct,
                                  nn = args$n_neighbours)
    
    flog.info(sprintf("Found a total of %d communities for sample %s",
                      lov_res$ncomms,idx))
    
    # Create pseudo data for joint analysis
    flog.info("Pool spots")
    # use relative frequencies in pooling
    pct <- pool_by_membership(nct,
                              lov_res$membership)
    
    # Classify spots
    flog.info('Initiating iC10 classification... ')
    res <- tryCatch(classify(t(pct)), error = function(e) F)
    
    # if classification is successfull save results    
    if (!is.logical(res)) {

      flog.info("Classification Successfull")
      
      
      posterior <- res$posterior[lov_res$membership,]
      class <- res$class[lov_res$membership]
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
          width = FIGWIDTH,
          height = FIGHEIGHT)
      
      # scale according to spot library size
      libsize <- apply(ct,1,sum)
      mxlib <- max(libsize)
      scale <- MNS + (libsize / mxlib) * (MXS - MNS)
      
      # store scaling factors
      slist[[num]] <-scale
      
      # visualize classification
      g <- ggplot(mt,
                  aes(x = xcoord,
                      y = ycoord)
                  ) +
        
        geom_point(shape = 21,
                   aes(fill = class,
                       color = tumor),
                   size = scale) +
        
        scale_colour_manual(values=c("tumor"="black",
                                     "non"="white")) +
        
        scale_fill_manual(values=CMAP) +
        
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

flog.info("Generate compressed visualization")
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

# map values to unit cube 
mx <- apply(clr, 2, max)
mn <- apply(clr,2,min)
clr <- sweep(clr,2,mn,'-')
clr <- sweep(clr,2,(mx-mn),'/')

flog.info("Save compressed visualization")
# split joint matrix and plot compressed visualization
for (pos in 1:(length(index)-1)) {
    srgb <- clr[(index[pos]+1):index[pos+1],]
    srgb = rgb(r = srgb[,1],g = srgb[,2],b = srgb[,3])
    
    filename <- paste(c(odir,gsub("\\.tsv","\\.umap\\.png",
                                  basename(mpth[pos]))),
                      collapse = "/")
    
    png(file = filename,
        width = FIGWIDTH,
        height = FIGHEIGHT)
    
    g <- ggplot(mlist[[pos]],
                aes(x = xcoord,
                    y = ycoord)) +
      
      geom_point(shape = 21, 
                 fill = srgb,
                 aes(color = tumor),
                 size = slist[[pos]]) +
      
      scale_colour_manual(values=c("tumor"="black",
                                   "non"="white")) +
      
      ggtitle(paste(c(tlist[[pos]],
                      'UMAP compression of posterior'),
                    collapse = ' '))
    
    print(g)
    dev.off()
}
flog.info("Analysis Complete")