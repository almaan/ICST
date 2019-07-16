#!/usr/bin/Rscript
sh <- suppressPackageStartupMessages
sh(library(iC10))
sh(library(argparse))
sh(library(ggplot2))
sh(library(futile.logger))

# Functions ------------------

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
  
  features <- matchFeatures(Exp=x,
                            Exp.by.feat=feat)

  features <- normalizeFeatures(features,
                                method)
  res <- iC10(features)
  
  return(res)
} 


# Parser -----------------------

parser <- ArgumentParser()

parser$add_argument("-c","--count_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-m","--meta_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-s","--select_for",
                    type = "character",
                    default = NULL)

parser$add_argument("-o","--output_dir",
                    type = "character",
                    default = "/tmp")

parser$add_argument("-gs","--gene_set",
                    type = "character",
                    default = NULL)


args <- parser$parse_args()

# Main -----------------------------

if (!(is.null(args$gene_set))) {
    gene_set <- read.csv(args$gene_set,
                         sep = '\t')
} else {
    gene_set <- NULL
}


# set variables
select_for <- args$select_for
cpth <- args$count_file
mpth <- args$meta_file
# set output name
ofilename <- paste(gsub(":|-| |","",
                        as.character(Sys.time())),
                        "bulk_analysis",
                        sep = '.')

opth <- paste(args$output_dir,ofilename, sep = '/')
print(opth)

# make sure count and meta data are matched
cpth <- sort(cpth)
mpth <- sort(mpth)

# to store information
clist <- list()
mlist <- list()
plist <- c()
unigenes <- c()

tstart <- as.numeric(Sys.time())
tend <- tstart

for (section in 1:length(cpth)) {
    
    # temporary count matrix
    ct <- read.csv(cpth[section],
	               sep = '\t',
	               row.names = 1,
	               header = 1) 

    # temporary meta file
    mt <- read.csv(mpth[section],
	               sep = '\t',
	               row.names = 1,
	               header = 1) 

    # only use spots present in both meta and count data
    inter <- intersect(rownames(ct), rownames(mt))
    ct <- ct[inter,]
    mt <- mt[inter,]
    
    # select specified spots
    if (!(is.null(select_for))) {
            flog.info(sprintf("Select only spots with annotation : %s",select_for))
            idx <- which(mt$tumor == select_for)
            ct <- ct[idx,]
            mt <- mt[idx,]
    }
    
    # get relative frequencies of counts
    ct <- sweep(ct,1,rowSums(ct),'/')

    if (!(is.null(args$gene_set))) {
            gene_set <- read.csv(args$gene_set,
                                 sep = '\t')
    
            cinter <- intersect(gene_set,colnames(ct)) 
            ct <- ct[,cinter]
    }

    # store count and meta data
    clist[[section]] <- ct
    mlist[[section]] <- mt
    
    # store patient information
    plist <- c(plist,as.character(mt$patient[1]))
    # update list of observed genes
    unigenes <- union(unigenes,colnames(ct))
    
    # print progress
    tend <- as.numeric(Sys.time())
    eta <- (tend - tstart) / section * (length(cpth) - section ) / 60.0 
    flog.info(sprintf(" Loaded section %s | Progress :  %d / %d |  ETA %s min",
                      cpth[[section]],
                      section,
                      length(cpth),
                      eta))

}


# number of genes present in data 
ngenes <- length(unigenes)
# patient names 
unipatient <- unique(plist)
# number of patients 
npatients <- length(unipatient)

# empty dataframe to fill for patient pooled data
pcnt <- as.data.frame(matrix(0,
                             nrow = npatients,
                             ncol = ngenes)
                     ) 

colnames(pcnt) <- unigenes

flog.info("Pooling sections from the same patient")
for (patient in 1:npatients) {
    # get which sections that belong to patient
    idx <- which(unipatient[patient] == plist)
    # number of sections for patient
    nsections <- length(idx)
    # pool sections
    for (section in idx ) {
       pcnt[patient,colnames(clist[[section]])] <- pcnt[patient,colnames(clist[[section]])] + 
                                                   (colMeans(clist[[section]]) / (nsections))
    }
}

# classify patient pooled data
ic10_res <- classify(t(pcnt))

# Visualize -----------------------------------

# get number of classes
nclasses <- ncol(ic10_res$posterior)

# prepare columns for long-format used in ggplot

# patient id
patient_id <-  (unlist(lapply(unipatient,
                             function(x) {rep(x,nclasses)})
                         ))
# posterior probabilities
posteriors <- sapply(t(ic10_res$posterior),
                     function(x) { c(x) })
# class names
classname <- rep(paste("IC", seq(1,nclasses), sep = ''),
             npatients)

# create long-formatted dataframe
plotdf <- data.frame(patient = patient_id,
                     IC = classname,
                     y = posteriors)

write.table(plotdf,
          paste(opth,"tsv",sep = '.'),
          sep = '\t',
          col.names = T,
          row.names = T,
          quote = F)

flog.info(sprintf("Saving plot to : %s", opth))

# plot data
g <- ggplot(plotdf,
            aes(fill=IC,
                y=y,
                x=patient)) + 

        geom_bar(stat = "identity",
                 position = "fill")

png(paste(opth,"png", sep = "."),
    width = 2100,
    height = 500)

print(g)

flog.info("Analysis successfully complete")

dev.off()


    
 

