#!/usr/bin/Rscript
library(argparse)
library(futile.logger)
library(future.apply)


source("spatial_funcs.r")
source("funcs.r")
source("utils.r")

# Parser ---------------------------

parser <- ArgumentParser()

parser$add_argument("-c","--count_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-m","--meta_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-o","--odir",
                    type = "character",
                    default = "/tmp/autocorr")

parser$add_argument("-ng","--n_genes",
                    type = "integer",
                    default = NULL)

parser$add_argument("-ls","--length_scale",
                    default = 0.5)

parser$add_argument("-mc","--max_clusters",
                    type = "integer",
                    default = NULL)

parser$add_argument("-nw","--n_workers",
                    type = "integer",
                    default = 1)

parser$add_argument("-pm","--pre_morans",
                    type = "character",
                    default = NULL
                    )

parser$add_argument("-tq","--top_quantile",
                    type = "double",
                    default = 0.95
                    )



args <- parser$parse_args()

# Session Information ------------

# generate unique identifier tag for analysis
tag <- getUniqueIdentifier()
# check if specified output directory
# exists. Create it if not

checkMakeDir(args$odir)
# Set number of workers
options(future.globals.maxSize= 10000*1024^2)
plan(multiprocess,
     workers = args$n_workers) 

flog.info(sprintf("Using %d workers",args$n_workers))


# Load Data ----------------

# load specified data
dta <- load_multiple(cpth = args$count_file,
                     mpth = args$meta_file)

flog.info("Loaded data")

# select top "ngenes" expressed genes if
# specified by user. Mainly for evaluation.
if (!(is.null(args$n_genes))) {
    ngenes <- min(c(ncol(dta$count_data),args$n_genes))
    selidx <- order(colSums(dta$count_data),
                    decreasing = T)[1:ngenes]
# else use all genes
} else {
    selidx <- c(1:ncol(dta$count_data))
}

dta$count_data <- dta$count_data[,selidx]

# Main ---------------------------

# transform observed data to relative frequencies
dta$count_data <- sweep(dta$count_data,
                        1,
                        rowSums(dta$count_data),'/')

# adjust for eventual zero expression
dta$count_data[is.na(dta$count_data)] <- 0.0

# get coordinates for all spots
crd <- getCoordinates(gsub(pattern = '[0-9]*_',
                           replace = '',
                           x = rownames(dta$count_data)))

flog.info("Computed Weights Matrix")
# generate weights matrix
wmat <- getSpatialWeights(crd = crd,
                          startpos = dta$startpos,
                          sigma = args$length_scale)

diag(wmat) <- 0.0

# compute Moran's I for all genes
dta$count_data <- sweep(dta$count_data,
                        2,
                        colMeans(dta$count_data),
                        '-'
                        )

wmat <- sweep(wmat,
              1,
              rowSums(wmat),
              '/'
              )

if (is.null(args$pre_morans)) {
    flog.info("Computing Moran's I")
    res <- getMoransIlarge(nfmat = as.matrix(dta$count_data),
                           wmat = as.matrix(wmat))
} else {
    res <- read.table(args$pre_morans,
                      sep = '\t',
                      row.names = 1,
                      header = 1,
                      stringsAsFactors = F
                      )
}

#res <- getMoransIv2(nfmat = as.matrix(dta$count_data),
#                    wmat = as.matrix(wmat))

# set rownames to gene names
rownames(res) <- colnames(dta$count_data,)
# remove eventual nans
res <- res[!(is.na(res$I)),]

# Path for Moran's I values
morans_opth <- paste(args$odir,
                     paste(tag,
                           "moransi",
                           'tsv',
                           sep = '.'),
                     sep = '/')

# save Moran's I values 
write.table(res,
            morans_opth,
            sep = '\t',
            quote = F,
            col.names = T)

# get values of top quantile
args$top_quantile <- ifelse(args$top_quantile < 1,
                            yes = args$top_quantile,
                            no = 0.95)

qv <- as.numeric(quantile(res$I,args$top_quantile))

# select only top genes with positive autocorrelation
ac_genes <- res[res$I > 0.0 & res$I > qv,]
ac_genes <- as.character(rownames(ac_genes))

# generate correlation based dissimilarity  matrix
flog.info("Generate dissimilarity matrix")
dmat <- getCorr(dta$count_data[,ac_genes])

# perform hierchical clustering and use
# Dunn's index for finding optimal cluster number
flog.info("Clustering data and evaluating optimal cluster number")
dres <- clusterWithDunn(dmat,
                        maxClusters = args$max_clusters)

cluster_labels <- dres$labels

# Save Results -------------------

# set base of filename for results
base_opth <- paste(args$odir,
                    paste(tag,
                          "autocorr.clustering",
                          sep = '.'),
                    sep = '/')
# image filename
image_opth <- paste(base_opth,
                    'png',
                    sep = '.')

# table filename
table_opth <- paste(base_opth,
                    'tsv',
                    sep = '.') 

flog.info(sprintf("Saving text based results to >> %s",
                  table_opth))

# save text format of results
write.table(data.frame(cluster = cluster_labels),
            table_opth,
             sep = '\t',
             quote = F,
             row.names = rownames(dmat),
             col.names = T)

flog.info(sprintf("Saving visualization of result to >> %s",
                  image_opth))

# assemble heatmap object
phm <- makeDistanceBasedHeatMap(dmat,
                                cluster_labels)


# save image format of results
png(image_opth,
    width = length(ac_genes) * 2 + 5,
    height = length(ac_genes) * 2 + 10)

print(phm)

dev.off()

flog.info(sprintf("Saved heatmap of result to >> %s",
                  image_opth))

