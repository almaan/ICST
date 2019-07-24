#!/usr/bin/Rscript
library(argparse)
library(futile.logger)
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

args <- parser$parse_args()

# Session Information ------------

# generate unique identifier tag for analysis
tag <- getUniqueIdentifier()
# check if specified output directory
# exists. Create it if not

checkMakeDir(args$odir)

# Load Data ----------------

# load specified data
dta <- load_multiple(cpth = args$count_file,
                     mpth = args$meta_file)

# select top "ngenes" expressed genes if
# specified by user. Mainly for evaluation.
if (!(is.null(args$n_genes))) {
    ngenes <- min(c(ncol(dta$count_data),args$n_genes))
    selidx <- order(colSums(dta$count_data),
                    decreasing = T)[1:ngenes]
    print(selidx)
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

# generate weights matrix
wmat <- getSpatialWeights(crd = crd,
                          startpos = dta$startpos,
                          sigma = args$length_scale)

# compute Moran's I for all genes
flog.info("Computing Moran's I")
res <- getMoransIlarge(fmat = dta$count_data,
                       wmat = wmat)

# set rownames to gene names
rownames(res) <- colnames(dta$count_data,)
# remove eventual nans
res <- res[!(is.na(res$I)),]

# get values of top quantile
qv <- as.numeric(quantile(res$I,0.95))

# select only top genes with positive autocorrelation
ac_genes <- res[res$I > 0.0 & res$I > qv,]
ac_genes <- as.character(rownames(ac_genes))

# generate correlation based distance matrix
flog.info("Computing correlations")
dmat <- getCorr(dta$count_data[,ac_genes])

# perform hierchical clustering and use
# Dunn's index for finding optimal cluster number
dres <- clusterWithDunn(dmat,
                        maxClusters = args$max_clusters)
cluster_labels <- dres$labels

# Save Results -------------------

# assemble heatmap object
phm <- makeDistanceBasedHeatMap(dmat,
                                cluster_labels)

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

# save text format of results
write.table(data.frame(cluster = cluster_labels),
            table_opth,
             sep = '\t',
             quote = F,
             row.names = rownames(dmat),
             col.names = T)

# save image format of results
png(image_opth,
    width = length(ac_genes) * 50 + 300,
    height = length(ac_genes) * 50 + 100)
print(phm)
dev.off()
