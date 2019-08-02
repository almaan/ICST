#!/usr/bin/Rscript

library(glmnet)
library(caret)
library(corrplot)

library(argparse)

source("utils.r")


source("utils.r")
source("funcs.r")

makeCorrDistanceMatrix <- function(res,
                               ascorr = F) {

    joint <- as.data.frame(matrix(0,
                                  nrow = length(res$genesets),
                                  ncol = length(res$tumorgenes)
                                  )
                            )
    rownames(joint) <- names(res$genesets)
    colnames(joint) <- as.character(res$tumorgenes)

    for (pat in names(res$genesets)) {
        joint[pat,as.character(res$genesets[[pat]]$genes)] <- res$genesets[[pat]]$val
    }

    dmat <- cor(t(joint),
            method = 'spearman')

    if (!(ascorr)){
        dmat <- 1.0 - dmat
    }

    return(dmat)
}


clusterPatients <- function(dmat) {
    cl <- cluster::agnes(dmat,
                         diss = T,
                         method = 'ward')

    return(cl)
}


makeFisherDistanceMatrix <- function(res) {
    npat <- length(res$genesets)
    pnames <- names(res$genesets)

    dmat <- as.data.frame(matrix(0,
                                 nrow = npat,
                                 ncol = npat
                                 )
                         )

    rownames(dmat) <- colnames(dmat) <- as.character(pnames)
    
    for ( ii in 1:(npat-1)) {
        for (jj in (ii+1):npat) {
            dmat[ii,jj] <- dmat[jj,ii] <- doHyperGeometric(res$genesets[[pnames[ii]]]$genes,
                                                           res$genesets[[pnames[jj]]]$genes,
                                                           res$allgenes)
        }
    }
    
    diag(dmat) <- NA

    dmat <- -1.0 / log(dmat)
    diag(dmat) <- 0.0

    
    return(dmat)
}

makeFisherDistanceMatrix <- function(res) {
    npat <- length(res$genesets)
    pnames <- names(res$genesets)

    dmat <- as.data.frame(matrix(0,
                                 nrow = npat,
                                 ncol = npat
                                 )
                         )

    rownames(dmat) <- colnames(dmat) <- as.character(pnames)
    
    for ( ii in 1:(npat-1)) {
        for (jj in (ii+1):npat) {
            dmat[ii,jj] <- dmat[jj,ii] <- doHyperGeometric(res$genesets[[pnames[ii]]]$genes,
                                                           res$genesets[[pnames[jj]]]$genes,
                                                           res$allgenes)
        }
    }
    
    diag(dmat) <- 0.0
    
    return(dmat)
}


# Parser -----------------

parser <- ArgumentParser()

parser$add_argument("-c","--count_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-m","--meta_file",
                    type = "character",
                    nargs = "+")

args <- parser$parse_args()

cpths  <- sort(args$count_file)
mpths <- sort(args$meta_file)

dta <- load_multiple(cpths,
                     mpths)

dta$count_data <- as.matrix(dta$count_data)

dta$count_data[is.na(dta$count_data)] <- 0.0

res <- doPatientWiseLogisticRegression(dta$count_data,
                                       dta$meta_data)

save(res, file = '/tmp/patLogReg.R')


unipat <- unique(dta$meta_data$patient)
dta$bulk_data <- sapply(unipat,
                          function(x) { colMeans(dta$count_data[(dta$meta_data$patient == x) &
                                                                (dta$meta_data$tumor == 'tumor'),
                                                 res$tumorgenes])
                                      }
                       )

std <- apply(dta$bulk_data, 2, sd)
mu <- apply(dta$bulk_data,2, mean)


#dta$bulk_data <- sweep(dta$bulk_data,
#                        2,
#                        mu,
#                        '-'
#                        )
#
#dta$bulk_data <- sweep(dta$bulk_data,
#                        2,
#                        std,
#                        '/'
#                        )
#

dmat <- cor(dta$bulk_data,
            method = "spearman"
            )

dmat <- 1 - dmat
clustering <- clusterPatients(dmat)
dmat <- 1 - dmat

save(dmat,file = '/tmp/dmat.R')

pdf('/tmp/corr.plot.pdf')
par(mfrow = c(1,2),
    mar = c(5,6,4,6)
    )

#corrplot(dmat[rev(clustering$order),clustering$order],
#         method = 'circle',
#         col = rev(colorspace::diverge_hcl(n = 100))
#         )

heatmap(dmat[rev(clustering$order),clustering$order],
        Colv = NA,
        Rowv = NA,
        scale = 'none',
        col = rev(colorspace::diverge_hcl(n = 100))
        )

cluster::pltree(clustering)
dev.off()
