#!/usr/bin/Rscript

# Computes average correlation between spots with a 
# paiwise distance d < r, where r is a user defined radius
# The correlation is compared between tumor and non tumor
# spots 

library(RANN)
library(dabestr)
library(argparse)
library(futile.logger)

source("utils.r")

# Parser -----------------

parser <- ArgumentParser()

parser$add_argument("-c","--count_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-m","--meta_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-o","--output_dir",
                    type = "character",
                    default = "/tmp/hetero")

parser$add_argument("-r","--radius",
                    default = '1.1',
                    type = "double"
                    )


args <- parser$parse_args()

# Prepare Data -------------------

cpths <- sort(args$count_file)
mpths <- sort(args$meta_file)

radius = args$radius
feature <- 'tumor'
fvals <- c("tumor","non")
plist <- c()
sublist <- c()
res <- list()

# Main -------------------------

for (section in 1:length(cpths) ) {
    dta <- readSTsection(cpths[section],
                         mpths[section]
                        )

    patient_id <- dta$meta_data$patient[1]
    
    dta$count_data <- sweep(dta$count_data,
                            1,
                            rowSums(dta$count_data),
                            '/')
    
    dta$count_data[is.na(dta$count_data)] <- 0.0

    uni <- table(dta$meta_data$tumor)

    if (any( as.numeric(uni) < 10 ) | length(uni) < 2) { 
        flog.info(sprintf("Sample %s will be discarded. nTumorspots : %d  | nNonTumorspots : %d",
                  cpths[section],
                  uni['tumor'],
                  uni['non']
                  )
                )
        next
    }

    plist <- c(plist,patient_id)
    sublist <- c(sublist,dta$meta_data$subtype[1])

    crd <- getCoordinates(rownames(dta$count_data))
    
    corref <- as.data.frame(matrix(0,
                                   nrow = nrow(dta$count_data),
                                   ncol = 2
                                   )
                            )
    
    colnames(corref) <- c("cdiff","type")
    rownames(corref) <- rownames(dta$count_data)
    corrmat <- cor(t(dta$count_data),
                   method = 'spearman')
    flog.info(sprintf("Correlation matrix for %s computed",
                      cpths[section])) 

    for (val in fvals ) {
        flog.info(sprintf("Computing neighbours within %s spots",
                          val))

        idx <- which(dta$meta_data[feature] == val)

        nbr <- nn2(crd[idx,],
                  searchtype = "radius", 
                  radius = radius,
                  k = 9,
                  )$nn.idx

        nbr <- nbr[,-1]
    
        wnbr <- which(!(rowSums(nbr) == 0))
        wonbr <- which(rowSums(nbr) == 0)
        corref <- corref[-wonbr,]
        nwnbr <- length(wnbr)
    
        for (spot in wnbr) {
            nz <- idx[nbr[spot,!(nbr[spot,] == 0)]]
            nnz <- length(nz)
            mureal <- mean(corrmat[idx[spot],nz])
            elegible <- idx[!(idx == idx[c(spot,nbr[spot,])])]

            if (length(elegible) < 1 ) {
                next
            }
    
            sname <- rownames(dta$count_data)[idx[spot]]
            corref[sname,"cdiff"] <- mureal 
            corref[sname,'type'] <- val
        }

        rownames(corref) <- paste(section,rownames(corref), sep = '_')
        res[[section]] <- corref
    }
}

plist <- as.character(plist)
sublist <- sublist[!(duplicated(plist))]
unipat <- unique(plist) 
names(sublist) <- unipat
txtres <- as.data.frame(matrix(0,
                               nrow = length(unipat),
                               ncol = 5
                               )
                        )

colnames(txtres) <- c("mean",
                      "ci_low",
                      "ci_high",
                      "has_zero",
                      "subtype"
                      )

rownames(txtres) <- unipat

# Visualize -----------------------------

for (pat in unipat) {
    
    tdf <- as.data.frame(dplyr::bind_rows(res[which(plist == pat)]))

    if (!(all(c('tumor','non') %in% tdf$type ))) {
            next
    }

    opth <- paste(args$output_dir,
         '/',
         "test.",
         as.character(pat),
         '.png',
         sep = '')

    
    
    pltme <- dabestr::dabest(tdf, type, cdiff,
                             idx = c("non", "tumor"),
                             paired = FALSE)

 
    txtres[pat,c('mean','ci_low','ci_high')] <- as.numeric(c(pltme$result$difference,
                                                             pltme$result$bca_ci_low,
                                                             pltme$result$bca_ci_high
                                                             ))

    haszero <- ( (txtres[pat,'ci_low'] <= 0) &
                (txtres[pat,'ci_high'] > 0)) 

    txtres[pat,'has_zero'] <- ifelse(haszero,
                                     yes =  1,
                                     no = 0) 


    txtres[pat,'subtype'] <- sublist[pat]

    p <- plot(pltme)
    png(opth)
    print(p,
          palette = 'Pastel1',
          theme = ggplot2::theme_minimal()
          )

    dev.off()

}

# Save Tabulated Result ---------
topth <- paste(args$output_dir,
               '/',
               'summary.tsv',
               sep = ''
               )

write.table(txtres,topth, sep = '\t', quote = F)

