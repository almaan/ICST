#!/usr/bin/Rscript

library(ggplot2)
library(argparse)
library(ggpubr)
library(gridExtra)


source("funcs.r")
source("utils.r")
source("spatial_funcs.r")

# Parser ---------------------------

parser <- ArgumentParser()

parser$add_argument("-c","--count_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-m","--meta_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-is","--immune_set",
                    type = "character",
                    default = "data/enr.immune.txt")

parser$add_argument("-ts","--tumor_set",
                    type = "character",
                    default = "data/enr.tumor.txt")

parser$add_argument("-o","--odir",
                    type = "character",
                    default = "/tmp/enrichment")

parser$add_argument("-pe","--plot_enrichment",
                    action = 'store_true',
                    default = F)

parser$add_argument("-pmh","--plot_mh",
                    action = 'store_true',
                    default = F)

args <- parser$parse_args()

# Variables ----------------

# get unique tag for output files
tag <- getUniqueIdentifier()

# set variables from input
cpth <- args$count_file
mpth <- args$meta_file
odir <- args$odir
tpth <- args$tumor_set
ipth <- args$immune_set

# get set of tumor associated genes
oritumst <- read.csv(tpth,
                 sep = '\n',
                 header = 1,
                 row.names = NULL)

oritumst <- as.character(unlist(oritumst))

# get set of immune associated genes
oritmest <- read.csv(ipth,
                     sep = '\n',
                     header = 1,
                     row.names = NULL)


oritmest <- as.character(unlist(oritmest))

# correct for character conversions
oritmest <- toupper(gsub("\\.","-",oritmest))
oritumst <- toupper(gsub("\\.","-",oritumst))

if (!(dir.exists(odir))) {
        dir.create(odir)
}
# set titles for plotting
titles <- c("tumor-realted",
            "immune-related")

# lists to store data
mhlist <- c()
plist <- c()
rlist <- c()

# Computations ----------------------

for (section in 1:length(cpth)) {

    # get section identifier for sample
    bname <- gsub("count_data\\-",
                  "",
                  basename(cpth[section]))
    # temporary count matrix
    ct <- read.csv(cpth[section],
                   sep = '\t',
                   row.names = 1,
                   header = 1) 

    # normalize frequencies 
    ct <- sweep(ct,1,rowSums(ct),'/')
    
    # temporary meta file
    mt <- read.csv(mpth[section],
                   sep = '\t',
                   row.names = 1,
                   header = 1,
                   stringsAsFactors = F) 

    # make covariates categorical not numerical
    mt[] <- lapply(mt, as.character)
    # get intersection of meta and count data
    inter <- intersect(rownames(mt),
                       rownames(ct))
    ct <- ct[inter,]
    mt <- mt[inter,]
    
    # remove sections with singular feature
    if (length(unique(mt$tumor)) == 1) {
        next
    }
    # store patient and replicate information 
    plist <- c(plist,mt$patient[1])
    rlist <- c(rlist,mt$replicate[1])

    # convert to uppercase to match gene sets
    colnames(ct) <- toupper(colnames(ct))
    # get intersections of sets and count data
    tumst <- intersect(colnames(ct),oritumst)
    tmest <- intersect(colnames(ct),oritmest)

    # get spotwise enrichment scores 
    tumevals <- getEnrichementScore(ct,
                                    tumst,
                                    mass =0.75)

    immevals <- getEnrichementScore(ct,
                                    tmest,
                                    mass =0.75)

    # compute and store Morisitian-Horn score
    mhlist <- c(mhlist,
                getMorisitiaHorn(tumevals,
                                 immevals)
                )
    
    # visualize enrichment
    if (args$plot_enrichment) {
        # get coordinates for plotting 
        crd <- getCoordinates(rownames(ct))
        # vector to iterate over
        feats <- list(tumevals,immevals) 
        # iterate over tumor and immune sets 
        for (num in 1:2) {
            # long-formated data frame
            plotdf <- data.frame(x = crd[,1],
                                 y = crd[,2],
                                 enr = feats[[num]],
                                 annot = mt$tumor)
            # spatial visualization 
            g1 <- ggplot(plotdf,
                         aes(x = x,
                             y = y)
                         ) + 
                
                    geom_point(size = 5,
                               stroke = 1,
                               pch = 21,
                               aes(fill = enr,
                                   color = annot)
                               ) +

                    scale_colour_manual(values=alpha(c("green",
                                                        "red"),
                                                     0.2)
                                        ) + 
            
                    scale_fill_gradient(low = "#FFFFFF",
                                        high = "#56B1F7",
                                        space = "Lab",
                                        na.value = "grey50",
                                        guide = "colourbar",
                                        aesthetics = "fill") +

                    ggtitle(titles[num])
            
            # get enrichment for different categories 
            tmrvals <- feats[[num]][mt$tumor == 'tumor']
            nonvals <- feats[[num]][mt$tumor == 'non']
            # generate long-formatted data frame 
            boxplotdf <- data.frame(annot = c(rep('tumor',length(tmrvals)),
                                              rep('non',length(nonvals))),
                                    vals = c(tmrvals,nonvals)
                                   )
            # boxplot of enrichment in tumor and non-tumor
            # spots. Including t-test result

            g2 <- ggplot(boxplotdf,aes(x = annot,
                                       y = vals)) + 

                         geom_boxplot(outlier.color = "red",
                                      aes(fill = annot)) +

                         stat_compare_means(method = "t.test")
            
            # save data to image 
            png(paste(odir,
                      paste("enrich",
                            titles[num],
                            bname,
                            tag,
                            'png',
                            sep = '.'),
                      sep = '/'),
                width = 1000,
                height = 500)

            # arrange in 1 x 2 grid
            grid.arrange(g1,
                         g2,
                         ncol = 2)

            dev.off()
         }
    }
}

# Plot Morisitian Horn -------------------

# Visualize Morisitian-Horn Index
if (args$plot_mh) {
    # set unique id for each section
    unid <- paste(plist,
                  rlist,
                  sep="_")
    # long-formatted data frame of all sections
    allplot <- data.frame(section = unid,
                          patient = plist,
                          MH = mhlist)
    # unique patient id's  
    unipatient <- unique(plist)
    # patient mean Morisitian-Horn index
    pat_mhlist <- c()
    # patient standard deviation of Morisitian-Horn index
    pat_sdlist <- c()
   
    # get statistics
    for (patient in unipatient) {
        idx <- which(plist == patient)  
        pat_mhlist <- c(pat_mhlist,
                        mean(mhlist[idx])) 

        pat_sdlist <- c(pat_sdlist,
                        sd(mhlist[idx]))
    }
    # long-foratted plot patient wise
    # including errors

    patplot <- data.frame(patient = unipatient,
                          MH = pat_mhlist,
                          err = pat_sdlist)
    
    g3 <- ggplot(allplot,
                 aes(x = section,
                     y = MH)
                 ) + 

            geom_bar(stat ="identity",
                     aes(fill = patient)
                    )
    
    g4 <- ggplot(patplot,
                 aes(x = patient,
                     y = MH)
                 ) + 

            geom_bar(stat ="identity",
                     aes(fill = patient)
                     ) +

            geom_errorbar(aes(ymin = MH - 2*err,
                              ymax = MH + 2*err))
    
    # save data to image 
    png(paste(odir,
              paste('morisitia_horn',
                    tag,
                    'png',
                    sep = '.'),
              sep = '/'),
        width = 2000,
        height = 500)
    
    # arrange in 1 x 2 grid
    grid.arrange(g3,g4,
                 ncol = 2)

    dev.off()
}
