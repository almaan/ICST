#!/usr/bin/Rscript

library(argparse)
library(iC10)
library(ggplot2)
library(futile.logger)


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

parser$add_argument("-bp","--bar_plot",
                    action = "store_true",
                    default = F)

parser$add_argument("-sp","--scatter_plot",
                    action = "store_true",
                    default = F)

args <- parser$parse_args()

iC10GENES <- read.csv("iC10.gene_symbols.tsv",
                        sep = '\t',
                        header = 1,
                        stringsAsFactors = F)

iC10GENES <- as.vector(unlist(iC10GENES))

cpth <- sort(args$count_file)
mpth <- sort(args$meta_file)

clist <- list()
mlist <- list()
plist <- c()

nsamples <- length(cpth)

for (sample in 1:nsamples) {
    
    flog.info(sprintf("Analyzing sample %s | %d/%d",cpth[sample],sample,length(cpth)))

    ct <- read.table(cpth[sample],
                     sep = '\t',
                     header = 1,
                     row.names = 1)
    #load meta data
    mt <- read.table(mpth[sample],
                     sep = '\t',
                     header = 1,
                     row.names = 1)
    
    
    # select only spots found in meta and count data
    rinter <- intersect(rownames(ct),rownames(mt))
    ct <- ct[rinter,]
    mt <- mt[rinter,]

    ct <- sweep(ct,1,rowSums(ct),'/')
    ct[is.na(ct)] <- 0.0
    cinter <- intersect(colnames(ct), iC10GENES)
    nct <- as.data.frame(matrix(0, nrow(ct), length(iC10GENES)))
    rownames(nct) <- rownames(ct)
    colnames(nct) <- as.character(iC10GENES)
    nct[,cinter] <- ct[,cinter]
    clist[[sample]] <- nct
    mlist[[sample]] <- mt
    plist <- c(plist,mt$patient[1])
}

    
if (args$bar_plot) {
    for (sample in 1:length(clist)) {
        flog.info("")
        for (class in c("tumor","non")) {
            ct <- clist[[sample]]
            mt <- mlist[[sample]]

            exprmean <- colMeans(ct[mt$tumor == class,])
            q90 <- quantile(exprmean,0.9)
            top <- rep(0,length(exprmean))
            top[exprmean > q90] <- 1
            
            spot <- rownames(ct)
            genename <- colnames(ct)
            
            stdep <-  exprmean - 2.0 * apply(ct[mt$tumor == class,],2,sd)
            stden <-  exprmean + 2.0 * apply(ct[mt$tumor == class,],2,sd)
            
            df <- data.frame( expr = exprmean,
                              top = top,
                              stdep =  stdep,
                              stden = stden,
                              genename = genename)
            
            g <- ggplot(data = df, aes(x = genename, y = expr)) +
                    geom_bar(stat = "identity", aes(fill = top)) + 
                    geom_errorbar(aes(x = genename, ymin = stden,ymax = stdep)) + 
                    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))
            
            oname <- paste("/tmp/ic10.genes",class,"png", sep = '.')
            png(oname, width = 8 * length(iC10GENES), height = 500)
            print(g)
            dev.off()
        }
    }
}


if (args$scatter_plot) {
    n_patients <- length(unique(plist))
    n_genes <- length(iC10GENES)
    exprvec <- c()
    varvec <- c()
    pvec <- c()
    cvec <- c()
    gvec <- c()
    
    for (patient in plist) {
        for (class in c("tumor","non")) {
            pexpr <- c()
            for (sample in which(plist == patient)) {
                pexpr <- rbind(pexpr,clist[[sample]][mlist[[sample]]$tumor == class,]) 
            }

            exprvec <- c(exprvec,colMeans(pexpr))
            varvec <- c(varvec,apply(pexpr,2,var))
            pvec <- c(pvec,rep(patient,dim(pexpr)[2])) 
            cvec <- c(cvec,rep(class, dim(pexpr)[2]))
            gvec <- c(gvec,colnames(clist[[1]]))
    
        }
    }

    df <- data.frame(expr = exprvec,
                     var = varvec,
                     class = as.factor(cvec),
                     patient = as.factor(pvec),
                     genes  = gvec
                     )
    
    
    goi <- df$expr >= sort(df$expr,decreasing = T)[min(length(df$expr),100)]
    df$goi <- goi
    
    df <- df[sample(nrow(df)),]    

    g <- ggplot(data = df, aes(x = var, y = expr)) + 
            geom_point(aes(fill = patient,
                                  shape = class
                                  ),
                       alpha = 0.5) + 

            scale_shape_manual(values=c(tumor = 23,
                                        non = 21)) +

            annotate("text", x = df$var[df$goi],
                             y = df$expr[df$goi],
                             label = df$genes[df$goi],
                             hjust = 0,
                             vjust = 0)
    
    imgname <- paste(args$output_dir,paste("expr_v_var","png", sep = '.'),sep = '/')

    png(imgname, width = 2000, height = 2000)
    print(g)
    dev.off()

    fname <- paste(args$output_dir,"expr_v_var.tsv", sep = '/')
    write.table(unique(df$gene[df$goi]), fname,  sep = '\t', row.names = F, col.names = F)
}



 
