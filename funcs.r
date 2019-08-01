#!/usr/bin/Rscript

chr_dist <- function(cnt,
                     gene_mapping,
                     nbins = 100,
                     edges = NULL,
                     binsize = NULL) {

    # Map observered transcript to binned chromosome
    # position
    #
    # Args:
    #   cnt - count matrix (n_spots x n_genes)
    #   gene_mapping - dataframe with genes as rows,
    #                   chromosomal start/end position
    #                   and chromosome number as columns
    #   nbins - number of bins to use
    #   edges - n_bins x 2 data frame with lower edges
    #               stored in column "ledge" and upper
    #               edges stored in column "uedge".
    #   binsize - desired binsize. Use as alternative to
    #               nbins
    #
    # Return:
    #   list with bin information (upper and lower edge)
    #   and a (n_spots x n_bins) matrix with the
    #   number of observed transcripts in each bin
    #   for each spot

    
    # set gene to chromosome position mapping object
    G2CHR <- gene_mapping

    inter <- intersect(colnames(cnt),rownames(G2CHR))
    success_map <- length(inter) / ncol(cnt)

    if (success_map >= 0.5) {
        print(sprintf("%f %% of all genes were successfully mapped to chromosomal position",
                      success_map * 100))
    } else {
        print("WARNING : less than 50% of all genes were successfully mapped")
    }
    
    cnt <- cnt[,inter]
    G2CHR <- G2CHR[inter,]
 
    suppressPackageStartupMessages(require(GenomicRanges))

    # load chromosomal information
    data(hg19IdeogramCyto,
         package = "biovizBase")
    data(hg19Ideogram,
         package = "biovizBase")

    chrs <- as.character(levels(seqnames(hg19IdeogramCyto)))
    seqlths <- seqlengths(hg19Ideogram)[chrs]    

    # order information
    chrs <- gsub('X','23',chrs)
    chrs <- gsub('Y','24',chrs)
    chrs <- chrs[order(as.numeric(gsub('chr','',chrs)))]
    chrs <- gsub('23','X',chrs)
    chrs <- gsub('24','Y',chrs)
    seqlths <- as.numeric(seqlths[chrs])
    names(seqlths) <- chrs
   
    cumseqlths <- cumsum(c(0,as.numeric(seqlths[-length(seqlths)])))
    names(cumseqlths) <- chrs

    # get chromosomal position for genes  
    startpos <- G2CHR[colnames(cnt),"start_position"] 
    endpos <- G2CHR[colnames(cnt),"end_position"]
    chrpos <- paste("chr",G2CHR[colnames(cnt),"chromosome_name"],sep = "")

    # transform to absolute position 
    startpos <- startpos + as.numeric(cumseqlths[chrpos])
    endpos <- endpos + as.numeric(cumseqlths[chrpos])
   
    # setup edges for bins
    xmin <- 0
    xmax <- as.numeric(tail(cumseqlths,1))
  
    if (is.null(edges)) {
        if (is.null(binsize)) {
        # generate bins
            ledge <- seq(xmin,xmax, length.out = (nbins + 1))
        } else {
            ledge <- seq(xmin,xmax, by = binsize)
            nbins <- length(ledge) - 1
        }
        uedge <- ledge + diff(ledge)[1]
        ledge <- ledge[-length(ledge)]
        uedge <- uedge[-length(uedge)]
     } else {
        # use user-defined edges
        nbins = nrow(edges)
        ledge <- as.numeric(unlist(edges$ledge))
        uedge <- as.numeric(unlist(edges$uedge))
     }
    # get number of transcripts falling into each bin
    bincount <- matrix(0,
                       nrow = nrow(cnt),
                       ncol = nbins)
    
    for (gene in 1:ncol(cnt)) {
        # conditions for bin membership
        c1 <- (startpos[gene] >= ledge) & (startpos[gene] < uedge)
        c2 <- (endpos[gene] >= ledge) & (endpos[gene] < uedge)
        c3 <- (startpos[gene] <= ledge) & (endpos[gene] >= uedge)
        # get which bins the gene cover
        inbin <- which( c1 | c2 | c3)
        bincount[,inbin] <- bincount[,inbin] + cnt[,gene]
    }
   
    # format bincount object
    bincount <- as.data.frame(bincount)
    colnames(bincount) <- paste(ledge,
                                uedge,
                                sep = '_')

    rownames(bincount) <- rownames(cnt)
    
    bins <- data.frame(start = ledge, end = uedge)
    
    result <- list(bins = bins,
                   bincount = bincount)
    
    return(result)
}


make_g2chr <- function(genes,
                       genes_sep = '\n',
                       dataset="hsapiens_gene_ensembl",
                       host="grch37.ensembl.org",
                       path="/biomart/martservice"){
    
    # Make reference data frame which allows for
    # mapping of a gene (by hgnc symbol) to relative
    # chromosome position
    # Args:
    #   genes : character vector containing genes 
    #               to be included in reference
    #   dataset : biomart dataset
    #   host    : biomart host  
    #   path    : biomart path
    # Return:
    #   dataframe with columns "hgnc_symbol","start_position", 
    #       "end_position" and "chromosome_name". To be used as
    #       a reference when mapping genes to chromosomes

    library(biomaRt)

    genes <- unlist(genes)
    
    # create mart object
    mart  <- useMart("ensembl",
                     dataset=dataset,
                     host=host,
                     path=path)
    
    # specify attributes and flters to use
    attribs <- c("hgnc_symbol",
                 "start_position",
                 "end_position",
                 "chromosome_name")

    filters <- c("hgnc_symbol")
    
    # retrieve attributes for genes
    mapped <- getBM(attributes = attribs,
                    filters = filters,
                    values = genes,
                    mart = mart)

    # remove duplicate hits. Keeps first
    mapped <- mapped[!(duplicated(mapped$hgnc_symbol)),]
    
    # use intersecting set of genes
    inter <- intersect(genes,
                       mapped$hgnc_symbol)

    mapped <- mapped[mapped$hgnc_symbol %in% inter,]
    genes <- genes[inter]
    
    # remove genes with unvalid chromosome name
    valnames <- c(as.character(c(1:22)),"X","Y")
    validchr <- mapped$chromosome_name %in% valnames
    mapped <- mapped[validchr,]
    
    # set rownames to symbols
    rownames(mapped) <- mapped$hgnc_symbol
    
    return(mapped)
}
 

getEnrichementScore <- function(cnt,
                                gset,
                                mass = 0.8,
                                tset = NULL) {
    # Computes enrichement score for a given
    #   set of spots with known expression.
    #
    # Args: 
    #   cnt - (n_spots x n_genes) coint matrix
    #   gset - vector with genes to assess
    #            enrichment of 
    #   mass - what fraction of the transcripts
    #           the top set should make up
    # 
    #   tset - total set of genes. If not provided
    #            use all genes present in count matrix
    # Returns:
    #   Enrichment score for each spot
    
    getTopNames <- function(sct,
                            names,
                            mass = 0.8) {
        # Get names of the genes whose counts
        #   contsitute a given (mass) proportion
        #   of the total number of transcripts
        #   in a given spot
        #
        # Args: 
        #   sct - (n_genes,) vector containing the
        #       expression vector of a certain spot
        #   names - (n_genes,) name of genes, matched
        #       with sct
        #   mass - what fraction of the transcripts
        #           the top set should make up
        # Returns:
        #   The set of genes containing a given
        #    proportion of the total number of
        #    observed transcripts
 

        # get ordered indices
        idx <- order(sct,decreasing = T) 
        # compute cummulative sum
        sel <- cumsum(sct[idx])
        # get value corresponding to desired
        # proportion of transcript mass
        q <- mass * sel[length(sel)]

        return(names[idx][sel <= q])
    
    }
    
    
    
    doHyperGeometric <- function(sptset,
                                 gset,
                                 tset) {

        # Conduct a hyper-geometric test
        # to obtain pvalue of intersection size
        # between two sets
        #
        # Args:
        #   sptset - vector with gene names
        #       associated to a spot. Top genes,
        #   gset - vector with set of names
        #           for genes of interest
        #   tset - vector with set of names
        #           for all observed genes
        # Returns:
        #   Pvalue of having an intersection 
        #       equal or larger than the
        #       observed value
        
        # get intersection of set of intereset
        # and spot associated set
        inter <- intersect(sptset,gset)
        # get cardinality of full set 
        N <- length(tset)
        # cardinality of spotset 
        m <- length(sptset)
        # cardinality of complement to spotset
        n <- N - m 
        # cardinality of set of intereset
        k <- length(gset)
        # cardinality of intersection
        x <- length(inter) 
        # do hypergeomatric test
        p <- phyper(x-1, m,n,k, lower.tail = F, log.p = F) 
        
        return (p)
    }
    
    if (is.null(tset)) {
        tset <- colnames(cnt)
    }
    
    # get spot associated sets
    spotsets <- apply(cnt,
                      1,
                      getTopNames,
                      tset,
                      mass)
    # compute pvalues
    pvals <- unlist(lapply(spotsets,
                           doHyperGeometric,
                           gset,
                           tset))
    # compute enrichement score
    escore <- -log(pvals)
    # correct for eventual zeros
    escore[is.na(escore)] <- max(escore,na.rm = T)

    return(escore)
    
}

build_bins <- function(bin.size=1e+6, chrs=NULL, outfile=NULL){


    # build bins for genome hg19
    # this function create bins by moving binsize window across the genome until meeting
    # chromosome centromere and telomere. When reach the midpoints and endpoints,
    # a mid/endpoints-bin.size bin will be add.
    #
    # Args:
    #   bin.size - size of the bins, default: 1e+6
    #   chrs - chromosome names that need to be binned
    #          if chrs=NULL, then the entire genome will be bined
    #   outfile - if NULL, the output will be returned as a dataframe
    #             otherwise, the output will be written into outfile
    # Return:
    #   dataframe with column names (chr, arm, start, end, start_gps, end_gps)
    #   start_gps and end_gps are the absoulte genomic positions

    library(GenomicRanges)

    # load chromosomal information
    data(hg19IdeogramCyto,
         package = "biovizBase")
    data(hg19Ideogram,
         package = "biovizBase")

    ordered.chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5",
                   "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11","chr12","chr13","chr14","chr15",
                   "chr16","chr17","chr18","chr19","chr20",
                   "chr21","chr22","chrX", "chrY")

    if (is.null(chrs)){
        chrs <- ordered.chrs
    }

    seqlths <- seqlengths(hg19Ideogram)[chrs]
    midpoints <- sapply(chrs, function(x) max(end(hg19IdeogramCyto[seqnames(hg19IdeogramCyto)==x & grepl("p",hg19IdeogramCyto$name),])))

    all.bins <- c()

    for (c in chrs){
        p <- seq(0, midpoints[c], by=bin.size)
        q <- seq(midpoints[c], seqlths[c], by=bin.size)

        bins <- data.frame(chr = rep(c, length(c(p,q))),
                            arm = c(rep("p",length(p)), rep("q",length(q))),
                            start = c(p[-length(p)], max(as.numeric(midpoints[c]-bin.size),p[1]), q[-length(q)], max(as.numeric(seqlths[c]-bin.size),q[1])),
                            end = c(p[-1], as.numeric(midpoints[c]), q[-1], as.numeric(seqlths[c])))
        bins$start_gps <- gps(chr=c, position=bins$start)
        bins$end_gps <- gps(chr=c, position=bins$end)
        all.bins <- rbind(all.bins, bins)
    }

    if (is.null(outfile)){
        return(all.bins)
    }else{
        write.table(all.bins, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
    }
}


clusterWithDunn <- function(dmat,
                            maxClusters = Inf,
                            intraclust = c("complete"),
                            interclust = c("average")) {
    
    # Performs hierchical clustering with
    #  optimal number of clusters selected
    #  based on Dunn's Index.
    #
    # Args:
    #   dmat - (n_samples x n_samples) distance matrix
    #   maxClusters - maximum number of clusters to
    #       to use in evaluation.
    #   intraclust - method(s) to use for intracluster
    #       distance estimation
    #   interclust - method(s) to use for intercluster
    #       distance estimation
    # Returns:
    #  list with cluster indices, Dunn's index for each
    #   number of cluster, method used in estimation of
    #   Dunn's index and the number of clusters used
    #   in the optimal configuration.

    # upper bound for number of clusters to evaluate
    nTryCluster <- min(c(maxClusters-1,
                         nrow(dmat)-2))
   
    # obtain hierchical clustering object
    cres <- cluster::agnes(as.dist(dmat))
    # vector to hold value of Dunn's Index
    dvec <- rep(x = 0,
                times = nTryCluster) 

    # vector to hold values for methods used
    mvec <- replicate(n = 2,
                      expr = rep(x = "0",
                                 times = length(dvec)
                      )
                     )

    # evaluate different number of clusters 
    for (k in 1:length(dvec)) {

        # cluster using  k + 1 clusters
        pred <- as.integer(cutree(cres,k+1))
        # evaluate intercluster and intracluster
        # distances

        calcs <- clv::cls.scatt.data(dmat,
                                pred,
                                dist = 'euclidean')

        # compute Dunn's Index from intra/intercluster
        # distances
        dunn <- clv::clv.Dunn(calcs,
                         intracls = intraclust,
                         intercls = interclust)
        
        # get index of best performing clustering
        # method
        maxp <-  which(dunn == max(dunn),
                       arr.ind = TRUE)
        
        # pick first combination if equal
        # value occur
        if (nrow(maxp) > 1) {
            maxp <- maxp[1,]
        }
        
        # store results
        dvec[k] <- as.numeric(dunn[maxp[1],maxp[2]])
        mvec[k,] <- c(rownames(dunn)[maxp[1]],
                      colnames(dunn)[maxp[2]])  
       
    }
    
    # get optimal number of clusters
    maxidx <- which(dvec == max(dvec)) 
    ncl <- as.numeric(maxidx + 1)
    # get cluster labels using optimal number of clusters
    lbls <- as.integer(cutree(cluster::agnes(as.dist(dmat),
                                    method = mvec[maxidx,1]),
                              ncl))
    
    return(list(labels = lbls,
                DunnsI = dvec,
                methodName = mvec,
                nclusters = ncl)
           )
}

makeDistanceBasedHeatMap <- function(dmat,
                                     cluster_labels,
                                     fontsize = 10) {
    
    # Creates a heatmap with rows and cols grouped by
    # provided labels. Requires pheatmap to work.
    #
    # Args:
    #   dmat - (n_samples x n_samples) distance matrix to
    #       visualize.
    #   cluster_labels - (n_samples) vector with a label for
    #       each sample
    # Returns:
    #   pheatmap object
    
    # get order of cluster labels
    ordr <- order(cluster_labels)
    
    # data frame for heatmap annotation
    annot <- data.frame(cluster = as.factor(cluster_labels[ordr]),
                        row.names = rownames(dmat)[ordr])

    # create pheatmap object
    phm <- pheatmap::pheatmap(dmat[ordr,
                                   ordr],
                              cluster_rows = F,
                              cluster_cols = F,
                              fontsize = 10,
                              annotation_row = annot,
                              annotation_col = annot)
    
    return(phm)
}

ReactomePaEnrichment <- function(glist) {

    # Perform enrichment analysis using
    # the Reactome database
    #
    # Args:
    #   glist - vector with genes to be 
    #       queried.
    # Returns:
    #   list with two elements; "tbl" containing a data
    #       frame with enrichment results and "visual"
    #       containing visualizations of the result.
    #       Can be provided as "method" argument
    #       to doDbEnrichment

    # convert symbols to entrez id's
    ezg <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                         keys=glist,
                                         column="ENTREZID",
                                         keytype="SYMBOL",
                                         multiVals="first")
   
    # remove genes where conversion failed
    ezg <- ezg[!(is.na(ezg))]
    # perform enrichment analysis
    res <- try( x <- ReactomePA::enrichPathway(gene=ezg,
                                  pvalueCutoff=0.05,
                                  readable=T)

                    )
   
    # if analysis was successfull store results
    if(!(class(res) == "try-error")) {
        # store visualizations
        viz <- list(emap = ReactomePA::emapplot(x),
                    bplot = graphics::barplot(x,
                                              showCategory = Inf),
                    dotplot = ReactomePA::dotplot(x,
                                                  showCategory=Inf)
                    )

        # convert to writeable format 
        x <- as.data.frame(x)
        # assign rownames as pathway ID
        rownames(x) <- x$ID
        return(list(tbl = x,
                    visual = viz))
    } else {
        return(NULL)
    }
}

doDbEnrichment <- function(glist,
                           label,
                           method) {

    # Perform Database based enrichment of
    # gene sets. 
    #
    # Args:
    #   glist - (n_genes, ) vector with genes to query for
    #       enrichment
    #   label - (n_genes,) vector indicating cluster
    #       membership of each gene
    #   method - function performing an enrichment
    #       analysis. Should return a list with
    #       two elements "tbl" and "visual".
    #       Where tbl holds a data frame of the
    #       enrichment result and visual any
    #       eventual visualizations to be saved.
    #   Returns:
    #       list of two lists; "tbl" containing 
    #       data frames with tabulated
    #       enrichment results and "visual" containing
    #       eventual visualizations.
    
    # get unique labels 
    unilab <- unique(label)
    # for storing results
    rlist <- list(tbl = list(),
                  visual = list()) 

    # iterate over all clusters
    for ( num in 1:length(unilab)) {
        idx <- which(label == unilab[num])
        # perform enrichment
        tmp <- method(glist[idx])
        
        # store tabulated results
        rlist$tbl[[num]] <- tmp$tbl 
        # add cluster identity to result
        rlist$tbl[[num]]$cluster <- as.character(unilab[num])
        # make rownames unique 
        rownames(rlist$tbl[[num]]) <- paste(as.character(unilab[num]),
                                        rownames(rlist$tbl[[num]]),
                                        sep='_')

        # store visualization of results
        rlist$visual[[num]] <- tmp$visual 
    } 

    return(rlist)
}

saveEnrichmentResult <- function(result,
                                 odir,
                                 tag = NULL) {
    
    # Save Enrichment results obtained
    # from doDbEnrichment. Saves tabulated
    # results as tsv files and visual representations
    # as png images.
    #
    # Args:
    #   result - list generated by doDbEnrichment
    #   odir - output directory
    #   tag - unique identifier
   
    # generate joint rownames
    rnames <- unlist(sapply(result$tbl,
                            rownames))
    
    # create joint dataframe of all tabulated results
    result$tbl <- as.data.frame(data.table::rbindlist(result$tbl))
    rownames(result$tbl) <- rnames
    # asseble output filename
    oname <- paste(odir,
                   paste('enrichment_analysis', 
                         tag,
                         'tsv',
                          sep = '.'),
                   sep = '/')

    # save results 
    write.table(result$tbl,
                oname,
                sep = '\t',
                col.names = T,
                row.names = T,
                quote = F
                )

    # save visualizations
    if (!(is.null(result$visual))) {
        for ( v in 1:length(result$visual)) {
            for ( w in 1:length(result$visual[[v]])) {
                imname <- paste(odir,
                                paste('enrichment_visualization',
                                      tag,
                                      names(result$visual[[v]])[w],
                                      as.character(v),
                                      'png',
                                      sep = '.'
                                      ),
                                sep = '/'
                                )

                 png(imname,
                     width = 1000,
                     height  = 2000)

                 print(result$visual[[v]][[w]])
                 dev.off()
             }
        }
    }
}

getClassAccuracy <- function(pred,
                             true,
                             labs) {
acc <- c()
unilabs <- unique(labs)

for (lab in unilabs) {
    idx <- which(labs == lab)
    acc <- c(acc,
             sum(pred[idx] == true[idx]) /
             length(idx)
            )
    }

    names(acc) <- ifelse(unilabs == 1,
                         'tumor',
                         'non')
    return(acc)

}

makeSets <- function(nsamples,
                     ptrain = 0.8,
                     labs = NULL,
                     balance = F) {

    idx <- seq(1:nsamples)

    if (balance) {
        minlen <- min(as.numeric(table(labs)))
    }

    if (is.null(labs)) {

        idx <- sample(idx)
        sep <- floor(length(idx)*0.8)
        train_idx <- idx[1:sep]
        test_idx <- idx[(sep + 1):length(idx)]

    } else {

        train_idx <- c()
        test_idx <- c()

        for (lab in unique(labs)) {
            lidx <- sample(idx[labs == lab])

            if (balance) {
                lidx <- lidx[1:minlen]
            }
            sep <- floor(length(lidx)*0.8)

            train_idx <- c(train_idx,
                           lidx[1:sep]
            )

            test_idx <- c(test_idx,
                          lidx[(sep+1):length(lidx)]
                          )

        }
    }

    return(list(train_idx = train_idx,
                test_idx = test_idx))

}


doLogisiticRegression <- function(dta,
                                  ptrain = 0.8,
                                  tumor_label = 'tumor'
                                  ) {

    y <- dta$meta_data$tumor

    y <- ifelse(y == tumor_label,
                1,
                0)
    
    sets <- makeSets(nrow(dta$count_data),
                          ptrain = 0.8,
                          labs = y,
                          balance =F
                          )
    
    
    cv_lasso <- cv.glmnet(x = dta$count_data[sets$train_idx,],
                          y = y[sets$train_idx],
                          alpha = 1,
                          family = 'binomial'
                          )
    
    
    model <- glmnet(x = dta$count_data[sets$train_idx,],
                    y = y[sets$train_idx],
                    alpha = 1,
                    lambda = cv_lasso$lambda.1se,
                    family = 'binomial'
                    )
    
    betas <- data.frame(gene = colnames(dta$count_data),
                        val = as.numeric(as.numeric(coef(model)))[-1]
                        )
    
    betas <- betas[order(betas$val,decreasing = T),]

    probs <- predict(model,
                     newx = dta$count_data[sets$test_idx,],
                     type = 'response')
    
    
    
    adj_probs <- ifelse(probs > 0.5,
                        1,
                        0
                        )
    
    acc <- getClassAccuracy(adj_probs,
                            y[sets$test_idx],
                            y[sets$test_idx]
                            )

    return(list(betas = betas,
                acc = acc
                )
          )
}
