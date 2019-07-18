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
    escore[is.na(escore)] <- 0.0

    return(escore)
    
}
