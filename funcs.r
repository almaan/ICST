#!/usr/bin/Rscript

chr_dist <- function(cnt,
                     nbins = 100,
                     gene_mapping = NULL,
                     edges = NULL) {

    # Map observered transcript to binned chromosome
    # position
    #
    # Args:
    #   cnt - count matrix (n_spots x n_genes)
    #   nbins - number of bins to use
    #   gene_mapping - dataframe with genes as rows,
    #                   chromosomal start/end position
    #                   and chromosome number as columns
    #   edges - n_bins x 2 data frame with lower edges
    #               stored in column "ledge" and upper
    #               edges stored in column "uedge".
    #
    # Return:
    #   list with bin information (upper and lower edge)
    #   and a (n_spots x n_bins) matrix with the
    #   number of observed transcripts in each bin
    #   for each spot

    
    # set gene to chromosome position mapping object
    if (!is.null(gene_mapping)) {

    # if provided mapping object then use this
        G2CHR <- gene_mapping
    } else {

    # otherwise read default mapping file
        G2CHR <- read.table("data/g2chr.tsv",
                        sep = '\t',
                        header = 1,
                        row.names = 1,
                        stringsAsFactors = F)
    }

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
        # generate bins
        ledge <- seq(xmin,xmax, length.out = (nbins + 1))
        uedge <- ledge + diff(ledge)[1]
        
        ledge <- ledge[-length(ledge)]
        uedge <- uedge[-length(uedge)]
    } else {
        # use user-defined edges
        cnames <- c("ledge","uedge")
        hascols <- all(apply(sapply(colnames(edges),grepl,cnames,1,any)))
        
        if (!hascols) {
            print("Provided edges are not compatible")
            return(NULL)
        }

        ledge <- as.numeric(edges$ledge)
        uedge <- as.numeric(edges$uedge)
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

