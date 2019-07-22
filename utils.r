#!/usr/bin/Rscript


load_multiple <- function(cpths,
                          mpths) {

    # Construct joint data frame from multiple ST-sections

    # Args:
    #   cpths : vector[character] - vector of full paths to count matrices
    #   mpths : vector[character] - vector of full paths to meta files
    # Returns:
    #  list with the joint count matrix (count_data) and 
    #   joint meta file (meta_data)


    # unlist
    cpths <- unlist(cpths)
    mpths <- unlist(mpths)

    # for storing data
    clist <- list()
    mlist <- list()
    index <- c(0)
    rnames <- c()
    genes <- c()

    for (sample in 1:length(cpths)) {

        # load count data
        ct <- read.table(cpths[sample],
                     sep = '\t',
                     header = 1,
                     row.names = 1)

        # load meta data
        mt <- read.table(mpths[sample],
                         sep = '\t',
                         header = 1,
                          )
        
        # convert entries to character to avoid wrong typecasting 
        mt[] <- lapply(mt,as.character)
       
        # select intersecting set of spots
        inter <- intersect(rownames(ct),rownames(mt)) 
        ct <- ct[inter,]
        mt <- mt[inter,]
        # add unique sample identifer for easy separation
        mt$unique_id <- rep(sample,nrow(mt))     

        # store loaded count files
        clist[[sample]] <- ct
        mlist[[sample]] <- mt
       
        # store number of spots in sample
        index <- c(index,nrow(ct))
        # add unique rownames for joint matrix
        rnames <- c(rnames,paste(sample,rownames(ct),sep = '_'))
        # modify set of observed genes
        genes <- union(genes,colnames(ct))

    }
    # compute start positions for each sample
    index <- cumsum(index)
    # prepare  empty count data frame to fill
    cjoint <- as.data.frame(matrix(0,nrow = index[length(index)], ncol = length(genes))) 
    rownames(cjoint) <- rnames
    colnames(cjoint) <- genes
    # prepare empty meta data frame to fill 
    mjoint <- as.data.frame(matrix(0,nrow = index[length(index)], ncol = dim(mlist[[1]])[2]))
    rownames(mjoint) <- rnames 
    colnames(mjoint) <- colnames(mlist[[1]])

    # fill joint data frames
    for ( pos in 1:(length(index)-1)) {
        cjoint[(index[pos]+1):(index[pos+1]),colnames(clist[[pos]])] <- clist[[pos]]  
        mjoint[(index[pos]+1):(index[pos+1]),] <- mlist[[pos]]
    }
  
    # list of joined marices
    joint <- list(count_data = cjoint,
                  meta_data = mjoint,
                  cumsum = cumsum)

    return(joint)

}



gps <- function(chr="chr3", position=c(12000000, 8000000)) {

    # Convert a genomic position to a (x-axis) plot position
    #
    # Args:
    #   chr: Chromosome (UCSC notation)
    #   position: Chromosomal position (hg19)
    #
    # Returns:
    #   x-axis position for genomic plots.

    require(GenomicRanges)
    # Reference data
    data(hg19IdeogramCyto, package = "biovizBase")
    data(hg19Ideogram, package = "biovizBase")
    chrs <- as.character(levels(seqnames(hg19IdeogramCyto)))
    seqlths <- seqlengths(hg19Ideogram)[chrs]
    new.order <- c("chr1", "chr2", "chr3", "chr4", "chr5",
                   "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11","chr12","chr13","chr14","chr15",
                   "chr16","chr17","chr18","chr19","chr20",
                   "chr21","chr22","chrX","chrY")
    seqlths <- seqlths[new.order]

    # Error handling
    if (chr %in% 1:22){
        chr <- paste("chr",chr,sep="")
    }else if (chr==23 || chr=="X"){
        chr <- "chrX"
    }else if (chr==24 || chr=="Y"){
        chr <- "chrY"
    }

    test <- which(new.order == chr)
    if (length(test) == 0) {
        stop("Chromosome not found: chr1-22/X/Y expected.")
    }
    position <- as.numeric(position)

    # Translation stage
    if(chr != "chr1") {
        hit = which(names(seqlths) == chr) - 1
        hit = sum(as.numeric(seqlths[1:hit])) + as.numeric(position)
    }else{
        hit = as.numeric(position)
    }
    return(hit)
} # end of gps function


genomePlot <- function(chr=paste("chr",c(1:22,"X","Y"),sep=""), YLIM = c(-1.5, 1.5), break.col="lightgrey", break.lty=1, main.title="TITLE", add.xlab=TRUE, YLAB=NA) {

    # Plot an empty log-R genome plot
    #
    # Args:
    #   chr: Chromosomes to plot.
    #   YLIM: y-axis range
    #   break.col: colour of the line breaks seperating consecutive chromosome plots
    #   break.lty: line type of the line breaks seperating consecutive chromosome plots
    #   main.title: Plot title
    #   add.xlab: add x-axis lables
    #
    # Returns:
    #   An empty plot.

    require(GenomicRanges)
    # Reference data
    data(hg19IdeogramCyto, package = "biovizBase")
    data(hg19Ideogram, package = "biovizBase")
    chrs <- as.character(levels(seqnames(hg19IdeogramCyto)))
    seqlths <- seqlengths(hg19Ideogram)[chrs]
    new.order <- c("chr1", "chr2", "chr3", "chr4", "chr5",
                   "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11","chr12","chr13","chr14","chr15",
                   "chr16","chr17","chr18","chr19","chr20",
                   "chr21","chr22","chrX", "chrY")
    seqlths <- seqlths[new.order]
    midpoints <- seqlths/2

    # Error handling
    test <- diff(match(chr, new.order))
    test <- which(test < 0)
    if (length(test) > 0) {
        stop("Please re-order chromosomes in your input vector (by common-sense ordering)")
    }

    # convert distances using gps()
    seqlths <- seqlths
    for (i in 1:length(seqlths)) {
        seqlths[i] <- gps(names(seqlths)[i], position=seqlths[i])
        midpoints[i] <- gps(names(midpoints)[i], position=midpoints[i])
    }


    # Define x-axis limits
    flags <- which(names(seqlths) %in% chr)
    flags <- range(flags)
    if (flags[1] == 1) {
        XLIM <- c(0, seqlths[flags[2]])
        seqlths <- c(0, seqlths[which(names(seqlths) %in% chr)])
    }else{
        XLIM <- c(seqlths[flags[1]-1], seqlths[flags[2]])
        seqlths <- c(seqlths[flags[1]-1], seqlths[which(names(seqlths) %in% chr)])
    }

    # create empty plot
    plot(x=0, type="n", xlim=XLIM,
        ylim=YLIM , xaxt="n", xlab=NA, main=main.title, frame.plot=FALSE, ylab=YLAB)
    # Add inter-chromosomal breaks
    abline(v=seqlths, col=break.col, lty=break.lty)
    # Add chromosome lables
    if (add.xlab) {
    mtext(text = gsub("chr", "", names(midpoints)[1:length(chr)]),
          side=1, at=midpoints[1:length(chr)], cex=0.6)
    }

   return(NULL);
}



# Filter and expression matrix
# This function filters an expression matrix by removing genes which are expressed at low levels and in only very few spots
# Args:
#   mat - spot (row) x genes (column) expression matrix
#   genelist - if not NULL, genes not in the genelist will be removed from the matrix
#   minExp - genes with colSums > minExp are being kept
#   minSpot - genes expressed in > minSpot spots will be kept
# Return:
#   filtered matrix

filterMatrix <- function(mat,genelist=NULL,minExp=50,minSpot=10,minNonZero=10){

    if (!is.null(genelist)){
        mat <- mat[intersect(genelist,colnames(mat)),]
    }
    mat <- mat[,which(colSums(mat)>minExp)]
    percentgenes=apply(mat,2,function(x){sum(x>0)})
    mat=mat[,percentgenes>minSpot]

    keep <- apply(mat, 2, function(x) sum(x!=0)>minNonZero)
    mat <- mat[,keep]
    return(mat)
}


# Normalize a expression matrix of raw counts to log2 (CPM/10+1) values
#
# This function normlizes a matrix of raw gene counts to log2(CPM/10+1).
# Args:
#    expmat - A spot X genes expression matrix of raw RNA counts
# Return:
#    normalized expression matrix

normMat = function (expmat, log=FALSE){
  sdepth=rowSums(expmat)
  rmat=expmat / sdepth
  if (log){
    rmat=rmat*1000000
    rmat=log2(rmat/10+1)
  }
  return(rmat)
}


calcNormFactors <- function (expmat){

    # This function calculates an adjustment factor for each spot in a normalized expression matrix.
    # Args:
    #   expmat - each row represnt a spot, each col represent a RNA feature
    # Returns:
    #   a vector of an adustment factor for each spot

    n=rowMeans(expmat)
    return(n)
}


getCoordinates <- function(spotnames,
                           delim = 'x') {

    # Get coordinates from a list of spot names
    #
    # Args:
    #   spotnames - list with name of spots
    #   delim - character separating x and y cooridnate
    #       in spotnames
    # Returns:
    #   (n_spots x 2 ) - matrix of coordinates
    
    # helper function for split
    splstr <- function(x,p) { strsplit(x,'x')[[1]][p]}
    
    # get coordinates
    xcrd <- as.numeric(unlist(sapply(spotnames,splstr,1)))
    ycrd <- as.numeric(unlist(sapply(spotnames,splstr,2)))
    # join coordinates 
    crd <- cbind(xcrd, ycrd)

    return(crd)
}
    



