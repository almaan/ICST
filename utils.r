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
    joint <- list(count_data = cjoint, meta_data = mjoint)

    return(joint)

}

