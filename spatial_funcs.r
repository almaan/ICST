#!/usr/bin/Rscript

getMorisitiaHorn <- function(v1,
                             v2) {
    #Get Morisitia-Horn Co-localization
    # measure for two vectors. As defined in
    # doi.org/10.1186/s13058-015-0638-4
    #
    # Args:
    #   v1 - (n_samples,) with value of feature 1
    #       for each region
    #   v2 - (n_samples,) with value of feature 2
    #       for each region
    # Returns:
    #   Morisitia-Horn colocalization measure
    #   for the two features
    #   

    # get proportions in each region
    p1 <- v1 / sum(v1)
    p1[is.na(p1)] <- 0
    p2 <- v2 / sum(v2)
    p2[is.na(p2)] <- 0

    # compute Morisitia-Horn value  
    M <- 2.0 * (p1 %*% p2)
    M <- M / ( sum(p1 ^2 ) + sum(p2 ^2) )
    
    return(M)

    }


getMoransIlarge <- function(nfmat,
                            wmat
                            ) {

    # Get Moran's I for spatial data
    #   Optimized for analysis with multiple features
    #   Moran's I is defined as
    #    I = N/W * numerator / denominator 
    #    numerator : [ sum_i sum_j w_ij*(x_i - mu)*(x_j-mu)]
    #    denominator :  sum_i (x_i-mu)^2
    # 
    # Args:
    #   fmat - (n_sample x n_features) feature matrix
    #   wmat - (n_sample x n_sample) weights matrix
    # Returns:
    #   Moran's I for each feature
    tstart <- Sys.time()
    # sum of all weights
    W <- sum(wmat)
    # number of samples 
    N <- nrow(nfmat)

    # compute denomoniator
    denom <- colSums(nfmat^2)
    nomin <- future.apply::future_apply(as.matrix(nfmat),
                                        2,
                                        function(x) {
                                            sum(wmat * (x %o% x))
                                                    }
                                        )

    # compute full Moran's I
    I <- (N / W) * (nomin / denom)
    
    # get expected value
    EI <- -1/ (  N - 1) 
   
    # stats for variance compuation
    s1 <- 0.5 * sum( (2*wmat)^2)
    s2 <- sum((2*rowSums(wmat))^2)

    s3 <- ( colSums(nfmat^4) / N ) / ( colSums(nfmat^2) / N ) ^2
    s4 <- (N^2 - 3*N + 3)*s1 - N*s2 + 3*(W^2)
    s5 <- (N^2-N)*s1 - 2*N*s2 + 6*(W^2)

    VarI <- ( N*s4 - s3*s5 ) / ( (N-1) * (N-2) * (N-3) * W^2 ) - EI^2
    SeI <- sqrt(VarI) 

    #z <- (I - EI ) / SeI 
    pvals <- pnorm(I,
                   mean = EI,
                   sd = SeI) 

    pvals <- sapply(pvals,
                    function(x) {
                        ifelse(x > 0,
                               1-x,
                               x)
                    }
                    )
    pvals <- 2*pvals
    
    #pvals <- 2*pnorm(-abs(z))
    padj <- p.adjust(pvals, method = 'BH')

    res <- data.frame(I = I,
                      sd = SeI,
                      p.value = pvals,
                      p.adjust = padj)

    print(Sys.time() - tstart)
    return(res)

    }   

getMoransIsmall <- function(fmat,wmat) {
    # Get Moran's I for spatial data
    # Optimized for analysis with few features
    # Moran's I is defined as
    #    I = N/W * numerator / denominator 
    #    numerator : [ sum_i sum_j w_ij*(x_i - mu)*(x_j-mu)]
    #    denominator :  sum_i (x_i-mu)^2
    #
    # Args:
    #   fmat - (n_sample x n_features) feature matrix
    #   wmat - (n_sample x n_sample) weights matrix
    # Returns:
    #   Moran's I for each feature

    fun <- function(v, w) {
        # center data
        vf <- v - mean(v) 
        # number of smaples
        N <- length(vf)
        # for summation
        nom <- 0.0
        # iterate over all 
        for (ii in 1:(N-1)) {
            for (jj in (ii+1):N) {
                nom <- nom +  vf[ii]*vf[jj]*w[ii,jj]
            }
        }
        # compute denomitator for morans index
        den <- sum(vf^2)
        # compute full Moran's I
        I <- 2 * N/sum(w) * nom / den
        return(I)
    }
    
    tstart <- Sys.time()
    return(apply(df, 2, fun, wmat))
    print(Sys.time()-tstart)
}

getMoransIv2 <- function(nfmat,
                         wmat) {
    tstart <- Sys.time() 
    res <- future.apply::future_apply(nfmat,
                                  2,
                                  ape::Moran.I,
                                  weight = wmat,
                                  alternative = 'two.sided'
                                  )

    res <- as.data.frame(do.call(rbind,
                                 res))
    res[] <- apply(res,
                   2,
                   as.numeric
                   )

    res$p.adjust <- p.adjust(res$p.value,
                             method = 'BH')
    colnames(res)[1] <- "I"
    print(Sys.time() - tstart)
    return(res)

}


getSpatialWeights <- function(crd,
                              startpos = NULL,
                              sigma = 1) {
    # Get spatially based weights
    #  uses a RBF kernel to generate weights between
    #  points. If multiple sections are used weights
    #  between sections are set to zero.
    #
    # Args:
    #   crd - (n_samples x 2 ) matrix of coordinates
    #   startpos - (n_sections) vector. Use if crd
    #       contains coordinates of multiple sections
    #       separated in z-direction. 
    #   sigma - (double) scalar giving the length scale
    #
    # Returns:
    #   A symmetric (n_samples x n_samples) matrix
    #   with the spatial weights between each pair
    #   of points
    
    # if only one section is used
    if (is.null(startpos)) {
        startpos <- c(1,nrow(crd) + 1)
    # if multiple sections are used
    } else {
        startpos <- c(startpos,nrow(crd) + 1)
    }
    # list of section weights
    wlist <- list()
    # compute weights for each section
    # independently
    for (pos in 1:(length(startpos) - 1)) {
        # get section indices
        idx <- c(startpos[pos]:(startpos[pos+1]-1))
        # use rbf kernel with distance as similaity measure 
        wlist[[pos]] <- as.matrix(exp(-0.5 * dist(crd[idx],
                                       diag = T)^2 / 
                                       sigma^2
                           )) 
    }

    # construct joint weight matrix
    wmat <- Matrix::bdiag(wlist) 

    return(as.matrix(wmat))

}

fromMat2List <- function(mat,crd) { 
    longformat <- c()
    for (ii in 1:nrow(mat)) {
       rowv <- sum(mat[ii,])
       if (rowv > 0) {
        longformat <- rbind(longformat,
                            t(replicate(rowv,
                                        crd[ii,]))) 
       }
    }

    return(longformat)
}

clusterByDensity <- function(longformat,
                             eps = 2,
                             minPts = 10) {
    
    clustered <- dbscan::dbscan(longformat,
                                eps = eps,
                                minPts = minPts)
    idx <- !(duplicated(longformat))
    crd <- longformat[idx,]
    lbls <- as.factor(clustered$cluster[idx])

    return(list(labels = lbls, crd = crd ))
}

clusterByExpression <- function(mat,
                                maxClusters=10,
                                criterion = 'BIC'
                             ) {
    print(dim(mat))
    opt_gmm = ClusterR::Optimal_Clusters_GMM(mat,
                                             max_clusters = maxClusters,
                                             criterion = criterion, 
                                             dist_mode = "eucl_dist",
                                             seed_mode = "random_subset",
                                             km_iter = 10,
                                             em_iter = 10,
                                             var_floor = 1e-10, 
                                             plot_data = F)
      
        ncomp <- which(opt_gmm == min(opt_gmm)) 
    
        gmm = ClusterR::GMM(mat,
                            ncomp,
                            dist_mode = "eucl_dist",
                            seed_mode = "random_subset",
                            km_iter = 10,
                            em_iter = 10,
                            verbose = F)        
        
        pr = ClusterR::predict_GMM(mat,
                                   gmm$centroids,
                                   gmm$covariance_matrices,
                                   gmm$weights) 
    
        return(as.factor(pr$cluster_labels))
}    

getCorr <- function(cnt) {
    nc <- ncol(cnt)
    dmat <- matrix(0,nc,nc)
    rownames(dmat) = colnames(cnt)
    colnames(dmat) = colnames(cnt)

    for (xx in 1:(nc-1)) {
        for (yy in (xx +1):nc) {
                dmat[xx,yy] <- 1 - cor(cnt[,xx],
                                       cnt[,yy],
                                       method = 'pearson') 
                dmat[yy,xx] <- dmat[xx,yy]
        }
    }

    return(dmat)
}

getCorr <- function(cnt) {
    nc <- ncol(cnt)
    dmat <- 1 - cor(as.matrix(cnt))
    colnames(dmat)  <- colnames(cnt)
    rownames(dmat) <- colnames(cnt)
    return(dmat)
}

computeGeneralG <- function(cmat,w) {
    # Compute the gradient of the general G
    #  statistic for all samples
    # 
    # Args:
    #   cmat - (n_samples x n_features) count matrix
    #   w - (n_samples x n_samples) weight matrix
    #
    # Returns:
    #   n_features vector of General G statistic
    #   for each feature
    #

    fun <- function(x) {
         # helper function
         # computes general G

         xixj <- (x %*% t(x))

         return(sum(w * xixj ) /
                sum(xixj))
    }

    W <- sum(w)
    n <- nrow(cmat) 
    expval <- W / (n*(n - 1))
    # iterate over all genes
    res <- apply(as.matrix(cmat),
                 2,
                 fun)
    res <- res - expval
    return(res)
}

computeGeneralGgrad <- function(cmat,w) {
    # Compute the gradient of the general G
    #  statistic for all samples
    # 
    # Args:
    #   cmat - (n_samples x n_features) count matrix
    #   w - (n_samples x n_samples) weight matrix
    #
    # Returns:
    #   (n_samples x n_features) matrix of partial
    #   derivatives
    #

    fun <- function(x) {
        # helper function
        # computes gradient w.r.t.
        # one gene

        xixj <- x %*% t(x)
        sm <- sum(x)
        p1 <- sum(xixj)
        p2 <- sum(w * xixj)
        p2 <- p2*(sm + x)
        p3 <- p1^2
        p1 <- p1*(rowSums(sweep(w,1,x,'*')) )
        y <- (p1 - p2) / p3
        y <- ifelse(is.na(y),0,y)
        return(y)

        }
    
    # iterate over all genes
    res <- apply(as.matrix(cmat),
                 2,
                 fun)
    return(res)
}
