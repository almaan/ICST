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
    p2 <- v2 / sum(v2)

    # compute Morisitia-Horn value  
    M <- 2.0 * (p1 %*% p2)
    M <- M / ( sum(p1 ^2 ) + sum(p2 ^2) )
    
    return(M)

    }

getMoransIlarge <- function(fmat,
                          wmat) {

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

    # sum of all weights
    W <- sum(wmat)
    # number of samples 
    N <- nrow(fmat)
    
    # center data
    nfmat <- sweep(fmat,
                   2,
                   colMeans(fmat),
                   '-') 

    # compute denomoniator
    denom <- colSums(nfmat ^ 2 )
    
    # get all combinations of pairs
    combs <- combn(nrow(fmat),
                   m = 2)

    # compute nominator
    nomin <- apply(nfmat,
                   2,
                   function(x) {sum(wmat[cbind(combs[1,],
                                               combs[2,])] * 
                       x[combs[1,]] *
                       x[combs[2,]]) }
                    ) 

    # compute full Moran's I
    I = 2.0 * N / W * nomin / denom  

    return(I)

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
    

    return(apply(df, 2, fun, wmat))
}


getSpatialWeights <- function(crd,
                              sigma = 1) {
    # Get spatially based weights
    #  uses a RBF kernel to generate weights between
    #  points.
    # Args:
    #   crd - (n_samples x 2 ) matrix of coordinates
    #   sigma - (double) scalar giving the length scale
    # Returns:
    #   A symmetric (n_samples x n_samples) matrix
    #   with the spatial weights between each pair
    #   of points
    
    # compute euclidean distance matrix
    wmat <- dist(crd, diag = T)
    # use rbf kernel with distance as similaity measure 
    wmat <- exp(-0.5 * wmat^2 / sigma^2 ) 
    
    return(as.matrix(wmat))

}

    

