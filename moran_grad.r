#!/usr/bin/Rscript

x <- c(1,2,3,4,5,6)
xm <- mean(x)
N <- length(x)
dx <- rep(0,N)


getMoranGrad <- function(x,w) {
    N <- length(x)
    xm <- mean(x)
    g <- sum((x- xm) ^ 2 )
    dx <- rep(0,N)
    
    for (i in 1:length(x)) {
        xi <- 0.0
        gp <- -2*(sum(x-xm) / N ) + 2*(x[i] - xm) 
        for (k in 1:length(x)) {
            f <- 0.0
            fp <- 0.0
            for (j in 1:length(x)) {
                f <- f + w[k,j] * (x[j]-xm)*(x[k] - xm) 
                fp <- fp -x[k]/N - x[j] / N + 2 / N * xm
    
                if (i == j | i == k ) {
                    fp = fp - x[i] / N - xm
                }
    
                fp <- fp +  w[k,j] * fp
                dx[i] <- dx[i] + (fp*g - f*gp ) / g^2
            
            }
        }
    }
    print('done')
    return(dx)
}
