#!/usr/bin/Rscript

library(RANN)
library(igraph)
library(RColorBrewer)
library(base)

cosine_similarity <- function(xx,yy) {
  nx <- xx / sqrt(xx %*% xx)
  ny <- yy / sqrt(yy %*% yy)
  theta <- acos(nx %*% ny)
  print(theta)
  return(theta)
}

l2_distance <- function(xx,yy) {
  delta <- xx - yy
  delta_norm <- sqrt(delta %*% delta)
  return(delta_norm)
}

louvain_clustering <- function(xx,yy,cnt) {
  
  crd <- data.frame(x= xx, y = yy)
  nbr <- nn2(crd,crd,
             k = 9,
             searchtype = 'radius',
             radius = 1.5)[[1]]
  
  tov <- c()
  frv <- c()
  weights <- c()
  
  for ( from in 1:dim(crd)[1]) {
    for (to in 1:dim(crd)[1]) {
        if ((to %in% nbr[from,]) &  !(to == from)) {
          verts <- sort(c(to,from))
          tov <- c(tov,verts[1])
          frv <- c(frv,verts[2])
          # delta <- 1
          delta <- base::norm(as.matrix(cnt[to,]-cnt[from,]),type = "2")
          weights <- c(weights,delta)
        }
    }
  }
  
  print('ett')
  tedges <- data.frame(from = tov, to = frv)
  rownames(tedges) <- paste('x',c(1:dim(tedges)[1]))
  edges <- unique(tedges)
  weights <- weights[which(rownames(edges) %in% rownames(tedges))] 
  
  graph <- graph_from_data_frame(edges, directed=FALSE, vertices=rownames(crd))
  lov <- cluster_louvain(graph, weights = weights)
  lo <- norm_coords(as.matrix(crd))

  print('here')
  
  result <- list(graph = graph,
                  layout = lo,
                  membership = lov$membership)
  
    
  return(result)
}
  
# Interactive (Temporary)--------


x <- as.numeric(unlist(mt['xcoord']))
y <- as.numeric(unlist(mt['ycoord']))
cnt <- sweep(t(ct),2,apply(t(ct),2,sum),'/')
cnt[is.na(cnt)] <- 0.0

res <- louvain_clustering(x,y,cnt)
print(sort(unique(res$membership)))

nclusters <- length(unique(res$membership))
print(sprintf("found %d clusters", nclusters))
png('/tmp/graph.png')
cmap <- viridis(nclusters * 10)
clr <- cmap[res$membership * 10]


gr <- plot(res[["graph"]],
           layout = res[['layout']],
           vertex.size = 5,
           vertex.label = NA,
           vertex.color = clr
           )

dev.off()

