#!/usr/bin/Rscript

library(RANN)
library(igraph)
library(RColorBrewer)

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
        if ((from %in% nbr[to,]) &  !(to == from)) {
          verts <- sort(c(to,from))
          tov <- c(tov,verts[1])
          frv <- c(frv,verts[2])
          weights <- c(weights,0.01/dist(as.matrix(cnt[c(verts[1],verts[2]),])))
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

  print('here')
  
  result <- list(graph = graph,
                  layout = lo,
                  membership = lov$membership)
  
    
  return(result)
}
  
# Interactive (Temporary)--------


x <- as.numeric(unlist(mt['xcoord']))
y <- as.numeric(unlist(mt['ycoord']))
cnt <- sweep(ct,2,apply(ct,2,sum),'/')

res <- louvain_clustering(x,y,cnt)

nclusters <- length(unique(res$membership))
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

