source("utils.r")
library(mixtools)

#' Receive genomic coordinates of a gene list
#'
#' This function allows to receive the genomic positions of a vector of genes in HUGO format.
#'  param gene_names A vector of gene names in HUGO format.
#'  param ensembl_version Version of the ENSEMBL database used to quantify gene expression data. Defaul: v87.
#'  examples
#' getGenePositions(gene_names=c("EGFR","PDGFRA"))

getGenePositions= function(gene_names,ensembl_version="grch37.ensembl.org",species="human"){
  if (species=="human"){
	ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host=ensembl_version)
	gene_positions <- biomaRt::getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters ='hgnc_symbol', values =gene_names, mart = ensembl)
  }
  else {
	ensembl = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host=ensembl_version)
	gene_positions <- biomaRt::getBM(attributes=c('ensembl_gene_id','mgi_symbol','chromosome_name','start_position','end_position'), filters ='mgi_symbol', values =gene_names, mart = ensembl)
  }
  gene_positions=gene_positions[!duplicated(gene_positions[,2]),]
  gene_positions[which(gene_positions[,3]=="X"),3]=23
  gene_positions[which(gene_positions[,3]=="Y"),3]=24
  gene_positions[which(gene_positions[,3]=="MT"),3]=0
  gene_positions[which(nchar(gene_positions[,3])>2),3]=0
  gene_positions=gene_positions[order(as.numeric(gene_positions[,3]),decreasing=F),]
  return(gene_positions)
}



#' Calculate a normalization factor for each column in an expression matrix
#'
#' This function calculates an adjustment factor for each sample in a expression matrix.

calcNormFactors = function (expmat){
  n=colMeans(expmat)
  return(n)
}





#estimates a two component Gaussian Mixture Model on expression matrix and visualizes the results.


cpths <- list.files("/Users/hangxu/Projects_scg/ST/projects/BC/data/count/under_tissue", pattern="23288", full.names=TRUE)
mpths <- list.files("/Users/hangxu/Projects_scg/ST/projects/BC/data/meta/feature_files", pattern="23288", full.names=TRUE)
data <- load_multiple(cpths,mpths)


gene_pos=getGenePositions(colnames(data$count), ensembl_version<- "www.ensembl.org")

region <- read.delim("data/chromosome_arm_positions_grch38.txt", stringsAsFactors=FALSE)


pdf("ggs.pdf")

par(mfrow=c(3,4))

for (i in 1:nrow(region)){

  chr <- region[i,2]
  start <- region[i,3]
  end <- region[i,4]

  chr_genes=gene_pos[which(gene_pos[,3]==chr & gene_pos[,4]>start &  gene_pos[,5]<end),2]
  if(length(chr_genes)<50){next}

  expmat <- t(data$count)
  normFactor <- calcNormFactors(expmat)

  
  chr_exp=scale(colMeans(expmat[intersect(chr_genes,row.names(expmat)),])-normFactor)


  mixmdl <- normalmixEM(chr_exp, k=2, maxit=1000, maxrestarts=10)

  bestlog <- (-Inf)
  bestmix <- NULL

  if (mixmdl$loglik>bestlog){
          bestlog=mixmdl$loglik
          bestmix=mixmdl
        }

  plot(bestmix,which=2,breaks=50,col1=c("red","green"),main2=paste("Chr: ",chr,":",start,":",end,"\n","Log likelihood ",round(bestmix$loglik,1),sep=""),lwd2=3,xlab2="Expression z-score")

}

dev.off()


# END #