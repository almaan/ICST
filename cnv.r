#!/usr/bin/Rscript


source("funcs.r")
source("utils.r")

library(tidyverse)
library(ggplot2)
library(boot)
library(copynumber)
library(cowplot)
library(argparse)


# Parser ----------------------

parser <- ArgumentParser()
parser$add_argument("-c","--count_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-m","--meta_file",
                    type = "character",
                    nargs = "+")

parser$add_argument("-o","--output_dir",
                    type = "character",
                    default = "/tmp/ic10")


## bin_size sets the size of bins for genome
parser$add_argument("-b","--bin_size",
                    type = "integer",
                    default = 1e+7)


args <- parser$parse_args()

########### functions ###########

## Function for CNV heatmap
## plot CNV heatmap, resolution -- 10 kbp
## Params: 
##		seg.data - need columns: sample, start.gps, end.gps, var
##				   sample: sample name/id
##				   start.gps and end gps: absolute gps position on genome which can be calculated using gps function
##				   var is the value of copynumber ratio
##		outfile - cnv heatmap output filename
##		char.   - chromosomes to include in the plot
##		bin.size - resolution of the chromsomes
##		limit - c(limit.low, limit.high) if var is smaller than low.limit, var is set to low.limit, if var is higher than high.limit, var is set to high.limit
##		row.order - order of the rows in heatmap 

cnv.heatmap <- function(seg.data, outfile, chr=paste("chr",c(1:22,"X","Y"),sep=""), bin.size=500000, limit=c(-1.5,1.5), row.order=NULL) {

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

	require(ggplot2)
	require(reshape)
	require(RColorBrewer)

	sample.names <- unique(seg.data[,"sample"])
	sample.num <- length(sample.names)
	dat.mat <- matrix(rep(rep(0,ceiling(XLIM[2]/bin.size)), sample.num), ncol=sample.num)
	colnames(dat.mat) <- sort(sample.names)

	apply(seg.data, 1, function(x){ dat.mat[ round(as.numeric(x["start.gps"])/bin.size) : round(as.numeric(x["end.gps"])/bin.size), x["sample"] ] <<- x["var"]; return(NULL) })

	bin <- seq(1, as.numeric(XLIM[2]), by=bin.size)
	dat.mat <- cbind(dat.mat, bin)
	dat.mat <- as.data.frame(dat.mat)
	dat.melt <- melt(dat.mat, id.vars=c("bin"))
	dat.melt$bin <- as.numeric(as.character(dat.melt$bin))
	dat.melt$value <- as.numeric(as.character(dat.melt$value))

	col.pal <- colorRampPalette(brewer.pal(6,"RdBu"))(100)
	fillcolor <- function(x){
		color <- "white"
		x <- as.numeric(x)
		if (x>limit[2]){
			x <- limit[2]
		}else if (x<limit[1]){
			x <- limit[1]
		}
		idx <- round((x-limit[1])/(limit[2]-limit[1]),2)*100
		if (idx==0){idx=1}
		color <- col.pal[idx]
	}
	col.fill <- apply(dat.melt, 1, function(x) fillcolor(x["value"]))
	dat.melt <- cbind(dat.melt, col.fill)

	#par(mfrow=c(1,1))

	p <- ggplot(dat.melt, aes(x=bin, y=variable, fill=col.fill)) + geom_tile()
	p <- p + scale_fill_identity("", expand=c(0,0))
	if (!is.null(row.order)){
		p <- p + scale_y_discrete(limits=row.order)
	}
	p <- p + scale_x_continuous(breaks=as.numeric(midpoints), labels=as.character(names(midpoints)))
	p <- p + theme_bw()
	p <- p + theme(panel.grid.major=element_blank(), 
				   panel.grid.minor=element_blank(), 
				   panel.border=element_blank(), 
				   panel.background=element_blank(), 
				   axis.title.x=element_blank(), 
				   axis.title.y=element_blank(), 
				   #axis.text.y=element_blank(), 
				   axis.ticks.y=element_blank(), 
				   axis.text.x=element_text(angle=60, vjust=0.5, hjust=0.5))
	#p <- p + geom_hline(yintercept = seq(1:(sample.num-1))+0.5, color="lightgrey")
	p <- p + geom_vline(xintercept=as.numeric(seqlths[-c(1,length(seqlths))]), color="lightgrey", size=0.3)

	ggsave(p, file=outfile, width=12, height=(sample.num/4))

}

#########################

out.dir <- args$output_dir
cpths <- sort(args$count_file)
mpths <- sort(args$meta_file)
genome.bin.size <- args$bin_size

data <- load_multiple(cpths,mpths)

count <- data$count_data
meta <- data$meta_data



## normalize the count data
expmat <- count/rowSums(count)
#expmat <- as.data.frame(t(scale(t(normMat(count, log=TRUE)))), stringsAsFactors=FALSE)


normal.spot <- rownames(meta)[which(meta$tumor=="non")]
normal.expmat <- expmat[normal.spot,]

tmp <- apply(normal.expmat, 2, function(x) c(mean(x),var(x)))
tmp2 <- apply(count[normal.spot,], 2, function(x) sum(x!=0)/length(x))
genes.stat <- cbind(gene=colnames(tmp),as.data.frame(t(tmp)), tmp2)
colnames(genes.stat) <- c("gene","exp.mean","exp.var","perc.non.zero")


## select reference genes
ref.genes <- subset(genes.stat, exp.mean>0 & exp.var<1e-5 & perc.non.zero>0.1)


genes.stat$is.ref <- apply(genes.stat, 1, function(x) x["gene"] %in% ref.genes$gene)

ggplot(genes.stat, aes(x=exp.mean, y=exp.var, size=perc.non.zero, color=is.ref)) + geom_point(alpha=.3)
ggsave(file.path(out.dir,"genes.stat.png"))
write.table(genes.stat, file=file.path(out.dir,"genes.stat.tsv"), row.names=FALSE, sep="\t", quote=FALSE)

mapped <- make_g2chr(genes=ref.genes$gene)

## bin the genome
bins <-build_bins(bin.size=genome.bin.size)

## build ref
## calculate mean expression for each bin from pooled normals

bin.ref <- c()

for (i in 1:nrow(bins)){

	chr <- gsub("chr","",bins[i,1])
  	arm <- as.character(bins[i,2])
  	start <- bins[i,3]
  	end <- bins[i,4]

  	print(sprintf("Calculating ref expression for chr%s%s:%d-%d",chr,arm,start,end))

  	bin.genes <- subset(mapped, chromosome_name==chr & start_position <= end & end_position >= start)
  	num.bin.genes <- nrow(bin.genes)

  	print(sprintf("%d genes within this bin", num.bin.genes))

  	if (num.bin.genes==0){next}

  	ref.bin.expmat <- subset(normal.expmat, select=bin.genes[,1])

  	mean.bin.exp <- function(data, indices){
  		dt <- data[indices,]
  		if (ncol(data)>1){
  			mean(rowMeans(dt))
  		}else{
  			mean(dt)
  		}
  	}

  	bs <- boot(ref.bin.expmat, mean.bin.exp, R=100)
  	bs.bias <- colMeans(bs$t) -bs$t0
  	bs.sd <- apply(bs$t, 2, sd)
  	bs.ci <- boot.ci(boot.out=bs,index=1,type="norm")$norm

  	bin.ref <- rbind(bin.ref, c(chr, arm, start, end, bins[i,5], bins[i,6], num.bin.genes, bs$t0, bs.bias, bs.sd, bs.ci[2], bs.ci[3]))
}

colnames(bin.ref) <- c("chr","arm","start","end","start.gps","end.gps","num.genes","mean.exp","bias","sd","CI95.low","CI95.high")
bin.ref <- as.data.frame(bin.ref, stringsAsFactors=FALSE)
lapply(3:ncol(bin.ref), function(x) bin.ref[,x] <<- as.numeric(bin.ref[,x]))
write.table(bin.ref, file=file.path(out.dir,"reference.tsv"), sep="\t", quote=FALSE, row.names=FALSE)


patients <- unique(meta$patient)


## calculate expression level for each bin, each tumor tissue slide
## calculate bin level logRatio(tumor/normal)
## do segmentation using copynumber pacakge
## plot and write segments for each replicate

for (i in 1:length(patients)){

	pa <- patients[i]
	sections <- unique(meta[which(meta[,"patient"]==pa),"replicate"])

	pa.bin.exp <- c()

	for (j in 1:length(sections)){

		re <- sections[j]

	  	print(sprintf("Patient %s replicate %s", pa, re))

		se.meta <- subset(meta, patient==pa & replicate==re & tumor=="tumor")
		se.expmat <- expmat[rownames(se.meta),]
		se.bin.exp <- c()

		for(b in 1:nrow(bin.ref)){

			chr <- gsub("chr","",bin.ref[b,1])
		  	arm <- as.character(bin.ref[b,2])
		  	start <- bin.ref[b,3]
		  	end <- bin.ref[b,4]

		  	#print(sprintf("Calculating tumor expression for chr%s%s:%d-%d",chr,arm,start,end))

		  	bin.genes <- subset(mapped, chromosome_name==chr & start_position <= end & end_position >= start)
		  	bin.expmat <- subset(se.expmat, select=bin.genes[,1])

		  	cal.fc <- function(data, indices, ref.mean){
		  		fc <- mean(data[indices])/ref.mean
		  		fc
		  	}

		  	if(ncol(bin.expmat)>1){
		  		bin.expmat.data <- rowMeans(bin.expmat)
		  	}else{
		  		bin.expmat.data <- bin.expmat[,1]
		  	}
		  	bs <- boot(bin.expmat.data, cal.fc, R=100, ref.mean=(bin.ref[b,"mean.exp"]+bin.ref[b,"bias"]))
		  	logfc <- log(mean(bs$t[,1]),base=2)
			se.bin.exp <- c(se.bin.exp, logfc)
		}
		pa.bin.exp <- cbind(pa.bin.exp, se.bin.exp)
	}
	colnames(pa.bin.exp) <- sections

	med.pos <- apply(bin.ref, 1, function(x) round((as.numeric(x[3])+as.numeric(x[4]))/2))
	cn.data <- as.data.frame(cbind(Chrom=gsub("X",23,bin.ref$chr), Median.pos=med.pos, pa.bin.exp), stringsAsFactors=F)
	cn.data <- as.data.frame(apply(cn.data, 2, as.numeric), stringsAsFactors=FALSE)
	cn.data[cn.data==-Inf] <- NA

	#logR.mat.win <- winsorize(data=cn.data, pos.unit="bp", arms=NULL, method="mad", tau=2.5, k=25,
    #            assembly="hg19", digits=4, return.outliers=FALSE, save.res=FALSE,
    #            file.names=NULL, verbose=TRUE)
	logR.mat.win <- imputeMissing(cn.data, method="pcf")

	segments <- multipcf(data=logR.mat.win, pos.unit="bp", arms=NULL, Y=NULL, gamma=40, normalize=TRUE,
                          w=1, fast=TRUE, assembly="hg19", digits=4, return.est=FALSE, save.res=FALSE,
                          file.names=NULL, verbose=TRUE)

	segments$start.gps <- apply(segments, 1, function(x) gps(gsub(" ","",x[1]),x[3]))
	segments$end.gps <- apply(segments, 1, function(x) gps(gsub(" ","",x[1]),x[4]))

	pdf(paste0(out.dir,"/", pa,".pdf"))
	par(mfrow=c(3,1))
	nm <- colnames(segments)
	for (i in 6:(ncol(segments)-2)){
		title <- nm[i]
		genomePlot(YLIM=c(-2,2), main.title=title)
		segments(x0=segments$start.gps,x1=segments$end.gps,y0=segments[,i],y1=segments[,i])
	}
	dev.off()

	write.table(segments, file=paste0(out.dir,"/",pa,"_segments.tsv"), sep="\t",quote=FALSE, row.names=FALSE)

}


## prepare input for plotting cnv heatmap 
seg.fns <- list.files(out.dir, pattern="tsv", full.names=TRUE)
segments <- c()

for (f in seg.fns){

	pa <- gsub("_segments.tsv","",basename(f))
	seg <- read.delim(f, stringsAsFactors=FALSE)

	seg.long <- seg %>%
				gather(key=sample, value=var, -c(chrom, arm, start.pos, end.pos, start.gps, end.gps, n.probes)) %>%
				mutate(sample=paste(pa, sample, sep="_"))
	segments <- rbind(segments, seg.long)
}


sample.subtype <- meta %>%
				  select(patient, replicate, subtype) %>%
				  distinct() %>%
				  mutate(sample=paste(patient, replicate, sep="_")) %>%
				  arrange(desc(subtype))

## here is plotting the subtype annotation for each replicate
# p.left <- ggplot(sample.subtype, aes(x=1, y=sample, fill=subtype)) + geom_tile() + scale_y_discrete(limits=row.order) + theme_classic() +
# 			 theme(axis.title.x=element_blank(),
# 			 	   axis.text.x=element_blank(), 
# 				   axis.ticks.x=element_blank(),
# 				   axis.title.y=element_blank(), 
# 				   axis.text.y=element_blank(), 
# 				   axis.ticks.y=element_blank())


## plot cnv heatmap
cnv.heatmap(segments, outfile=file.path(out.dir,"cnv.heatmap.scaled.pdf"), row.order=sample.subtype$sample)






















# END #