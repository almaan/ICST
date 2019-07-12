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
	#	break.col: colour of the line breaks seperating consecutive chromosome plots
	#	break.lty: line type of the line breaks seperating consecutive chromosome plots
	#	main.title: Plot title
	#	add.xlab: add x-axis lables
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

# probLogR column order:
# chromosome, pos, logR
# "chr1", 12345, -0.5

plotProbeLogR <- function (probeLogR, color_probe="#e1e1e1"){

	gps_pos <- rep(NA, nrow(probeLogR))
	chr_n <- unique(probeLogR[,1])
	for (n in chr_n) {
		idx <- which(probeLogR[, 1] == n)
		gps_pos[idx] <- gps(chr = n, position=as.numeric(probeLogR[idx, 2]))
	}
	probeLogR <- cbind(probeLogR, gps_pos)

	######## plot ########

	## add bin logR
	points(x=probeLogR[, "gps_pos"], y=probeLogR[,3], pch=".",col=color_probe, cex=0.3)
}

# segLogR column order:
# chromosome, start, end, logR
# "chr1", 12345, 67891, logR

plotSegLogR <- function(segLogR, color_seg="#e60000"){

	gps_start <- rep(NA, nrow(segLogR))
	gps_end <- rep(NA, nrow(segLogR))
	chr_n <- unique(segLogR[, 1])
	for (n in chr_n) {
		idx <- which(segLogR[, 1] == n)
		gps_start[idx] <- gps(chr = n, position=as.numeric(segLogR[idx, 2]))
		gps_end[idx] <- gps(chr = n, position=as.numeric(segLogR[idx, 3]))
	}
	segLogR <- cbind(segLogR, gps_start, gps_end)

	## add seg logR
	segments(x0=segLogR$gps_start, x1=segLogR$gps_end, y0=segLogR[,4], y1=segLogR[,4], col=color_seg, lwd=2)
}
