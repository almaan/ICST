source("utils.r")

cpths <- list.files("/Users/hangxu/Projects_scg/ST/projects/BC/data/count/under_tissue", pattern="tsv", full.names=TRUE)
out.dir <- "/Users/hangxu/Projects_scg/ST/projects/BC/data/count_for_cibersort"
f <- cpths[1]
#for(f in cpths){

	r<- read.delim(f, stringsAsFactors=FALSE)

	r[] <- lapply(r, as.numeric) 

	keep.gene <- apply(r, 2, function(x) sum(x)>50)
	keep.spot <- apply(r, 1, function(x) sum(x!=0)>100)

	r <- r[keep.spot, keep.gene]

	r <- t(r)

	r <- cbind(GeneSymbol=rownames(r),r)

	write.table(r, file=file.path(out.dir, basename(f)), sep="\t", quote=FALSE, row.names=FALSE)

#}