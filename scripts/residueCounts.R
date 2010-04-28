infile <- commandArgs(TRUE)[1]
outfile <- commandArgs(TRUE)[2]
bp <- as.integer(commandArgs(TRUE)[3])
chunkSize <- as.integer(commandArgs(TRUE)[4])

temp <- read.csv(gzfile(infile), header=FALSE, row.names=NULL, quote="", colClasses=c("character", "NULL", "integer", "character"))
colnames(temp) <- c("chr", "position", "sequence")
temp$chunk <- paste(temp$chr, trunc(temp$position/chunkSize), sep="-")

temp.split <- split(temp[,1:3], temp$chunk)

out <- gzfile(outfile, open="wt")
thisChr <- "" 

for (chunk in names(temp.split)[order(gsub("-.*","",names(temp.split)), as.integer(gsub(".*-","",names(temp.split))))]) {
	this <- temp.split[[chunk]]
	#have we gone to a new chromosome?
	if (thisChr!=this$chr[1]) {
		thisChr <- this$chr[1]
		cat("\n", thisChr)
		writeLines(thisChr, out)
	}
	from <- min(this$position)
	to <- max(this$position)+bp-1
	tempMatrix <- cbind(pos=from:to, A=0, C=0, G=0, T=0, N=0)

	thisSplit <- do.call(rbind, strsplit(this$sequence, split=""))
	thisPos <- matrix(0:(bp-1), ncol=bp, nrow=nrow(thisSplit), byrow=TRUE) + this$position
	for (res in c("A","C","G","T","N")) {
		tempRes <- table(match(thisPos[thisSplit==res], tempMatrix[,"pos"]))
		tempIndex <- as.integer(names(tempRes))
		tempMatrix[tempIndex,res] <- tempMatrix[tempIndex, res] + tempRes
	}	
	
	tempOut <- which(rowSums(tempMatrix[,2:5])>0)
	writeLines(apply(tempMatrix[tempOut,], 1, paste, collapse=","), out)
	cat(".")
}
cat("\n")
close(out)

