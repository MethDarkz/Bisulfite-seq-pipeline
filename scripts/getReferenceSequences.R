infile <- commandArgs(TRUE)[1]
outfile <- commandArgs(TRUE)[2]
width <- as.integer(commandArgs(TRUE)[3])
chunkSize <- 500000

require(BSgenome.Hsapiens.UCSC.hg18)
f1 <- file(infile, open="rt")
fOut <- file(outfile, open="wt")
repeat {
	temp <- as.data.frame(do.call(rbind, strsplit(readLines(f1, n=chunkSize), split=",")), stringsAsFactors=FALSE)
	if (nrow(temp)==0) break
	colnames(temp) <- c("id","chr","fw","rv","strand", "MM")
	temp$fw <- as.integer(temp$fw)
	temp$rv <- as.integer(temp$rv)
	temp$fwSeq=getSeq(Hsapiens, temp$chr, temp$fw, width=width, strand=temp$strand)
	temp$rvSeq=getSeq(Hsapiens, temp$chr, temp$rv, width=width, strand=ifelse(temp$strand=="+", "-", "+"))
	
	writeLines(gsub(" ","",apply(temp, 1, paste, collapse=",")), fOut)
}

close(f1)
close(fOut)
