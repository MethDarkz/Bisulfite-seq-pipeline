infile <- commandArgs(TRUE)[1]
outfile <- commandArgs(TRUE)[2]

temp <- read.csv(gzfile(infile), header=FALSE, row.names=NULL, quote="", colClasses=c("NULL", "character", "integer", "character"))
colnames(temp) <- c("chr", "position", "strand")

alignLocs <- split(temp[,2:3], temp$chr)

require(BSgenome.Hsapiens.UCSC.hg18)
rs <- GenomeData(lapply(alignLocs, function(df) with(df, list(
	"+"=df$position[df$strand=="+"],
	"-"=df$position[df$strand=="-"]))))

save(rs, file=outfile)
