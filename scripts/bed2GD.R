infile <- commandArgs(TRUE)[1]
outfile <- commandArgs(TRUE)[2]
bp <- as.integer(commandArgs(TRUE)[3])

temp <- read.csv(gzfile(infile), header=FALSE, row.names=NULL, quote="", colClasses=c("character", "integer", "integer"))
colnames(temp) <- c("chr", "positionFW", "positionRV")

#reverse FW & RV if on minus strand
if (temp[1,2]>temp[1,3]) cols <- 3:2 else cols <- 2:3
alignLocs <- split(temp[,cols], temp$chr)

suppressMessages(require(BSgenome.Hsapiens.UCSC.hg18))
rs <- GenomeData(lapply(alignLocs, function(x) IRanges(start=x[,1], end=x[,2]+bp-1)))

save(rs, file=outfile)
