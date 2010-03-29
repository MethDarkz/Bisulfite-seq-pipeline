infile <- commandArgs(TRUE)[1]
outfile <- commandArgs(TRUE)[2]
width <- as.integer(commandArgs(TRUE)[3])

require(BSgenome.Hsapiens.UCSC.hg18)
temp <- read.csv(infile, header=FALSE, quote="", colClasses=c("character","character","integer","integer","character","integer"))
colnames(temp) <- c("id","chr","fw","rv","strand")

temp$fwSeq=getSeq(Hsapiens, temp$chr, ifelse(temp$strand=="+", temp$fw, temp$rv), width=width, strand=temp$strand)
temp$rvSeq=getSeq(Hsapiens, temp$chr, ifelse(temp$strand=="+", temp$rv, temp$fw), width=width, strand=ifelse(temp$strand=="+", "-", "+"))

write.table(temp, outfile, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")