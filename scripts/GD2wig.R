filename <- commandArgs(TRUE)[1]
fragSize <- commandArgs(TRUE)[2]

require(chipseq)
require(Repitools)
load(filename)
rs.total <- laneCounts(rs)

rs.extend <- extendReads(rs, seqLen=fragSize)
rs.coverage <- gdapply(rs.extend, coverage)

require(Repitools)
writeWig(rs.coverage, gsub(".Rdata","_UCSC.wig.gz", filename), normalise=rs.total/1000000)


