#!/bin/bash -e
#$1 = fastq filename w/o .fastq.gz
#$2 = Sample name
#$3 = Read Length
#$4 = Fragment length

PIPELINE_PATH=/home/aarsta/analysis-utils/ChIP-seq-pipeline
BOWTIE_PARAMS="-p 7 --nomaqround -n 3 -e 80 --chunkmbs 128 --solexa1.3-quals --best -y -k 2 --strata"
BOWTIE_BS_INDEX="/home/aarsta/data/genome/hg18/bowtie/bis.hg18"
MAX_MM=3
MIN_MM_DIFF=3

#Convert the GenomeData to wiggle and bigwig files
echo "* Converting to wiggle & bigwig";
R --vanilla --slave --quiet --args "$2".Rdata $4 < $PIPELINE_PATH/scripts/GD2wig.R;
wigToBigWig -clip "$2"_UCSC.wig.gz $PIPELINE_PATH/data/hg18.chrom.sizes "$2".bw
gunzip -c "$2"_UCSC.wig.gz | sed -e 's/variableStep/track type=wiggle_0 name="'"$2"'" description="'"$2"'"\
variableStep/' | gzip -c > "$2".wig.gz
rm "$2"_UCSC.wig.gz


