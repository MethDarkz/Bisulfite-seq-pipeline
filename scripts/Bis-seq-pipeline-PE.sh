#!/bin/bash -e 

#Check for existence of all necessary files and tools
if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  exit 1;
fi

source "$1";

if [ ! -e "$FASTQ1" ]; then
  echo "$FASTQ1"" does not exist!";
  exit 1;
fi

if [ ! -e "$FASTQ2" ]; then
  echo "$FASTQ2"" does not exist!";
  exit 1;
fi

if [ ! -e "$PIPELINE_PATH"/sql/Bis-seq-PE.schema ]; then
  echo "$PIPELINE_PATH""/sql/Bis-seq-PE.schema does not exist!";
  exit 1;
fi

#init the database
echo `date`" - Initialising the database";
rm -f "$PROJECT".db;
sqlite3 "$PROJECT".db < "$PIPELINE_PATH"/sql/Bis-seq-PE.schema;

#Convert fastq into csv
echo `date`" - Importing reads into the database";
#Import FW & RV seqs in db, one row per read id
gunzip -c "$FASTQ1" |awk -v file2="$FASTQ2" 'BEGIN {
  FS="/"
  cmd="gunzip -c " file2
} {
  readname=$1
  getline FWseq
  getline
  getline FWqual

  cmd | getline
  cmd | getline RVseq
  cmd | getline
  cmd | getline RVqual
  
  tmp1=$1
  tmp2=$3
  tmp3=$4
  getline < file2
  print readname "|" FWseq "|" RVseq "|" FWqual "|" RVqual
}' | sed 's/^@//' | sort -t "|" -k1,1 | sqlite3 "$PROJECT".db '.import /dev/stdin reads';

#Convert C residues in reads to T
echo `date`" - Bisulfite converting reads";
#setup named pipes
rm -f fastq1 fastq2;
mkfifo fastq1 fastq2;
gunzip -c "$FASTQ1" | sed -e 's/ .*//' -e '2~4s/C/T/g' > "$PROJECT".conv.fastq1;
gunzip -c "$FASTQ2" | sed -e 's/ .*//' -e '2~4s/G/A/g' > "$PROJECT".conv.fastq2;

#Map against forward strand, and filter out reads with too many MMs
echo `date`" - Bowtie mapping against forward strand";
"$BOWTIE_PATH"/bowtie --norc "$BOWTIE_PARAMS" "$GENOME_PATH"/plus -1 "$PROJECT".conv.fastq1 -2 "$PROJECT".conv.fastq2 2> mapping.plus.log |  sed -e 'N' -e 's/\n/\t/' | awk -v maxmm=$(($MAX_MM+$MIN_MM_DIFF)) 'BEGIN {FS="\t"}{
  numFW=split($8, tmpFW, ":")-1
  if (numFW==-1) numFW=0
  numRV=split($16, tmpRV, ":")-1
  if (numRV==-1) numRV=0
  if ((numFW+numRV)<maxmm) {
    print substr($1,1,length($1)-2) "|" $3 "|" ($4+1) "|" ($12+1) "|+|" (numFW+numRV)
  }
}' | sed 's/plusU_//'> "$PROJECT".both.map;

#re-setup named pipes for mapping against reverse strand
gunzip -c "$FASTQ1" | sed -e 's/ .*//' -e '2~4s/C/T/g' > fastq1 &
gunzip -c "$FASTQ2" | sed -e 's/ .*//' -e '2~4s/G/A/g' > fastq2 &

#Same for reverse strand
echo `date`" - Bowtie mapping against reverse strand";
"$BOWTIE_PATH"/bowtie --nofw "$BOWTIE_PARAMS" "$GENOME_PATH"/minus -1 "$PROJECT".conv.fastq1 -2 "$PROJECT".conv.fastq2 2> mapping.minus.log | sed -e 'N' -e 's/\n/\t/' | awk -v maxmm=$(($MAX_MM+$MIN_MM_DIFF)) 'BEGIN {FS="\t"}{
  numFW=split($8, tmpFW, ":")-1
  if (numFW==-1) numFW=0
  numRV=split($16, tmpRV, ":")-1
  if (numRV==-1) numRV=0
  if ((numFW+numRV)<maxmm) {
    print substr($1,1,length($1)-2) "|" $3 "|" ($12+1) "|" ($4+1) "|-|" (numFW+numRV)
  }
}' | sed 's/minusU_//'>> "$PROJECT".both.map;

gzip -f "$PROJECT".conv.fastq1;
gzip -f "$PROJECT".conv.fastq2;

#adjust no of mismatches for C's in read that are T's in the reference
echo `date`" - Getting the reference sequence of reads mapping positions"
sort -t "|" -k2,2 -k3,3n "$PROJECT".both.map | awk -v readLength="$READ_LENGTH" -v pipeline="$PIPELINE_PATH" -v genomePath="$GENOME_PATH" 'BEGIN {FS="|"}
function revComp(temp) {
    for(i=length(temp);i>=1;i--) {
        tempChar = substr(temp,i,1)
        if (tempChar=="A") {printf("T")} else
        if (tempChar=="C") {printf("G")} else
        if (tempChar=="G") {printf("C")} else
        if (tempChar=="T") {printf("A")} else
        {printf("N")}
    }
} { #entry point
    if ($2!=chr) { #if hit a new chromosome, read it into chrSeq
        chr=$2
        "awk -f "pipeline"/scripts/readChr.awk "genomePath"/"chr".fa" | getline chrSeq
    }
    FW=toupper(substr(chrSeq,$3,readLength)) #retrieve forward sequence
    RV=toupper(substr(chrSeq,$4,readLength)) #retrieve reverse sequence
    printf("%s|%s|%s|%s|%s|%s|",$1,$2,$3,$4,$5,$6)
    if ($5=="+") {
        printf("%s|",FW)
        revComp(RV)
        printf("\n")
    } else {
        revComp(FW)
        printf("|%s\n",RV)
    }
}' | sort -t "|" -k1,1 |  sqlite3 "$PROJECT".db '.import /dev/stdin mappingBoth'
gzip -f "$PROJECT".both.map;

echo `date`" - Adjusting number of mismatches for C->T errors in mapping"
sqlite3 "$PROJECT".db "SELECT mappingBoth.*, reads.sequenceFW, reads.sequenceRV
  FROM mappingBoth JOIN reads ON mappingBoth.id=reads.id;" | awk -v bp="$READ_LENGTH" 'BEGIN {FS="|"} {
	mm=0;
    for(i=1;i<=bp;i++) {
        	temp1 = substr($7,i,1)
    	temp2 = substr($9,i,1)
    	if (temp1=="C") {
    		if (temp2!="C"&&temp2!="T") mm++;
    	} else if (temp1!=temp2) mm++;
    	temp1 = substr($8,i,1)
    	temp2 = substr($10,i,1)
    	if (temp1=="G") {
    		if (temp2!="G"&&temp2!="A") mm++;
    	} else if (temp1!=temp2) mm++;
    }
    print $1"|"$2"|"$3"|"$4"|"$5"|"mm"|"$7"|"$8
}' | gzip -c > "$PROJECT".both.adjust.csv.gz;
gunzip -c "$PROJECT".both.adjust.csv.gz | sort -t "|" -k1,1 | sqlite3 "$PROJECT".db '.import /dev/stdin mappingAdjust'

echo `date`" - Combining  forward and reverse mappings and filtering out duplicate mappings";
gunzip -c "$PROJECT".both.adjust.csv.gz | sort -t "|" -k 1,1 -k 6,6n  | awk -v maxmm=$MAX_MM -v mindiff=$MIN_MM_DIFF 'BEGIN {FS = "|"} {
    s = $1
    if (s != prevs) {
        if ( FNR > 1 ) {
			if ( prevval<=maxmm && valdiff>=mindiff) print prevline
        }
        prevval = $6
        prevline = $0
        valdiff = mindiff
    }
    else valdiff = prevval-$6
    prevs = s
}
END {
	if ( prevval<=maxmm && valdiff>=mindiff) print prevline
}' | sort -t "|" -k1,1 |  sqlite3 "$PROJECT".db '.import /dev/stdin mapping'

echo `date`" - Exporting database to BAM files"
cp "$GENOME_PATH"/samdir "$PROJECT".mappings.plus.sam;
sqlite3 "$PROJECT".db "SELECT
mapping.id, mapping.chr, mapping.strand, mapping.positionFW, reads.sequenceFW, mapping.positionRV, reads.sequenceRV, reads.qualityFW, reads.qualityRV
FROM mapping LEFT JOIN reads ON mapping.id=reads.id WHERE mapping.strand='+';" | awk -v readLength="$READ_LENGTH" 'BEGIN {FS = "|"} 
	function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
	function revComp(temp) {
		for(i=length(temp);i>=1;i--) {
			tempChar = substr(temp,i,1)
			if (tempChar=="A") {printf("T")} else
			if (tempChar=="C") {printf("G")} else
			if (tempChar=="G") {printf("C")} else
			if (tempChar=="T") {printf("A")} else
			{printf("N")}
		}
		printf("\t")
	}
	function rev(temp) {
		for(i=length(temp);i>=1;i--) printf(substr(temp,i,1))
		printf("\n")
	} {
	print $1 "\t99\t" $2 "\t" $4 "\t255\t"readLength"M\t=\t" $6 "\t" (abs($6-$4)+readLength) "\t" $5 "\t" $8
	printf("%s\t147\t%s\t%s\t255\t"readLength"M\t=\t%s\t%s\t",$1,$2,$6,$4,(abs($6-$4)+readLength))
	revComp($7)
	rev($9)
}' >> "$PROJECT".mappings.plus.sam;
"$SAMTOOLS_PATH"/samtools import "$GENOME_PATH"/reflist "$PROJECT".mappings.plus.sam "$PROJECT".mappings.plus.bam;
"$SAMTOOLS_PATH"/samtools sort "$PROJECT".mappings.plus.bam "$PROJECT".mappings.plus.sort;
"$SAMTOOLS_PATH"/samtools index "$PROJECT".mappings.plus.sort.bam;
rm "$PROJECT".mappings.plus.bam "$PROJECT".mappings.plus.sam;

cp "$GENOME_PATH"/samdir "$PROJECT".mappings.minus.sam;
sqlite3 "$PROJECT".db "SELECT
mapping.id, mapping.chr, mapping.strand, mapping.positionFW, reads.sequenceFW, mapping.positionRV, reads.sequenceRV, reads.qualityFW, reads.qualityRV
FROM mapping LEFT JOIN reads ON mapping.id=reads.id WHERE mapping.strand='-';" | awk -v readLength="$READ_LENGTH" 'BEGIN {FS = "|"} 
	function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
	function revComp(temp) {
		for(i=length(temp);i>=1;i--) {
			tempChar = substr(temp,i,1)
			if (tempChar=="A") {printf("T")} else
			if (tempChar=="C") {printf("G")} else
			if (tempChar=="G") {printf("C")} else
			if (tempChar=="T") {printf("A")} else
			{printf("N")}
		}
		printf("\t")
	}
	function rev(temp) {
		for(i=length(temp);i>=1;i--) printf(substr(temp,i,1))
		printf("\n")
	} {
	printf("%s\t83\t%s\t%s\t255\t"readLength"M\t=\t%s\t%s\t",$1,$2,$4,$6,(abs($6-$4)+readLength))
	revComp($5)
	rev($8)
	print $1 "\t163\t" $2 "\t" $6 "\t255\t"readLength"M\t=\t" $4 "\t" (abs($6-$4)+readLength) "\t" $7 "\t" $9
}' >> "$PROJECT".mappings.minus.sam;
"$SAMTOOLS_PATH"/samtools import "$GENOME_PATH"/reflist "$PROJECT".mappings.minus.sam "$PROJECT".mappings.minus.bam;
"$SAMTOOLS_PATH"/samtools sort "$PROJECT".mappings.minus.bam "$PROJECT".mappings.minus.sort;
"$SAMTOOLS_PATH"/samtools index "$PROJECT".mappings.minus.sort.bam;
rm "$PROJECT".mappings.minus.bam "$PROJECT".mappings.minus.sam;

#create bed & Rdata (GD) file for coverage mapping for each strand
echo `date`" - Creating coverage bed and GenomeData files";
sqlite3 -csv "$PROJECT".db "SELECT chr, positionFW, positionRV FROM mapping WHERE strand='+';" | gzip -c > "$PROJECT".plus.bed.gz;
sqlite3 -csv "$PROJECT".db "SELECT chr, positionFW, positionRV FROM mapping WHERE strand='-';" | gzip -c > "$PROJECT".minus.bed.gz;

"$R_PATH"/R --vanilla --slave --args "$PROJECT".plus.bed.gz "$PROJECT".plus.Rdata "$READ_LENGTH" < "$PIPELINE_PATH"/scripts/bed2GD.R;
"$R_PATH"/R --vanilla --slave --args "$PROJECT".minus.bed.gz "$PROJECT".minus.Rdata "$READ_LENGTH" < "$PIPELINE_PATH"/scripts/bed2GD.R;

#Are C's found in CpG sites?
echo `date`" - Determining context of C residues"
sqlite3 -csv "$PROJECT".db "SELECT mapping.chr, mapping.positionFW, mapping.positionRV, mapping.strand, reads.sequenceFW, reads.sequenceRV
  FROM mapping JOIN reads ON mapping.id=reads.id;" | sort -t "," -k1,1 | awk 'BEGIN {FS = ","} {
	num=split($5, tmp, "C"); 
	upto=0
	for (i = 2; i <= num; i++) {
		upto = upto + 1 + length(tmp[i-1]);
		if ($4=="+") print $1 "," ($2+upto-1)
		else print $1 "," ($2+length($5)-(upto)-1)
	}
	num=split($6, tmp, "G"); 
	upto=0
	for (i = 2; i <= num; i++) {
		upto = upto + 1 + length(tmp[i-1]);
		if ($4=="+") print $1 "," ($3+length($5)-(upto))			
		else print $1 "," ($3+upto-2)
	}
}' | awk -v pipeline="$PIPELINE_PATH" -v genomePath="$GENOME_PATH" 'BEGIN {FS=","}
function revComp(temp) {
    for(i=length(temp);i>=1;i--) {
        tempChar = substr(temp,i,1)
        if (tempChar=="A") {printf("T")} else
        if (tempChar=="C") {printf("G")} else
        if (tempChar=="G") {printf("C")} else
        if (tempChar=="T") {printf("A")} else
        {printf("N")}
    }
} { #entry point
    if ($1!=chr) { #if hit a new chromosome, read it into chrSeq
        chr=$1
        "awk -f "pipeline"/scripts/readChr.awk "genomePath"/"chr".fa" | getline chrSeq
    }
    print(toupper(substr(chrSeq,$2,2)))
}' | sort | uniq -c > "$PROJECT".context;

#How many Cs in read vs reference?
echo `date`" - Determining conversion %"
sqlite3 -csv "$PROJECT".db "SELECT reads.sequenceFW, reads.sequenceRV, mapping.referenceFW, mapping.referenceRV
  FROM mapping JOIN reads ON mapping.id=reads.id;" | awk 'BEGIN {
	FS = ","
	readC=0
	refC=0
} {
	readC=readC+split($1,tmp,"C")-1
	readC=readC+split($2,tmp,"G")-1
	refC=refC+split($3,tmp,"C")-1
	refC=refC+split($4,tmp,"G")-1
} END {
	print readC","refC
}' > "$PROJECT".conversion;

NUM_READS=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM reads;"`;
NUM_UNIQUE=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM mapping;"`;
NUM_UNMAPPABLE=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM reads \
WHERE id NOT IN(SELECT id FROM mappingBoth);"`;
NUM_MULTIPLE=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM reads \
WHERE id IN(SELECT id FROM mappingBoth) \
AND id NOT IN (SELECT id FROM mapping);"`;

#Creating mapping log
echo `date`" - Creating mapping log"
echo "# reads processed: $NUM_READS" > mapping.log;
echo "# reads with at least one reported alignment: $NUM_UNIQUE" >> mapping.log;
echo "# reads that failed to align: $NUM_UNMAPPABLE" >> mapping.log;
echo "# reads with alignments suppressed due to -m: $NUM_MULTIPLE" >> mapping.log;
echo "Reported $NUM_UNIQUE alignments to 1 output stream(s)" >> mapping.log;
