#!/bin/bash -e

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

#init the database
echo `date`" - Initialising the database";
rm -f "$PROJECT".db;
sqlite3 "$PROJECT".db < "$PIPELINE_PATH"/sql/Bis-seq-PE.schema;

#Convert fastq into csv
echo `date`" - Importing reads into the database";
gunzip -c "$FASTQ1"  | awk 'NR%4==1 || NR%4==2 || NR%4==0' | sed -e 's/ .*//' | sed -e 'N' -e 'N' -e 's/\//,/' -e 's/\n/,/g' -e 's/^@//' >"$PROJECT".fastq1.csv;
gunzip -c "$FASTQ2"  | awk 'NR%4==1 || NR%4==2 || NR%4==0' | sed -e 's/ .*//' | sed -e 'N' -e 'N' -e 's/\//,/' -e 's/\n/,/g' -e 's/^@//' >"$PROJECT".fastq2.csv;
#Import FW & RV seqs in db, one row per read id
awk -v file2="$PROJECT".fastq2.csv 'BEGIN {FS=","} {
  tmp1=$1
  tmp2=$3
  tmp3=$4
  getline < file2
  print tmp1","tmp2","$3","tmp3","$4}' "$PROJECT".fastq1.csv | gzip -c > "$PROJECT".reads.csv.gz
gzip -f "$PROJECT".fastq1.csv;
gzip -f "$PROJECT".fastq2.csv;
gunzip -c "$PROJECT".reads.csv.gz | sort -t "," -k1,1 | sed -e 's/,/|/g' | sqlite3 "$PROJECT".db '.import /dev/stdin reads'


#Convert C residues in reads to T
echo `date`" - Bisulfite converting reads";
gunzip -c "$FASTQ1" | sed -e 's/ .*//' -e '2~4s/C/T/g' > "$PROJECT".conv.fastq1;
gunzip -c "$FASTQ2" | sed -e 's/ .*//' -e '2~4s/G/A/g' > "$PROJECT".conv.fastq2;

#Map against forward strand, and filter out reads with too many MMs
echo `date`" - Bowtie mapping against forward strand";
bowtie --norc "$BOWTIE_PARAMS" "$BOWTIE_BS_INDEX".plus -1 "$PROJECT".conv.fastq1 -2 "$PROJECT".conv.fastq2 2> mapping.plus.log | gzip -c > "$PROJECT".plus.map.gz;

gunzip -c "$PROJECT".plus.map.gz | sed -e 'N' -e 's/\n/\t/' | awk -v maxmm=$(($MAX_MM+$MIN_MM_DIFF)) 'BEGIN {FS="\t"}{
  numFW=split($8, tmpFW, ":")-1
  if (numFW==-1) numFW=0
  numRV=split($16, tmpRV, ":")-1
  if (numRV==-1) numRV=0
  if ((numFW+numRV)<maxmm) {
    print substr($1,1,length($1)-2) "," $3 "," ($4+1) "," ($12+1) ",+," (numFW+numRV)
  }
}' | sed 's/plusU_//'> "$PROJECT".both.map.csv

#Same for reverse strand
echo `date`" - Bowtie mapping against reverse strand";
bowtie --nofw "$BOWTIE_PARAMS" "$BOWTIE_BS_INDEX".minus -1 "$PROJECT".conv.fastq1 -2 "$PROJECT".conv.fastq2 2> mapping.minus.log | gzip -c > "$PROJECT".minus.map.gz;

gunzip -c "$PROJECT".minus.map.gz | sed -e 'N' -e 's/\n/\t/' | awk -v maxmm=$(($MAX_MM+$MIN_MM_DIFF)) 'BEGIN {FS="\t"}{
  numFW=split($8, tmpFW, ":")-1
  if (numFW==-1) numFW=0
  numRV=split($16, tmpRV, ":")-1
  if (numRV==-1) numRV=0
  if ((numFW+numRV)<maxmm) {
    print substr($1,1,length($1)-2) "," $3 "," ($12+1) "," ($4+1) ",-," (numFW+numRV)
  }
}' | sed 's/minusU_//'>> "$PROJECT".both.map.csv

gzip -f "$PROJECT".conv.fastq1;
gzip -f "$PROJECT".conv.fastq2;

#adjust no of mismatches for C's in read that are T's in the reference
echo `date`" - Getting the reference sequence of reads mapping positions"
sort --field-separator="," -k2,2 -k3,3n "$PROJECT".both.map.csv | awk -v readLength="$READ_LENGTH" -v pipeline="$PIPELINE_PATH" -v genomePath="$GENOME_PATH" 'BEGIN {FS=","}
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
    printf("%s,%s,%s,%s,%s,%s,",$1,$2,$3,$4,$5,$6)
    if ($5=="+") {
        printf("%s,",FW)
        revComp(RV)
        printf("\n")
    } else {
        revComp(FW)
        printf(",%s\n",RV)
    }
}' | gzip -c > "$PROJECT".both.reference.csv.gz;
gunzip -c "$PROJECT".both.reference.csv.gz | sort -t "," -k1,1 | sed -e 's/,/|/g' |  sqlite3 "$PROJECT".db '.import /dev/stdin mappingBoth'
gzip -f "$PROJECT".both.map.csv;

echo `date`" - Adjusting number of mismatches for C->T errors in mapping"
sqlite3 -csv "$PROJECT".db "SELECT mappingBoth.*, reads.sequenceFW, reads.sequenceRV
  FROM mappingBoth JOIN reads ON mappingBoth.id=reads.id;" | awk -v bp="$READ_LENGTH" 'BEGIN {FS=","} {
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
    print $1","$2","$3","$4","$5","mm","$7","$8
}' | gzip -c > "$PROJECT".both.adjust.csv.gz;
gunzip -c "$PROJECT".both.adjust.csv.gz | sort -t "," -k1,1 | sed -e 's/,/|/g' |  sqlite3 "$PROJECT".db '.import /dev/stdin mappingAdjust'

#Combine forward and reverse strand mappings, filter out duplicates
echo `date`" - Combining and filtering forward and reverse mappings";
gunzip -c "$PROJECT".both.adjust.csv.gz | sort -t "," -k 1,1 -k 6,6n  | awk -v maxmm=$MAX_MM -v mindiff=$MIN_MM_DIFF 'BEGIN {FS = ","} {
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
}' | gzip -c > "$PROJECT".map.csv.gz;

#import into the db
gunzip -c "$PROJECT".map.csv.gz | sort -t "," -k1,1 | sed -e 's/,/|/g' |  sqlite3 "$PROJECT".db '.import /dev/stdin mapping'

#export db into SAM files, keep strands separate
echo `date`" - Exporting database to BAM files"
cp "$GENOME_PATH"/hg18.samdir "$PROJECT".mappings.plus.sam;
sqlite3 -csv "$PROJECT".db "SELECT
mapping.id, mapping.chr, mapping.strand, mapping.positionFW, reads.sequenceFW, mapping.positionRV, reads.sequenceRV, reads.qualityFW, reads.qualityRV
FROM mapping LEFT JOIN reads ON mapping.id=reads.id WHERE mapping.strand='+';" | awk -v readLength="$READ_LENGTH" 'BEGIN {FS = ","} 
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
samtools import "$GENOME_PATH"/hg18.reflist "$PROJECT".mappings.plus.sam "$PROJECT".mappings.plus.bam;
samtools sort "$PROJECT".mappings.plus.bam "$PROJECT".mappings.plus.sort;
samtools index "$PROJECT".mappings.plus.sort.bam;
rm "$PROJECT".mappings.plus.bam "$PROJECT".mappings.plus.sam;

cp "$GENOME_PATH"/hg18.samdir "$PROJECT".mappings.minus.sam;
sqlite3 -csv "$PROJECT".db "SELECT
mapping.id, mapping.chr, mapping.strand, mapping.positionFW, reads.sequenceFW, mapping.positionRV, reads.sequenceRV, reads.qualityFW, reads.qualityRV
FROM mapping LEFT JOIN reads ON mapping.id=reads.id WHERE mapping.strand='-';" | awk -v readLength="$READ_LENGTH" 'BEGIN {FS = ","} 
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
samtools import "$GENOME_PATH"/hg18.reflist "$PROJECT".mappings.minus.sam "$PROJECT".mappings.minus.bam;
samtools sort "$PROJECT".mappings.minus.bam "$PROJECT".mappings.minus.sort;
samtools index "$PROJECT".mappings.minus.sort.bam;
rm "$PROJECT".mappings.minus.bam "$PROJECT".mappings.minus.sam;

#create bed & Rdata (GD) file for coverage mapping for each strand
echo `date`" - Creating coverage bed and GenomeData files"
sqlite3 -csv "$PROJECT".db "SELECT chr, positionFW, positionRV FROM mapping WHERE strand='+';" | gzip -c > "$PROJECT".plus.bed.gz
sqlite3 -csv "$PROJECT".db "SELECT chr, positionFW, positionRV FROM mapping WHERE strand='-';" | gzip -c > "$PROJECT".minus.bed.gz

R --vanilla --slave --args "$PROJECT".plus.bed.gz "$PROJECT".plus.Rdata "$READ_LENGTH" < "$PIPELINE_PATH"/scripts/bed2GD.R
R --vanilla --slave --args "$PROJECT".minus.bed.gz "$PROJECT".minus.Rdata "$READ_LENGTH" < "$PIPELINE_PATH"/scripts/bed2GD.R

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
}' | sort | uniq -c > "$PROJECT".context
