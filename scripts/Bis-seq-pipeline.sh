#!/bin/bash -e 

#set default paths - if these tools are within the $PATH then this script will work fine without them being defined in the .config
BOWTIE_PATH="bowtie"
SAMTOOLS_PATH="samtools"
R_PATH="R"

#Load the config file
if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  echo "Usage: Bis-seq-pipeline.sh [configfile]";
  exit 1;
fi

echo `date`" - Reading config file: $1";
source "$1";

#Check for existence of necessary files
if [ ! -e "$FASTQ" ]; then echo "$FASTQ does not exist!"; exit 1; fi
if [ ! -e "$PIPELINE_PATH"/sql/Bis-seq.schema ]; then echo "$PIPELINE_PATH""/sql/Bis-seq-PE.schema does not exist!"; exit 1; fi
if [ ! -e "$PIPELINE_PATH"/scripts/readChr.awk ]; then echo "$PIPELINE_PATH""/scripts/readChr.awk does not exist!"; exit 1; fi
if [ ! -e "$GENOME_PATH"/plus.1.ebwt ]; then echo "Forward bowtie index not found at $GENOME_PATH"; exit 1; fi
if [ ! -e "$GENOME_PATH"/minus.1.ebwt ]; then echo "Reverse bowtie index not found at $GENOME_PATH"; exit 1; fi

#check for existence of necessary tools 
if [ ! type -P sqlite3 &>/dev/null ]; then echo "sqlite3 command not found."; exit 1; fi
if [ ! type -P "$BOWTIE_PATH" &>/dev/null ]; then echo "$BOWTIE_PATH command not found."; exit 1; fi
if [ ! type -P "$SAMTOOLS_PATH" &>/dev/null ]; then echo "$SAMTOOLS_PATH command not found."; exit 1; fi
if [ ! type -P "$R_PATH" &>/dev/null ]; then echo "$R_PATH command not found."; exit 1; fi

#init the database
echo `date`" - Initialising the database";
rm -f "$PROJECT".db;
sqlite3 "$PROJECT".db < "$PIPELINE_PATH"/sql/Bis-seq.schema;

echo `date`" - Importing reads into the database";
#Import seqs into the db, one row per read id
gunzip -c "$FASTQ" |awk 'BEGIN {
  FS="/"
} {
  readname=$1
  getline seq
  getline
  getline qual
  print readname "|" seq "|" qual
}' | sed 's/^@//' | sort -t "|" -k1,1 | sqlite3 "$PROJECT".db '.import /dev/stdin reads';

#Convert C residues in reads to T
echo `date`" - Bisulfite converting reads";
gunzip -c "$FASTQ" | sed -e 's/ .*//' -e '2~4s/C/T/g' | gzip -c > "$PROJECT".conv.fastq.gz;

#Map against forward strand, and filter out reads with too many MMs
echo `date`" - Bowtie mapping against forward strand";
gunzip -c "$PROJECT".conv.fastq | "$BOWTIE_PATH" --norc "$BOWTIE_PARAMS" "$GENOME_PATH"/plus - 2> mapping.plus.log | awk -v maxmm=$(($MAX_MM+$MIN_MM_DIFF)) 'BEGIN {FS="\t"}{
  num=split($8, tmp, ":")-1
  if (num==-1) num=0
  if (num<maxmm) {
    print substr($1,1,length($1)-2) "|" $3 "|" ($4+1) "|+|" num
  }
}' | sed 's/plusU_//'> "$PROJECT".both.map;

#Same for reverse strand
echo `date`" - Bowtie mapping against reverse strand";
gunzip -c "$PROJECT".conv.fastq | "$BOWTIE_PATH" --nofw "$BOWTIE_PARAMS" "$GENOME_PATH"/minus - 2> mapping.plus.log | awk -v maxmm=$(($MAX_MM+$MIN_MM_DIFF)) 'BEGIN {FS="\t"}{
  num=split($8, tmp, ":")-1
  if (num==-1) num=0
  if (num<maxmm) {
    print substr($1,1,length($1)-2) "|" $3 "|" ($4+1) "|-|" num
  }
}' | sed 's/minusU_//'> "$PROJECT".both.map;

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
    seq=toupper(substr(chrSeq,$3,readLength)) #retrieve forward sequence
    printf("%s|%s|%s|%s|%s|",$1,$2,$3,$4,$5)
    if ($4=="+") {
        printf("%s\n",seq)
    } else {
        revComp(seq)
        printf("\n")
    }
}' | sort -t "|" -k1,1 |  sqlite3 "$PROJECT".db '.import /dev/stdin mappingBoth'
gzip -f "$PROJECT".both.map;

echo `date`" - Adjusting number of mismatches for C->T errors in mapping"
sqlite3 "$PROJECT".db "SELECT mappingBoth.*, reads.sequence
  FROM mappingBoth JOIN reads ON mappingBoth.id=reads.id;" | awk -v bp="$READ_LENGTH" 'BEGIN {FS="|"} {
    mm=0;
    for(i=1;i<=bp;i++) {
        temp1 = substr($6,i,1)
        temp2 = substr($7,i,1)
        if (temp1=="C") {
            if (temp2!="C"&&temp2!="T") mm++;
        } else if (temp1!=temp2) mm++;
    }
    print $1"|"$2"|"$3"|"$4"|"mm"|"$6
}' | gzip -c > "$PROJECT".both.adjust.csv.gz;
gunzip -c "$PROJECT".both.adjust.csv.gz | sort -t "|" -k1,1 | sqlite3 "$PROJECT".db '.import /dev/stdin mappingAdjust'

echo `date`" - Combining  forward and reverse mappings and filtering out duplicate mappings";
gunzip -c "$PROJECT".both.adjust.csv.gz | sort -t "|" -k 1,1 -k 5,5n  | awk -v maxmm=$MAX_MM -v mindiff=$MIN_MM_DIFF 'BEGIN {FS = "|"} {
    s = $1
    if (s != prevs) {
        if ( FNR > 1 ) {
            if ( prevval<=maxmm && valdiff>=mindiff) print prevline
        }
        prevval = $5
        prevline = $0
        valdiff = mindiff
    }
    else valdiff = prevval-$5
    prevs = s
}
END {
    if ( prevval<=maxmm && valdiff>=mindiff) print prevline
}' | sort -t "|" -k1,1 |  sqlite3 "$PROJECT".db '.import /dev/stdin mapping'

echo `date`" - Exporting database to BAM files"
cp "$GENOME_PATH"/samdir "$PROJECT".mappings.plus.sam;
sqlite3 "$PROJECT".db "SELECT
mapping.id, mapping.chr, mapping.strand, mapping.position, reads.sequence, reads.quality
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
"$SAMTOOLS_PATH" import "$GENOME_PATH"/reflist "$PROJECT".mappings.plus.sam "$PROJECT".mappings.plus.bam;
"$SAMTOOLS_PATH" sort "$PROJECT".mappings.plus.bam "$PROJECT".plus;
"$SAMTOOLS_PATH" index "$PROJECT".plus.bam;
rm "$PROJECT".mappings.plus.bam "$PROJECT".mappings.plus.sam;

cp "$GENOME_PATH"/samdir "$PROJECT".mappings.minus.sam;
sqlite3 "$PROJECT".db "SELECT
mapping.id, mapping.chr, mapping.strand, mapping.position, reads.sequence, reads.quality
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
"$SAMTOOLS_PATH" import "$GENOME_PATH"/reflist "$PROJECT".mappings.minus.sam "$PROJECT".mappings.minus.bam;
"$SAMTOOLS_PATH" sort "$PROJECT".mappings.minus.bam "$PROJECT".minus;
"$SAMTOOLS_PATH" index "$PROJECT".minus.bam;
rm "$PROJECT".mappings.minus.bam "$PROJECT".mappings.minus.sam;

#create bed & Rdata (GD) file for coverage mapping for each strand
echo `date`" - Creating coverage bed and GenomeData files";
sqlite3 -csv "$PROJECT".db "SELECT chr, positionFW, positionRV FROM mapping WHERE strand='+';" | gzip -c > "$PROJECT".plus.bed.gz;
sqlite3 -csv "$PROJECT".db "SELECT chr, positionFW, positionRV FROM mapping WHERE strand='-';" | gzip -c > "$PROJECT".minus.bed.gz;

"$R_PATH" --vanilla --slave --args "$PROJECT".plus.bed.gz "$PROJECT".plus.Rdata "$READ_LENGTH" < "$PIPELINE_PATH"/scripts/bed2GD.R;
"$R_PATH" --vanilla --slave --args "$PROJECT".minus.bed.gz "$PROJECT".minus.Rdata "$READ_LENGTH" < "$PIPELINE_PATH"/scripts/bed2GD.R;

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
    print "No of C residues in the read sequences: " readC
    print "No of C residues in the reference sequences: " refC
}' >> "$PROJECT".context;

NUM_READS=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM reads;"`;
NUM_REPORTED=`sqlite3 "$PROJECT".db "SELECT COUNT(DISTINCT id) FROM mappingBoth;"`;
NUM_UNIQUE=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM mapping;"`;
NUM_UNMAPPABLE=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM reads \
WHERE id NOT IN(SELECT id FROM mappingBoth);"`;
NUM_MULTIPLE=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM reads \
WHERE id IN(SELECT id FROM mappingBoth) \
AND id NOT IN (SELECT id FROM mapping);"`;

#Creating mapping log
echo `date`" - Creating mapping log"
echo "# reads processed: $NUM_READS" > mapping.log;
echo "# reads with at least one reported alignment: $NUM_REPORTED" >> mapping.log;
echo "# reads that failed to align: $NUM_UNMAPPABLE" >> mapping.log;
echo "# reads with alignments suppressed due to -m: $NUM_MULTIPLE" >> mapping.log;
echo "Reported $NUM_UNIQUE alignments to 1 output stream(s)" >> mapping.log;
#!/bin/bash -e

if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  exit 1;
fi

source "$1";

if [! -e "$FASTQ"]; then
  echo "$FASTQ"" does not exist!";
  exit 1;
fi


#init the database
echo "* Initialising the database";
rm -f "$PROJECT".db;
sqlite3 "$PROJECT".db < $PIPELINE_PATH/sql/Bis-seq.schema;

#Convert fastq into csv, import into sqlite db
echo "* Importing reads into the database";
gunzip -c "$FASTQ"  | awk 'NR%4==1 || NR%4==2' | sed -e 's/ .*//' | sed -e 'N' -e 's/\n/,/' -e 's/^@//' >"$PROJECT".fastq.csv;
echo -e ".separator \",\"\n.import $PROJECT.fastq.csv reads" | sqlite3 "$PROJECT".db

#Record position of C residues in reads
echo "* Importing C residues into the database";
awk 'BEGIN {FS = ","}; {\
num=split($2, tmp, "C"); \
upto=0
for (i = 2; i <= num; i++) {\
   upto = upto + 1 + length(tmp[i-1]);
   print $1 "," upto;  \
}
}' "$PROJECT".fastq.csv > "$PROJECT".C.csv;
echo -e ".separator \",\"\n.import $PROJECT.C.csv readsC" | sqlite3 "$PROJECT".db
gzip -f "$PROJECT".fastq.csv
gzip -f "$PROJECT".C.csv

#Convert C residues in reads to T
echo "* Bisulfite converting reads";
gunzip -c "$FASTQ" | sed -e 's/ .*//' -e '2~4s/C/T/g' > "$PROJECT".conv.fastq;

#Map against forward strand, and filter out reads with too many MMs
echo "* Bowtie mapping against forward strand";
bowtie --norc "$BOWTIE_PARAMS" "$BOWTIE_BS_INDEX".plus "$PROJECT".conv.fastq 2> mapping.plus.log | gzip -c > "$PROJECT".plus.map.gz;

gunzip -c "$PROJECT".plus.map.gz | awk '{ \
num=split($8, tmp, ":"); \
if (num<'"$(($MAX_MM+$MIN_MM_DIFF))"') {print $1 "," $3 "," ($4+1) "," $2 "," num}}' | sed 's/plusU_//'> "$PROJECT".both.map.csv

#Same for reverse strand
echo "* Bowtie mapping against reverse strand";
bowtie --nofw "$BOWTIE_PARAMS" "$BOWTIE_BS_INDEX".minus "$PROJECT".conv.fastq 2> mapping.minus.log | gzip -c > "$PROJECT".minus.map.gz;
gzip "$PROJECT".conv.fastq;

gunzip -c "$PROJECT".minus.map.gz | awk '{ \
num=split($8, tmp, ":"); \
if (num<'"$(($MAX_MM+$MIN_MM_DIFF))"') {print $1 "," $3 "," ($4+50) "," $2 "," num}}' | sed 's/minusU_//'>> "$PROJECT".both.map.csv

#Combine forward and reverse strand mappings, filter out duplicates
echo "* Combining and filtering forward and reverse mappings";
sort -t "," -k 1,1 -k 5,5n "$PROJECT".both.map.csv  | awk 'BEGIN {FS = ","} {
	line1 = $0;
	split(line1, line1s, ",");
	while ((getline line2) > 0) {
		split(line1, line1s, ",");
		split(line2, line2s, ",");
		if (line2s[1]==line1s[1]) {
			if ((line2s[5]-line1s[5])>='"$MIN_MM_DIFF"') {
				if (line1s[5]<='$MAX_MM') print line1;
			}
			while(line2s[1]==line1s[1]) {
				getline line1;
				split(line1, line1s, ",");
			}
		} else {
			if (line1s[5]<='$MAX_MM') print line1;
			line1 = line2;
		}
	}
	print line1;
}' | cut -d "," --complement -f5  > "$PROJECT".map.csv

#import into the db
echo -e ".separator \",\"\n.import ""$PROJECT"".both.map.csv mappingBoth" | sqlite3 "$PROJECT".db
gzip -f "$PROJECT".both.map.csv

echo -e ".separator \",\"\n.import ""$PROJECT"".map.csv mapping" | sqlite3 "$PROJECT".db
gzip -f "$PROJECT".map.csv

#create mappingC table
echo "* Creating mappingC table";
sqlite3 -csv "$PROJECT".db "SELECT mapping.id, mapping.chr, mapping.position, mapping.strand, readsC.position \
FROM readsC JOIN mapping ON readsC.id=mapping.id" | awk 'BEGIN {FS = ","}; {\
	printf "%s,%s,",$1,$2;
	if ($4=="+") {printf "%s,",($3+$5-1);} else {printf "%s,",($3-$5);}
	print $4
	}' > "$PROJECT".mapC.csv
echo -e ".separator \",\"\n.import ""$PROJECT"".mapC.csv mappingC" | sqlite3 "$PROJECT".db
gzip -f "$PROJECT".mapC.csv

#Create GenomeData R instances of the mapped reads and their Cs
echo "* Converting to GD Rdata";
R --vanilla --slave --quiet --args "$PROJECT".map.csv.gz "$PROJECT".Rdata < $PIPELINE_PATH/scripts/csvmap2GD.R &> /dev/null;
R --vanilla --slave --quiet --args "$PROJECT".mapC.csv.gz "$PROJECT".C.Rdata < $PIPELINE_PATH/scripts/csvmap2GD.R &> /dev/null;

#create mapping.log
echo "* Creating mapping.log";
NUM_READS=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM reads;"`;
NUM_UNIQUE=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM mapping;"`;
NUM_UNMAPPABLE=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM reads \
WHERE id NOT IN(SELECT id FROM mappingBoth);"`;
NUM_MULTIPLE=`sqlite3 "$PROJECT".db "SELECT COUNT(id) FROM reads \
WHERE id IN(SELECT id FROM mappingBoth) \
AND id NOT IN (SELECT id FROM mapping);"`;

echo "# reads processed: $NUM_READS" > mapping.log;
echo "# reads with at least one reported alignment: $NUM_UNIQUE" >> mapping.log;
echo "# reads that failed to align: $NUM_UNMAPPABLE" >> mapping.log;
echo "# reads with alignments suppressed due to -m: $NUM_MULTIPLE" >> mapping.log;
echo "Reported $NUM_UNIQUE alignments to 1 output stream(s)" >> mapping.log;


