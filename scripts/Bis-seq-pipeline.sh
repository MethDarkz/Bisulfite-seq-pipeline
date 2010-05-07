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


