#!/bin/bash -e

if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  exit 1;
fi

source "$1";

if [ ! -e "$FASTQ" ]; then
  echo "$FASTQ"" does not exist!";
  exit 1;
fi

#init the database
echo "* Initialising the database";
rm -f "$PROJECT".db;
sqlite3 "$PROJECT".db < $PIPELINE_PATH/sql/ChIP-seq.schema;

#Convert fastq into csv, import into sqlite db
echo "* Importing reads into the database";
gunzip -c "$FASTQ"  | awk 'NR%4==1 || NR%4==2' | sed -e 's/ .*//' | sed -e 'N' -e 's/\n/,/' -e 's/^@//' >"$PROJECT".fastq.csv;
echo -e ".separator \",\"\n.import $PROJECT.fastq.csv reads" | sqlite3 "$PROJECT".db

gzip -f "$PROJECT".fastq.csv

#Map against the genome and filter out reads with too many MMs
echo "* Bowtie mapping";
gunzip -c "$PROJECT".fastq.gz | bowtie "$BOWTIE_PARAMS" "$BOWTIE_INDEX" - 2> mapping.log | gzip -c > "$PROJECT".map.gz;

gunzip -c "$PROJECT".map.gz | awk '{ \
num=split($8, tmp, ":"); \
if (num<'"$(($MAX_MM+$MIN_MM_DIFF))"') {print $1 "," $3 "," ($4+1) "," $2 "," num}}' > "$PROJECT".both.map.csv

#Filter out duplicates
echo "* Filtering mappings";
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


#Create GenomeData R instances of the mapped reads
echo "* Converting to GD Rdata";
R --vanilla --slave --quiet --args "$PROJECT".map.csv.gz "$PROJECT".Rdata < $PIPELINE_PATH/scripts/csvmap2GD.R &> /dev/null;

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


