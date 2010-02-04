#!/bin/bash
#$1 = filename of reads without fastq.gz
#$2 = tag size

if [ ! -e "$1".fastq.gz ]; then
  echo "$1"".fastq.gz does not exist!";
  exit 1;
fi

#init the database
rm -f "$1".db;
sqlite3 "$1".db < /home/aarsta/analysis-utils/ChIP-seq-pipeline/sql/ChIP-seq.schema;

#Convert fastq into csv, import into sqlite db
gunzip -c "$1".fastq.gz  | awk 'NR%4==1 || NR%4==2' | sed -e 'N' -e 's/\n/,/' -e 's/^@//'>"$1".fastq.csv;
echo -e ".separator \",\"\n.import $1.fastq.csv reads" | sqlite3 "$1".db

#Record position of C residues in reads
awk 'BEGIN {FS = ","}; {\
num=split($2, tmp, "C"); \
upto=0
for (i = 2; i <= num; i++) {\
   upto = upto + 1 + length(tmp[i-1]);
   print $1 "," upto;  \
}
}' "$1".fastq.csv > "$1".C.csv;
echo -e ".separator \",\"\n.import $1.C.csv readsC" | sqlite3 "$1".db
gzip "$1".fastq.csv
gzip "$1".C.csv

#Convert C residues in reads to T
gunzip -c "$1".fastq.gz | sed -e '2~4s/C/T/g' | gzip -c > "$1".conv.fastq.gz;

#map against plus bisulfite genome (no rc)
gunzip -c "$1".conv.fastq.gz | bowtie --norc --best -y -p 7 -v 3 --un unaligned.plus.fastq -m 1 --max multiple.plus.fastq --strata /home/aarsta/data/genome/hg18/bowtie/bis.hg18.plus - 2> mapping.plus.log | gzip -c > "$1".plus.map.gz;
bowtie --norc --best -y -p 7 -v 3 -k 1 /home/aarsta/data/genome/hg18/bowtie/bis.hg18.plus multiple.plus.fastq | gzip -c > "$1"_multiple.plus.map.gz;
gzip unaligned.plus.fastq;
gzip multiple.plus.fastq;

#map against minus bisulfite genome (no fw)
gunzip -c "$1".conv.fastq.gz | bowtie --nofw --best -y -p 7 -v 3 --un unaligned.minus.fastq -m 1 --max multiple.minus.fastq --strata /home/aarsta/data/genome/hg18/bowtie/bis.hg18.minus - 2> mapping.minus.log | gzip -c > "$1".minus.map.gz;
bowtie --nofw --best -y -p 7 -v 3 -k 1 /home/aarsta/data/genome/hg18/bowtie/bis.hg18.minus multiple.minus.fastq | gzip -c > "$1"_multiple.minus.map.gz;
gzip unaligned.minus.fastq;
gzip multiple.minus.fastq;

#convert maps into csv and import into db
gunzip -c "$1".plus.map.gz | awk '{print $1 "," $3 "," $4 "," $2}' > "$1".plus.map.csv
echo -e ".separator \",\"\n.import ""$1"".plus.map.csv mappingFW" | sqlite3 "$1".db
gzip "$1".plus.map.csv

gunzip -c "$1".minus.map.gz | awk '{print $1 "," $3 "," $4 "," $2}' > "$1".minus.map.csv
echo -e ".separator \",\"\n.import ""$1"".minus.map.csv mappingRV" | sqlite3 "$1".db
gzip "$1".minus.map.csv

gunzip -c "$1"_multiple.plus.map.gz | awk '{print $1 "," $3 "," $4 "," $2}' > "$1"_multiple.plus.map.csv
echo -e ".separator \",\"\n.import ""$1""_multiple.plus.map.csv mappingFWMulti" | sqlite3 "$1".db
gzip "$1"_multiple.plus.map.csv

gunzip -c "$1"_multiple.minus.map.gz | awk '{print $1 "," $3 "," $4 "," $2}' > "$1"_multiple.minus.map.csv
echo -e ".separator \",\"\n.import ""$1""_multiple.minus.map.csv mappingRVMulti" | sqlite3 "$1".db
gzip "$1"_multiple.minus.map.csv

#Suck out unique FW strand mappings
sqlite3 -csv "$1".db "SELECT * \
FROM mappingFW \
WHERE id NOT IN(SELECT id FROM mappingRV) \
AND id NOT IN(SELECT id FROM mappingFWMulti) \
AND id NOT IN(SELECT id FROM mappingRVMulti);" | sed 's/plusU_//' > "$1".map.csv

#Suck out unique RV strand mappings
sqlite3 -csv "$1".db "SELECT * \
FROM mappingRV \
WHERE id NOT IN(SELECT id FROM mappingFW) \
AND id NOT IN(SELECT id FROM mappingFWMulti) \
AND id NOT IN(SELECT id FROM mappingRVMulti);" | sed 's/minusU_//' >> "$1".map.csv

#import into the db
echo -e ".separator \",\"\n.import ""$1"".map.csv mapping" | sqlite3 "$1".db
gzip "$1".map.csv

#create mappingC table (minus strand C residues are shifted 1bp left
sqlite3 -csv "$1".db "SELECT mapping.id, mapping.chr, mapping.position, mapping.strand, readsC.position \
FROM readsC JOIN mapping ON readsC.id=mapping.id" | awk 'BEGIN {FS = ","}; {\
	printf "%s,%s,",$1,$2;
	if ($4=="+") {printf "%s,",($3+$5);} else {printf "%s,",($3+'"$1"'-$5);}
	print $4
	}' > "$1".mapC.csv
echo -e ".separator \",\"\n.import ""$1"".mapC.csv mappingC" | sqlite3 "$1".db
gzip "$1".mapC.csv


