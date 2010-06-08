#!/bin/bash -e
if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  exit 1;
fi

source "$1";

mkdir -p plusU minusU;
rm -f samdir reflist CpGsites.csv;
grep -h "^>" *.fa | sed 's/^>//' > reflist;
for file in *.fa;
do
	echo -e "Converting ""$file";
#On the plus strand of every chromosome, convert all C's to Ts'
	sed -e '/^>/s/>/>plusU_/' -e '/^>/!s/[Cc]/t/g' "$file" > plusU/"$file";

#On the minus strand of every chromosome, convert all G's to As'
	sed -e '/^>/s/>/>minusU_/' -e '/^>/!s/[Gg]/a/g' "$file" > minusU/"$file";
	
#Record the size of each chromosome in basepairs
	echo -e "@SQ\tSN:""${file%.fa}""\tLN:$(( $(awk -f "$PIPELINE_PATH"/scripts/readChr.awk "$file" | wc -c) - 1 ))" >> samdir;

#Record the position of every CpG site
	awk -f ~/analysis-utils/ChIP-seq-pipeline/scripts/readChr.awk "$file" | grep -bo [cC][gG] | awk -v chr="${file%.fa}" 'BEGIN {FS=":"} {print chr "," ($1+1)}' >> CpGsites.csv;
done

echo -e "Creating bowtie index of the plus strand";
ls plusU/*.fa | sed -e :a -e '$!N;s/\n/,/;ta' -e 's/$/ plus/' | xargs -L 2 "$BOWTIE_PATH"/bowtie-build;

echo -e "Creating bowtie index of the minus strand";
ls minusU/*.fa | sed -e :a -e '$!N;s/\n/,/;ta' -e 's/$/ minus/' | xargs -L 2 "$BOWTIE_PATH"/bowtie-build;
