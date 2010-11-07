#!/bin/bash -e
if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  exit 1;
fi

source "$1";

mkdir -p plusU minusU;
rm -f reflengths reflist CpGsites.csv Csites.csv Gsites.csv;
grep -h "^>" *.fa | sed 's/^>//' > reflist;
for file in *.fa;
do
	echo -e "Converting ""$file";
#On the plus strand of every chromosome, convert all C's to Ts'
	sed -e '/^>/s/>/>plusU_/' -e '/^>/!s/[Cc]/t/g' "$file" > plusU/"$file";

#On the minus strand of every chromosome, convert all G's to As'
	sed -e '/^>/s/>/>minusU_/' -e '/^>/!s/[Gg]/a/g' "$file" > minusU/"$file";
	
#Record the size of each chromosome in basepairs
	echo -e "${file%.fa}""\t$(( $(awk -f "$PIPELINE_PATH"/scripts/readChr.awk "$file" | wc -c) - 1 ))" >> reflengths;

#Record the position of every CpG site
    awk -f "$PIPELINE_PATH"/scripts/readChr.awk $file | grep -bo [cCgG] | awk -v pipeline="$PIPELINE_PATH" -v chr="${file%.fa}" 'BEGIN {
        "awk -f "pipeline"/scripts/readChr.awk "chr".fa" | getline chrSeq
        FS=":";
    } {
        temp = toupper(substr(chrSeq, $1, 3))
        if (substr(temp, 2, 1)=="C") {
            if (substr(temp, 3, 1)=="G") print chr "," $1+1 ",CG" >> "CpGsites.csv"
            else print chr "," $1+1 ",CH" >> "Csites.csv"
        } else if (substr(temp, 1, 1)!="C") print chr "," $1 ",DG" >> "Gsites.csv"
    }'
done

echo -e "Creating bowtie index of the plus strand";
ls plusU/*.fa | sed -e :a -e '$!N;s/\n/,/;ta' -e 's/$/ plus/' | xargs -L 2 "$BOWTIE_BUILD";

echo -e "Creating bowtie index of the minus strand";
ls minusU/*.fa | sed -e :a -e '$!N;s/\n/,/;ta' -e 's/$/ minus/' | xargs -L 2 "$BOWTIE_BUILD";
