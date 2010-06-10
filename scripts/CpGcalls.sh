#!/bin/bash -e

#set default paths - if these tools are within the $PATH then this script will work fine without them being defined in the .config
SAMTOOLS_PATH="samtools"

#Load the config file
if [ ! -e "$1" ]; then
  echo "$1"" does not exist!";
  echo "Usage: CpGcalls.sh [configfile]";
  exit 1;
fi

echo `date`" - Reading config file: $1";
source "$1";

if [ ! -e "$GENOME_PATH/CpGsites.csv" ]; then echo "$GENOME_PATH/CpGsites.csv not found"; exit 1; fi
if [ ! -e "$GENOME_PATH/reflist" ]; then echo "$GENOME_PATH/reflist not found"; exit 1; fi
if [ ! type -P "$SAMTOOLS_PATH" &>/dev/null ]; then echo "$SAMTOOLS_PATH command not found."; exit 1; fi

echo `date`" - Exporting plus strand to $PROJECT.plus.CpG.gz";
$SAMTOOLS_PATH pileup "$PROJECT".plus.bam | cut -f 1,2,5 | awk -v genomeDir="$GENOME_PATH" -v strand="+" '
function nextCpG() {
    if((getline tmp < (genomeDir"/CpGsites.csv"))==0) exit
    split(tmp, tmp2, ",")
    chr=tmp2[1]
    chrNo=chrTable[chr]
    pos=tmp2[2]
    if (strand=="-") pos++
} BEGIN {
    i=1
    while (getline temp < (genomeDir"/reflist")) chrTable[temp] = i++
    nextCpG()
} {
    while (chrTable[$1]>chrNo) nextCpG()
    while ($2>pos&&$1==chr) nextCpG()
    if ($2==pos) { #print chr, pos, A,C,G,T,N
        residues=toupper($3)
        A=C=G=T=N=0
        for (i=1;i<=length(residues);i++) {
            residue=substr(residues,i,1)
            if (residue=="A") A++;
            if (residue=="C") C++;
            if (residue=="G") G++;
            if (residue=="T") T++;
            if (residue=="N") N++;
        }
        print chr","pos","A","C","G","T","N
    }
}' | gzip -c > "$PROJECT".plus.CpG.gz;

echo `date`" - Exporting minus strand to $PROJECT.minus.CpG.gz";
$SAMTOOLS_PATH pileup "$PROJECT".minus.bam | cut -f 1,2,5 | awk -v genomeDir="$GENOME_PATH" -v strand="-" '
function nextCpG() {
    if((getline tmp < (genomeDir"/CpGsites.csv"))==0) exit
    split(tmp, tmp2, ",")
    chr=tmp2[1]
    chrNo=chrTable[chr]
    pos=tmp2[2]
    if (strand=="-") pos++
} BEGIN {
    i=1
    while (getline temp < (genomeDir"/reflist")) chrTable[temp] = i++
    nextCpG()
} {
    while (chrTable[$1]>chrNo) nextCpG()
    while ($2>pos&&$1==chr) nextCpG()
    if ($2==pos) { #print chr, pos, A,C,G,T,N
        residues=toupper($3)
        A=C=G=T=N=0
        for (i=1;i<=length(residues);i++) {
            residue=substr(residues,i,1)
            if (residue=="A") A++;
            if (residue=="C") C++;
            if (residue=="G") G++;
            if (residue=="T") T++;
            if (residue=="N") N++;
        }
        print chr","pos","A","C","G","T","N
    }
}' | gzip -c > "$PROJECT".minus.CpG.gz;
