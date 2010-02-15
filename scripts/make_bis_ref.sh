#!/bin/bash
mkdir -p plusU plusM minusU minusM;
for file in *.fa;
do
	sed -e '1 s/chr/plusU_chr/' -e '2,$ s/[Cc]/t/g' "$file" > plusU/"$file";
	sed -e '1 s/chr/plusM_chr/' -e '2,$ s/[Cc][Gg]/Mg/g' -e '/[Cc]$/{
N
s/[Cc]\n[Gg]/M\
g/
}' -e '2,$ s/[Cc]/t/g' -e '2,$ s/M/c/g' "$file" > plusM/"$file";
	sed -e '1 s/chr/minusU_chr/' -e '2,$ s/[Gg]/a/g' "$file" > minusU/"$file"; 
	sed -e '1 s/chr/minusM_chr/' -e '2,$ s/[Cc][Gg]/cM/g' -e '/[Cc]$/{
N
s/[Cc]\n[Gg]/c\
M/
}' -e '2,$ s/[Gg]/a/g' -e '2,$ s/M/g/g' "$file" > minusM/"$file";
done


