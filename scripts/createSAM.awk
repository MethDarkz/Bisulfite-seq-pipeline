BEGIN {
    FS = "|"
    #setup ASCII to numeric conversion table
    for (i=0; i<=255; i++) {
        t = sprintf("%c", i)
        ord[t] = i
    }
    if (strand=="+") {
        FWflag = 99
        RVflag = 147
        contigFlag = 0
    } else {
        FWflag = 83
        RVflag = 163
        contigFlag = 16   
    }
} 

function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}

function revComp(temp) {
    tempRev = ""
	for(i=length(temp);i>=1;i--) {
		tempChar = substr(temp,i,1)
		if (tempChar=="A") {tempRev = tempRev "T"} else
		if (tempChar=="C") {tempRev = tempRev "G"} else
		if (tempChar=="G") {tempRev = tempRev "C"} else
		if (tempChar=="T") {tempRev = tempRev "A"} else
		{tempRev = tempRev "N"}
	}
	return tempRev
}

function rev(s) {
  p = ""
  for(i=length(s); i > 0; i--) { p = p substr(s, i, 1) }
  return p
}

function printContig(FW, RV, FWq, RVq, offset) {
    newSeq = substr(FW, 1, offset)
    newQ = substr(FWq, 1, offset)
	for(i=1;i<(readLength-offset);i++) {
        if (ord[substr(FWq, offset+i, 1)]>ord[substr(RVq, i, 1)]) {fw=1} else #FW quality score higher
        if (ord[substr(FWq, offset+i, 1)]<ord[substr(RVq, i, 1)]) {fw=0} else #RV quality score higher
        fw = int(2*rand()) #choose randomly
        if (fw==1) {
            newSeq = newSeq substr(FW, offset+i, 1)
            newQ = newQ substr(FWq, offset+i, 1)
        } else {
            newSeq = newSeq substr(RV, i, 1)
            newQ = newQ substr(RVq, i, 1)
        }
    }
    newSeq = newSeq substr(RV, readLength-offset)
    newQ = newQ substr(RVq, readLength-offset)
    printf("%s\t%s\n", newSeq, newQ)
}

#Program entry point
{ 
    if (abs($6-$4)>=readLength) { #non-overlapping reads
        if (strand=="+") {
	        print $1 "\t" FWflag "\t" $2 "\t" $4 "\t255\t"readLength"M\t=\t" $6 "\t" (abs($6-$4)+readLength) "\t" $5 "\t" $8
	        printf("%s\t%s\t%s\t%s\t255\t"readLength"M\t=\t%s\t%s\t",$1,RVflag,$2,$6,$4,(abs($6-$4)+readLength))
	        printf("%s\t", revComp($7))
	        print(rev($9))
        } else {
	        print $1 "\t" FWflag "\t" $2 "\t" $4 "\t255\t"readLength"M\t=\t" $6 "\t" (abs($6-$4)+readLength) "\t" revComp($5) "\t" rev($8)
	        printf("%s\t%s\t%s\t%s\t255\t"readLength"M\t=\t%s\t%s\t",$1,RVflag,$2,$6,$4,(abs($6-$4)+readLength))
	        printf("%s\t", $7)
	        print($9)
        }
    } else { #ok merge reads into a contig
        if (strand=="+") {contigStart=$4} else {contigStart=$6}
        printf("%s\t%s\t%s\t%s\t255\t%sM\t*\t0\t0\t", $1, contigFlag, $2, contigStart, (abs($6-$4)+readLength))
        if (strand=="+") {
            printContig($5, revComp($7), $8, rev($9), abs($6-$4))
        } else {
            printContig($7, revComp($5), $9, rev($8), abs($6-$4))
        }
    }
}
