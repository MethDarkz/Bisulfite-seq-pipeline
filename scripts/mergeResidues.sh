#!/bin/bash -e
awk -v f1=s_1/s_1.residues.plus.gz -v f2=s_3/s_3.residues.plus.gz 'BEGIN {
	r1="/home/aarsta/analysis-utils/ChIP-seq-pipeline/scripts/residueChr.sh " f1
	r2="/home/aarsta/analysis-utils/ChIP-seq-pipeline/scripts/residueChr.sh " f2
	rs1=(r1 | getline var1)
	split(var1, spl1, ",")
	rs2=(r2 | getline var2)
	split(var2, spl2, ",")
	lastchr1=spl1[1]
	lastchr2=spl2[1]
	while((rs1==1)&&(rs2==1)) {
		if (spl1[1]!=spl2[1]) {
			if (lastchr1==spl1[1]) while (spl1[1]!=spl2[1]) {
				print var1
				lastchr1=spl1[1]
				rs1=(r1 | getline var1)
				split(var1, spl1, ",")				
			} else while (spl1[1]!=spl2[1]) {
				print var2
				lastchr2=spl2[1]
				rs2=(r2 | getline var2)
				split(var2, spl2, ",")
			}
		}
		if (spl1[2]==spl2[2]) {
			print spl1[1] "," spl1[2] "," (spl1[3]+spl2[3]) "," (spl1[4]+spl2[4]) "," (spl1[5]+spl2[5]) "," (spl1[6]+spl2[6]) "," (spl1[7]+spl2[7])
			lastchr1=spl1[1]
			rs1=(r1 | getline var1)
			split(var1, spl1, ",")
			lastchr2=spl2[1]
			rs2=(r2 | getline var2)
			split(var2, spl2, ",")
		} else if (spl1[2]>spl2[2]) {
			print var2
			lastchr2=spl2[1]
			rs2=(r2 | getline var2)
			split(var2, spl2, ",")

		} else {
			print var1
			lastchr1=spl1[1]
			rs1=(r1 | getline var1)
			split(var1, spl1, ",")
		}
	}
	while (rs1==1) {
		print var1
		rs1=(r1 | getline var1)
	}
	while (rs2==1) {
		print var2
		rs2=(r2 | getline var2)
	}
}' > testjoin.csv

