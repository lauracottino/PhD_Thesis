#!/bin/sh
#SBATCH --mem=50gb

ref_seq="Homo_sapiens_assembly38.fasta"
sig1="CFF0X.svsig.gz"
sig2="CWF0F.svsig.gz"
sig3="CNV0L.svsig.gz"
sig4="AB2543.svsig.gz"
sig5="DSX0S.svsig.gz"
sig6="AJK0K.svsig.gz"
sig7="AB0330.svsig.gz"
sig8="DEC0V.svsig.gz"
vcf="joint_8.vcf"

module load smrtlink
hostname
/usr/bin/time -f "%e %M" pbsv call -j 4 --ccs ${ref_seq} ${sig1} ${sig2} ${sig3} ${sig4} ${sig5} ${sig6} ${sig7} ${sig8} ${vcf}
bgzip ${vcf} 
tabix ${vcf}.gz 

