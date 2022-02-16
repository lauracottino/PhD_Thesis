#!/bin/sh
#SBATCH --mem=50gb

ref_seq="/dataB/aux/38/Homo_sapiens_assembly38.fasta"
sig1="/spaces/lcottino/Pacbio/pbsv/vcfs/CFF0X_combined/CFF0X.svsig.gz"
sig2="/spaces/lcottino/Pacbio/pbsv/vcfs/CWF0F/CWF0F.svsig.gz"
sig3="/spaces/lcottino/Pacbio/pbsv/vcfs/CNV0L/CNV0L.svsig.gz"
sig4="/spaces/lcottino/Pacbio/pbsv/vcfs/AB2543/AB2543.svsig.gz"
sig5="/spaces/lcottino/Pacbio/pbsv/vcfs/DSX0S/DSX0S.svsig.gz"
sig6="/spaces/lcottino/Pacbio/pbsv/vcfs/AJK0K/AJK0K.svsig.gz"
sig7="/spaces/lcottino/Pacbio/pbsv/vcfs/AB0330/AB0330.svsig.gz"
sig8="/spaces/lcottino/Pacbio/pbsv/vcfs/DEC0V/DEC0V.svsig.gz"
vcf="joint_8_new.vcf"

module load smrtlink
hostname
/usr/bin/time -f "%e %M" pbsv call -j 4 --ccs ${ref_seq} ${sig1} ${sig2} ${sig3} ${sig4} ${sig5} ${sig6} ${sig7} ${sig8} ${vcf}
bgzip ${vcf} 
tabix ${vcf}.gz 

