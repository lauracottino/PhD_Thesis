#!/bin/bash
#SBATCH --mem=1gb

VCF=/spaces/emma/Build38Run3/PCRfreeRun/VCFs/
classpath="/opt/exp_soft/bioinf/svtoolkit/lib/gatk/GenomeAnalysisTK.jar"
java -cp ${classpath} org.broadinstitute.gatk.tools.CatVariants \
    -R /dataB/aux/38/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta \
    -V ${VCF}/chr1.vcf \
    -V ${VCF}/chr2.vcf \
    -V ${VCF}/chr3.vcf \
    -V ${VCF}/chr4.vcf \
    -V ${VCF}/chr5.vcf \
    -V ${VCF}/chr6.vcf \
    -V ${VCF}/chr7.vcf \
    -V ${VCF}/chr8.vcf \
    -V ${VCF}/chr9.vcf \
    -V ${VCF}/chr10.vcf \
    -V ${VCF}/chr11.vcf \
    -V ${VCF}/chr12.vcf \
    -V ${VCF}/chr13.vcf \
    -V ${VCF}/chr14.vcf \
    -V ${VCF}/chr15.vcf \
    -V ${VCF}/chr16.vcf \
    -V ${VCF}/chr17.vcf \
    -V ${VCF}/chr18.vcf \
    -V ${VCF}/chr19.vcf \
    -V ${VCF}/chr20.vcf \
    -V ${VCF}/chr21.vcf \
    -V ${VCF}/chr22.vcf \
    -out merged.vcf \
    -assumeSorted
