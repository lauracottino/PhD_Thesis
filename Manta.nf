   
#!/usr/bin/env nextflow

ref = file("/dataB/aux/38/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta", type:'file')
refi = file("/dataB/aux/38/Homo_sapiens_assembly38/Homo_sapiens_assembly38.fasta.fai", type:'file')


mixed = "/dataJ/build38/datasets/*/bam/*/*"


mixed_pairs = Channel.fromFilePairs("${mixed}{.bam,.bai}", size:2)
mixed_bams  = Channel.fromPath("${mixed}{.bam,.bai}")

kg = "/external/diskC/build38/datasets/1000G/*"
kg_pairs = Channel.fromFilePairs("${kg}{.bam,.bai}", size:2) { file -> file.simpleName }
kg_bams  = Channel.fromPath("${kg}{.bam,.bai}")



pairs = kg_pairs.mix(mixed_pairs)
bams  = kg_bams.mix(mixed_bams)

params.variant = 'diploidSV'
params.out     = "out"
params.manta_parallel = 15
params.graphtyper_parallel = 12
params.first_chrom=1
out = params.out


process runManta {
	cpus params.manta_parallel
        memory '5GB'
	tag {idSample}
	maxForks 60
	input:
	  file(ref) 
	  file(refi)
	  set idSample, file(bam) from pairs
        publishDir "run_temp/manta_vcf", mode:'copy'
	output:
	  set idSample, file("${idSample}.vcf.gz"),
      	                file("${idSample}.vcf.gz.tbi") \
            into mantaOutput 
	"""
	configManta.py \
	    --bam ${idSample}.bam \
	    --referenceFasta ${ref} \
	    --runDir .
	./runWorkflow.py -j ${params.manta_parallel} -g 5
        cp results/variants/${params.variant}.vcf.gz ${idSample}.vcf.gz
        cp results/variants/${params.variant}.vcf.gz.tbi ${idSample}.vcf.gz.tbi
	"""
}

     

chroms = (params.first_chrom..22).collect { "chr$it" }

chrom_list = chroms.join(" ")


// See graphtyper bug https://github.com/DecodeGenetics/graphtyper/issues/55
//    for explanation of the sed line
process combineNewSiteVCFs {
  cpus  8 
  input:
     file(inputs) from mantaOutput.collect()
  output:
     set file("${out}.vcf.gz"), file("${out}.vcf.gz.tbi") into mergedVCF
  script:
       
    """
      ls *vcf.gz > input_vcfs
      svimmer --threads 8 input_vcfs $chrom_list --output  ${out}.vcf
      sed -i "/chr6.*BND/d" ${out}.vcf 
      bgzip ${out}.vcf
      tabix ${out}.vcf.gz
    """
}


process graphTyper {
  errorStrategy 'finish'
  cpus params.graphtyper_parallel
  memory '300GB'
  input:
    file(bams) from bams.collect()     
    file(ref)
    file(refi)
    set file(merged), file(tbi) from mergedVCF
    each chrom from chroms
  output:
    file("sv_results/${chrom}") into calledIndivVCFs 
  script:
     """
     #!/bin/bash
     ls *bam > bamlist
     graphtyper genotype_sv $ref $merged --sams=bamlist --region=$chrom \
                 --threads=${params.graphtyper_parallel}
     """

}

process  mergeCalledByChromosome {
  input:
    file(chrom) from calledIndivVCFs
  output:
    set file("${chrom}.vcf.gz"), file("${chrom}.vcf.gz.tbi") into chromResult
  script:
     """
      vcf-concat $chrom/*vcf.gz | bgzip -c > ${chrom}.vcf.gz
      tabix ${chrom}.vcf.gz
     """
}


chrom_vcf  = (params.first_chrom..22).collect { "chr${it}.vcf.gz" }.join(" ")

process mergeAllChroms {
  input:
    file(all) from chromResult.collect()
  output:
    set file("${out}.vcf.gz"), file("${out}.vcf.gz.tbi") into allResult
  publishDir "results", mode:'copy'
  script:
    """
      vcf-concat ${chrom_vcf} | bgzip -c > ${out}.vcf.gz
      tabix ${out}.vcf.gz
    """
}
