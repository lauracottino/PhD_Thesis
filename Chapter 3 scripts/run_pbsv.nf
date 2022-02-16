// (c) University of the Witwatersrand, Johannesburg
// Author Scott Hazelhurst
// Released under https://creativecommons.org/licenses/by-sa/4.0/




params.output = "sample"
params.has_bam = false
trf = file(params.tandem_example)

params.chrom_prefix="chr"

fnames = params.input.tokenize(",")


ref_seq = file(params.ref_seq)
ref_fai = file(params.ref_fai)
ref_mmi = file(params.ref_mmi)



if (params.has_bam) {
  got_bams = true
  fnames.each { 
     base=file(it).baseName; 
     if (! file("${params.bam}/${base}.bam").exists()) {
       println "No file ${params.bam}/${base}.bam";
       got_bams=false
     } 
     if (! file("${params.bam}/${base}.bam.bai").exists())  {
       println "No file ${params.bam}/${base}.bam.bai";
       got_bams=false
     }
  }
  if (! got_bams) {
   println "You promised me bams but at least some of them or their index files missing"
   System.exit(18)
  }
  globs = fname .collect { "${it}*" }.join(",")
  bf = Channel.fromFilePairs("${params.bam}/{${globs}}",size:3) { file -> file.simpleName }.map { b, f -> [b,f[0],f[1]] }
  (bam_file1, bam_file3) = bf.into (2)
} else {

  input_ch = Channel.fromPath(fnames)



   // Create a BAM file from the reads aligning to the reference genome
   // This is done per sample
   process pbio_bamify {
     cpus params.bamify_cpus
     memory params.bamify_mem
     errorStrategy 'finish'
     input:
       path(fq) from input_ch
     output:
       set val(the_id), file("${the_id}.bam"), file("${the_id}.bam.bai") into (bam_file1, bam_file2, bam_file3)
     publishDir params.bam
     script:
       the_id = fq.baseName       
       """
       module load smrtlink
       ls $fq/*.fa*gz   > files.fofn
       pbmm2 align $ref_mmi  files.fofn ${the_id}.bam --preset HIFI --rg '@RG\tID:$the_id\tSM:$the_id' --unmapped --sort 
       """
   }
}

autosomes = 1..22
chroms    = autosomes +  ['X','Y','M'] 


process discover {

  input:
     set val(base), file(bam), file(bai) from bam_file3
     file(trf)
  output:
     file("${base}.svsig.gz") into svsig_ch
  memory 4.GB
  cpus 1
  errorStrategy 'finish'
  script:
     base = bam.simpleName
     """
     module load smrtlink
     pbsv discover --tandem-repeats  $trf $bam ${base}.svsig.gz
     """
}

process pbsvcall {
  input:
     file(sigs) from svsig_ch.collect()
     file(ref_seq)
  output:
     set file("${vcf}.gz"), file("${vcf}.gz.tbi") into pbs_call_ch
  memory "12GB"
  errorStrategy 'finish'
  publishDir "${params.vcf}/"
  script:
  vcf = "${params.out}.vcf"
     """
       module load smrtlink
       hostname
       /usr/bin/time -f "%e %M" pbsv call -j 8 --ccs  $ref_seq $sigs $vcf
       bgzip $vcf
       tabix ${vcf}.gz 
     """
}
