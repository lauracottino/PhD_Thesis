#!/usr/bin/env nextflow

reports = Channel.fromPath("/*.csv")
prefix = "/SplitReports/"
suffix = ".txt"
script = "/opt/exp_soft/bioinf/PennCNV/split_illumina_report_NEW.pl"

process runSplitIlluminaReport {
	maxForks 10

   	 input:
    	 file(report) from reports 

   	 output:
    	 
	 """
	 ${script} -prefix ${prefix} -suffix ${suffix} -c ${report}
	 """
}
