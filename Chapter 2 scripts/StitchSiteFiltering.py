import pandas as pd
import sys
fname = sys.argv[1]
template_vcf=sys.argv[2]
output_vcf = sys.argv[3]
df = pd.read_csv(fname,delimiter="\t")
want = df[df['stitch_site_already_gt']==False]
# pull the prefatory part of the VCF out                                                                                             
template=""
f = open(template_vcf)
for line in f:
    if line[0] != "#": break
    template = template+line
fout = open(output_vcf,"w")
fout.write(template)
for i, cnv in want.iterrows():
     fout.write("\t".join([cnv['chrom'],str(cnv['start']),cnv['stitched_cnv_site_ID'],'A'\
,'<CNV>','.','PASS','END=%d;SVTYPE=CNV'%cnv['end']])+"\n")
