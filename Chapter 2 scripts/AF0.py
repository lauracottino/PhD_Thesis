import pandas as pd
import sys
fname = sys.argv[1]

df = pd.read_csv(fname,delimiter="\t")
want = df[df['AF']==0]
to_print = (want[['CHROM','POS']])
to_print.to_csv('/spaces/emma/MantaResults/AF0_to_exclude.tsv',sep='\t')
#print(df)
