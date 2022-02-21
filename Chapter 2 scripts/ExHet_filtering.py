import pandas as pd
import sys
fname = sys.argv[1]

df = pd.read_csv(fname,delimiter="\t")
want = df[df['P_HET_EXCESS']<0.00001]
to_print = (want[['CHR', 'POS']])
to_print.to_csv('ExHet_to_exclude.tsv',sep='\t')
