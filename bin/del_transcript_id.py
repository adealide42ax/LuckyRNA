import pandas as pd
import numpy as np
from pandas import DataFrame
import sys
path = sys.argv[1]
df = pd.read_csv(path,sep = '\t')
del df['transcript_id(s)']
df = df.set_index(df.columns[0])
df = df.round(0)

for col in df.columns:
    if col == 'gene_id':
        pass
    else:
        df[col] = df[col].astype("int")

df.to_csv(path, sep='\t')
