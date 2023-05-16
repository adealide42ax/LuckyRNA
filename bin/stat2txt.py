import glob
import pandas as pd
import re
import sys

filenames = glob.glob('*/2*.txt')

data = []

for filename in filenames:
    with open(filename, "r") as file:
        contents = file.read()
        lines = contents.splitlines()
        data.append({'filename': filename, 'gene number': len(lines)})
        
def process_filename(filename):
    parts = re.split(r'[/_]', filename)
    condition1, condition2, fc,fdr, _,change, method, _ = parts[-8:]
    return f"{condition1},{condition2},{fc},{fdr},{change},{method}"

def merge_columns(row):
    return row["Condition1"] + "_vs_" + row["Condition2"]

df = pd.DataFrame(data)
df['filename'] = df['filename'].apply(process_filename)
df[["Condition1", "Condition2", "FC", "FDR", "Change", "Method"]] = df["filename"].str.split(",", expand=True)
df["Compare"] = df.apply(merge_columns, axis=1)
columns_to_drop = ["filename", "Condition1", "Condition2"]
df = df.drop(columns_to_drop, axis=1)
df = df.reindex(columns=['Compare', 'FC','FDR','Change','Method','gene number'])

df.to_csv("stat_2_txt.csv", index=False, sep=",")