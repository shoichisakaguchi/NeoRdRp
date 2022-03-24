# %%
import os, sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# %%
file_seq = sys.argv[1] #"output/16.alnta.divvy_RDRP_0.25_1081-1104.fasta"
file_list = sys.argv[2] #"original/merged-seq.acc-num.tsv"
dir = sys.argv[3] # "data/"

df_name = pd.read_csv(file_list, sep="\t", header=None)
df_name.columns = ["num", "seq_name"]

df_name["num"] = df_name["num"].astype(str)

l = []
for record in SeqIO.parse(file_seq, "fasta"):
    num = record.id.split("_")[0]
    pos = record.id.split("_")[1]
    seq = record.seq
    l.append([num, pos, seq])

df_seq = pd.DataFrame(l, columns=["num", "pos", "seq"])

df_rename = pd.merge(
    df_name,
    df_seq,
    on="num",
    how="inner"
)

# 拡張子なしのファイル名を取得  
out = dir + "/" + os.path.splitext(os.path.basename(file_seq))[0] + ".full.fasta"

l_save = []
for i in df_rename.iterrows():
    seq_r = SeqRecord(i[1].seq, id = i[1].seq_name + " " + i[1].pos)
    l_save.append(seq_r)

SeqIO.write(l_save, out, "fasta")

# %%



