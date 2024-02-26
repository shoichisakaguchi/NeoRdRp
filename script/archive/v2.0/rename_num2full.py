import os, sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

file_seq = os.path.basename(sys.argv[1]) 
file_list = sys.argv[2] 
dir = sys.argv[3] 

df_name = pd.read_csv(file_list, sep="\t", header=None)
df_name.columns = ["num", "seq_name"]

df_name["num"] = df_name["num"].astype(str)

l = []
for record in SeqIO.parse(dir + "/" + file_seq, "fasta"):
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

# Get file name without extension  
out = dir + "/" + os.path.splitext(os.path.basename(file_seq))[0] + ".full.fasta"

l_save = []
for i in df_rename.iterrows():
    seq_r = SeqRecord(i[1].seq, id = str(i[1].seq_name) + " " + str(i[1].pos) )
    l_save.append(seq_r)

SeqIO.write(l_save, out, "fasta")