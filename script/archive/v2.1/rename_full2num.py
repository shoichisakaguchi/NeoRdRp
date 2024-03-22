#!/usr/bin/env python

# Last update: 2023/06/16
# Shoichi Sakaguchi

import csv
import argparse

# parser
parser = argparse.ArgumentParser(description='Process RdRp sequences into HMM profiles.')
parser.add_argument("-i", "--input", required=True, help="Input file name, required. Expecting amino acid sequence of RdRps in FASTA format.")
args = parser.parse_args()

# location of files
fasta_1 = args.input

### Rename sequences of input fasta
d = {} #{acc:seq}
d_acc_num = {} #{acc:num}
i = 1
with open(fasta_1) as f:
    for line in f:
        if line.startswith(">"):
            acc = line.rstrip().replace(">", "")
            seq = ""
            d_acc_num[acc] = i
            i += 1
        else:
            seq = seq + line.rstrip()
            d[acc] = seq

# out_fasta2 = fasta_1.replace("fasta", "").replace("fa", "") + "num.fasta"
out_fasta2 = "input.num.fasta"

save_list = []
with open(out_fasta2, "w") as f:
    i = 1
    for k, v in d.items():
        name = ">" + str(i)
        l = [name, v]
        save_list.append(l)
        i += 1
    writer = csv.writer(f, delimiter="\n")
    writer.writerows(save_list)

### save acc num tsv
# out_rename2 = fasta_1.replace("fasta", "").replace("fa", "") + "num-name.tsv"
out_rename2 = "input.num-name.tsv"

save_list = []
with open(out_rename2, "w") as f:
    for k, v in d_acc_num.items():
        l = [v, k]
        save_list.append(l)
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(save_list)

print("Renamed input file and created full_name-(tab)-number list.")