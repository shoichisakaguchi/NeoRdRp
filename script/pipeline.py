#!/usr/bin/env python

# NeoRdRp pipeline
# Last update: 2022/03/24
# Shoichi Sakaguchi

import csv, sys
import argparse
import subprocess as sp
from subprocess import PIPE
import glob

# parser
parser = argparse.ArgumentParser(description='Process RdRp sequences into HMM profiles.')
parser.add_argument("-i", "--input", required=True, help="Input file name, required. Expecting amino acid sequence of RdRps in FASTA format.")
parser.add_argument("-o", "--output", default="1st_round.aa.rdrp.hmm", help="Output filename, option. Extension should be '.hmm'.")
args = parser.parse_args()

# location of files
fasta_1 = args.input

# parameter
param_cdhit_threshould = str(0.6)
param_cdhit_wordsize = str(4)
param_cdhit_cluster = str(3)

# print parameter
print("Making HMM profiles from", fasta_1)
print("CD-HIT wordsize:", param_cdhit_wordsize)
print("CD-HIT minimum cluster size:", param_cdhit_cluster)

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
            
out_fasta2 = fasta_1.replace("fasta", "").replace("fa", "") + "num.fasta"

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
out_rename2 = fasta_1.replace("fasta", "").replace("fa", "") + "num-name.tsv"
save_list = []
with open(out_rename2, "w") as f:
    for k, v in d_acc_num.items():
        l = [v, k]
        save_list.append(l)
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(save_list)

print("Renamed input file and created acc-num list.")

### Clustering

print("Start clustering...")

cdhit_out = out_fasta2.replace(".fasta", ".cd-hit.fasta")

res = sp.run(
    [
        "cd-hit",
        "-i", 
        out_fasta2,
        "-o",
        cdhit_out,
        "-c",
        param_cdhit_threshould,
        "-n",
        param_cdhit_wordsize,
        "-M",
        "0",
        "-d",
        "0",
        "-T",
        "0"
    ],
    stdout=PIPE, stderr=PIPE, text=True)

print(res.stderr)

print("Finished clustering.")

### Clusters containing more than N sequences

print("Making each fasta from clusters.")

cluster_out = cdhit_out + ".clstr"

res = sp.run(
    [
        "make_multi_seq.pl",
        out_fasta2,
        cluster_out,
        "1st_round_multi-seq",
        param_cdhit_cluster
    ],
    stdout=PIPE, stderr=PIPE, text=True)

### Alignment

print("Making alignments.")

files = glob.glob("1st_round_multi-seq/*")
for file in files:
    out = file + ".aln.fasta"

    res = sp.run(
        [
        "mafft-linsi", 
        "--thread",
        "-1",
        file
        ],
    stdout=PIPE, stderr=PIPE, text=True)

    with open(out, "w") as f:
        f.write(res.stdout)

### Replace uncertain amino acids with gaps

print("Divvying...")

files = glob.glob("1st_round_multi-seq/*.aln.fasta")
for file in files:
    res = sp.run(
        [
        "divvier", 
        "-mincol",
        "3",
        "-divvygap",
        file
        ],
    stdout=PIPE, stderr=PIPE, text=True)

### split alignments: if 6 continued 1/4 gap in a region

print("Splitting sequences by gappy region.")

res = sp.run(
    [
        "python3",
        "cutgap.py",
        "1st_round_multi-seq/"
    ]
)

print(res.stdout)
print(res.stderr)

### rename

print("Renaming...")

files = glob.glob("1st_round_multi-seq/*RDRP*")
for file in files:
    res = sp.run(
        [
            "python3",
            "rename_num2full.py",
            file,
            out_rename2,
            "1st_round_multi-seq"
        ]
    )

print("Making HMM profiles.")

files = glob.glob("1st_round_multi-seq/*.full.fasta")
for file in files:
    out = file + ".hmm"

    res = sp.run(
        [
        "hmmbuild", 
        out,
        file,
        ],
    stdout=PIPE, stderr=PIPE, text=True)

print("Merging HMM profiles...")

path = "1st_round_multi-seq/"
files = glob.glob(path + "*.hmm")

with open(args.output, "wb") as f_new:
    for f in files:
        print(f)
        with open(f, "rb") as f_org:
            f_new.write(f_org.read())

print("Doing hmmpress...")

res = sp.run(
    [
        "hmmpress",
        args.output
    ],
    stdout=PIPE, stderr=PIPE, text=True
    )

print("Process completed.")