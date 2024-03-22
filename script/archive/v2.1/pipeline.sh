#!/bin/sh

# Last update: 2024/02/26
# Shoichi Sakaguchi

#
# usage
#
usage () {
  echo "usage:"
  echo "${0} -i INPUT [-p program directory] [-h]"
  echo " "
  echo "CD-HIT, biopython are necessary to run this script."
  exit 1
}

#
# CD-HIT parameter
#

threshould=0.6
wordsize=4
cluster=3

echo "CD-HIT wordsize: $wordsize"
echo "CD-HIT minimum cluster size: $cluster"

#
# path
#

dir_program="."

#
# maximum threads
#

CPUs=$((`lscpu -p=CPU | wc -l`-4))

echo "Maximum threads used: $CPUs"

#
# getopts
#
FLG_I=false
FLG_O=false
FLG_H=false

while getopts i:o:p:h OPT
do
  case $OPT in
    "i" ) INPUT="$OPTARG" ; FLG_I=true ;;
    "o" ) OUTPUT="$OPTARG" ; FLG_O=true ;;
    "p" ) dir_program="$OPTARG" ;;
    "h" ) FLG_H=true ;;
  esac
done

#
# Show usage without input or with help option.
#

if ! "${FLG_I}" || "${FLG_H}"; then usage; fi

#
# Rename sequences of input fasta
#
echo "Renaming..."

python3 $dir_program/rename_full2num.py -i $INPUT;

#
# CD-HIT to remove duplicates
#
echo "Removing duplicates..."

cd-hit -i input.num.fasta -o input.num.99.cdhit -c 0.99 -n 5 -M 0 -d 0 -T 0

#
# CD-HIT to cluster RdRps
#
echo "Clustering..."

cd-hit -i input.num.99.cdhit -o input.num.cdhit -c $threshould -n $wordsize -M 0 -d 0 -T 0

echo "Finished clustering."

#
# Parse clusters containing more than N sequences
#

echo "Making each fasta from clusters."

make_multi_seq.pl input.num.fasta input.num.cdhit.clstr multi-seq $cluster

#
# Alignment
#

echo "Making alignments..."

for i in `ls multi-seq | grep -v .aln.fasta | grep -v cdhit | grep -v multi-seq`; do mafft --auto --thread -1 multi-seq/$i > multi-seq/${i}.aln.fasta; done

#
# split alignments: if 6 continued 1/4 gap in a region
#

echo "Splitting sequences by gappy region..."

python3 $dir_program/cutgap.py multi-seq/

#
# rename
#

echo "Renaming..."

find multi-seq/ -name "*RDRP*" | xargs -P $CPUs -IZZZ python3 $dir_program/rename_num2full.py multi-seq/ZZZ input.num-name.tsv multi-seq

#
# make HMM profiles
#

echo "Making HMM profiles..."

find multi-seq/ -name "*.full.fasta" | xargs -P $CPUs -IYYY hmmbuild YYY.hmm YYY

echo "Merging HMM profiles..."

awk 1 multi-seq/*.hmm > input.hmm

echo "Doing hmmpress..."

hmmpress input.hmm

echo "Process completed."