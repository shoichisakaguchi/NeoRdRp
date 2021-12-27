## NeoRdRp

NeoRdRp is a comprehensive dataset of RNA-dependent RNA polymerase (RdRp) domains.

NeoRdRp consists of NCBI RefSeq and NCBI GenBank sequence collections updated as new RdRp are added. The files contained in this directory include FASTA files and HMM profiles and annotation data for NeoRdRp. A paper describing NeoRdRp is currently under submission.

## Download
| NeoRdRp version | HMM profiles     | Amino acid sequences | Annotation table         |
|-----------------|------------------|----------------------|--------------------------|
| 1.0 <br> (2021-12-27) | [NeoRdRp-HMM.v1 <br> (2,234 profiles)](https://github.com/shoichisakaguchi/NeoRdRp/blob/main/NeoRdRp-HMM.v1.hmm.gz)  | [NeoRdRp-seq.v1  <br> (12,316 sequences)](https://github.com/shoichisakaguchi/NeoRdRp/blob/main/NeoRdRp-seq.v1.fasta.gz) | [NeoRdRp-annotation.v1.txt](https://github.com/shoichisakaguchi/NeoRdRp/blob/main/NeoRdRp_annotation.v1.txt) |

## NeoRdRp-HMM.vX.hmm.gz
This file is the HMM profile created from the RdRp collection collected by NeoRdRp. Specify the unzipped NeoRdRp-HMM file as the hmmsearch option. An example of the command is shown as follows:

```hmmsearch --domtblout OUTPUT.domtblout --cpu 20 -E 1e-10 NeoRdRp-HMM.vX.hmm AMINO_ACID.FASTA```

## NeoRdRp-seq.vX.fasta.gz
This file is the amino acid sequence of the RdRp collection of NeoRdRp. By adding any RdRp sequence to this FASTA file, NeoRdRp can be trained.

## NeoRdRp_annotation.vX.txt
This file shows the results of annotation and score of the HMM profiles used in the UniProtKB search with NeoRdRp. Detected RdRp of interest can be verified by referring to the annotation information of the HMM profile.
