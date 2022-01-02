## NeoRdRp

NeoRdRp is a comprehensive dataset of RNA-dependent RNA polymerase (RdRp) domains.

NeoRdRp consists of NCBI RefSeq and NCBI GenBank sequence collections updated as new RdRp are added. The files contained in this directory include FASTA files and HMM profiles and annotation data for NeoRdRp. The paper describing NeoRdRp is published (Sakaguchi, Shoichi, Syun-Ichi Urayama, Yoshihiro Takaki, Hong Wu, Youichi Suzuki, Takuro Nunoura, Takashi Nakano, and So Nakagawa. 2022. “NeoRdRp: A Comprehensive Dataset for Identifying RNA-Dependent RNA Polymerase of Various RNA Viruses from Metatranscriptomic Data.” bioRxiv. https://doi.org/10.1101/2021.12.31.474423).

## Download
| NeoRdRp version | HMM profiles     | Amino acid sequences | Annotation table         | Score table |
|-----------------|------------------|----------------------|--------------------------|-------------|
| 1.0 <br> (2021-12-27) | [NeoRdRp-HMM.v1 <br> (2,234 profiles)](https://github.com/shoichisakaguchi/NeoRdRp/blob/main/NeoRdRp-HMM.v1.hmm.gz)  | [NeoRdRp-seq.v1  <br> (12,316 sequences)](https://github.com/shoichisakaguchi/NeoRdRp/blob/main/NeoRdRp-seq.v1.fasta.gz) | [NeoRdRp-annotation.v1.txt](https://github.com/shoichisakaguchi/NeoRdRp/blob/main/NeoRdRp_annotation.v1.txt) | [NeoRdRp-score.v1.xlsx](https://github.com/shoichisakaguchi/NeoRdRp/blob/main/NeoRdRp-score.v1.xlsx) |

## NeoRdRp-HMM.vX.hmm.gz
This file is the HMM profile created from the RdRp collection collected by NeoRdRp. Specify the unzipped NeoRdRp-HMM file as the hmmsearch option. An example of the command is shown as follows:

```hmmsearch --domtblout OUTPUT.domtblout --cpu 20 -E 1e-10 NeoRdRp-HMM.vX.hmm AMINO_ACID.FASTA```

## NeoRdRp-seq.vX.fasta.gz
This file is the amino acid sequence of the RdRp collection of NeoRdRp. By adding any RdRp sequence to this FASTA file, NeoRdRp can be trained.

## NeoRdRp_annotation.vX.txt
This file shows the results of annotation and score of the HMM profiles used in the UniProtKB search with NeoRdRp. Detected RdRp of interest can be verified by referring to the annotation information of the HMM profile.

## NeoRdRp-score.vX.xlsx
This file shows the following information for each HMM profile: the number and length of amino acid sequences used to build the HMM profile, the number of times they were used in the UniProtKB search, the mean score, and the standard deviation.
