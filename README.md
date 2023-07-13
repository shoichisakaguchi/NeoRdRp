## NeoRdRp

NeoRdRp is a comprehensive dataset of RNA-dependent RNA polymerase (RdRp) domains.

NeoRdRp consists of NCBI RefSeq and NCBI GenBank sequence collections updated as new RdRp are added. The files contained in this directory include FASTA files and HMM profiles and annotation data for NeoRdRp. 

## Download
| NeoRdRp version | HMM profiles     | Amino acid sequences | Annotation table         | Score table | Multiple Sequence Alignments |
|-----------------|------------------|----------------------|--------------------------|-------------|------------------------------|
|1.1<br>(2022-03-24)|[NeoRdRp-HMM.v1.1<br>(1,182 profiles)](https://github.com/shoichisakaguchi/NeoRdRp/blob/fe1070534804d9c9f9725421538cd92390a78ea8/NeoRdRp-HMM.v1.1.hmm.gz)|[NeoRdRp-seq.v1.1<br>(12,502 sequences)](https://github.com/shoichisakaguchi/NeoRdRp/blob/9f0e4f0ffdd8ae83ab8c333fb9a0d85fa0c0f39b/NeoRdRp-seq.v1.1.fasta.gz)|[NeoRdRp-annotation.v1.1.txt](https://github.com/shoichisakaguchi/NeoRdRp/blob/6214aeaf1f4566783c061ee4555c12671fe7404e/NeoRdRp-annotation.v1.1.txt)|[NeoRdRp-score.v1.1.xlsx](https://github.com/shoichisakaguchi/NeoRdRp/blob/629ff302ffeb75cf2e47449aef64e43647caa213/NeoRdRp-score.v1.1.xlsx)|[NeoRdRp-msa.v1.1.tar.gz](https://github.com/shoichisakaguchi/NeoRdRp/blob/main/archive/v1.1/NeoRdRp-msa.v1.1.tar.gz)|
| 1.0* <br> (2021-12-27) | [NeoRdRp-HMM.v1 <br> (2,234 profiles)](https://github.com/shoichisakaguchi/NeoRdRp/blob/9b7aa7e258b91866cb531d3184acdc79f3b9d6dc/archive/v1.0/NeoRdRp-HMM.v1.hmm.gz)  | [NeoRdRp-seq.v1  <br> (12,316 sequences)](https://github.com/shoichisakaguchi/NeoRdRp/blob/9b7aa7e258b91866cb531d3184acdc79f3b9d6dc/archive/v1.0/NeoRdRp-seq.v1.fasta.gz) | [NeoRdRp-annotation.v1.txt](https://github.com/shoichisakaguchi/NeoRdRp/blob/9b7aa7e258b91866cb531d3184acdc79f3b9d6dc/archive/v1.0/NeoRdRp_annotation.v1.txt) | [NeoRdRp-score.v1.xlsx](https://github.com/shoichisakaguchi/NeoRdRp/blob/9b7aa7e258b91866cb531d3184acdc79f3b9d6dc/archive/v1.0/NeoRdRp-score.v1.xlsx) | |

*Note that there was a small bug in the construction of the HMM profiles, which has been fixed in version 1.1. In addition, some parameters including the minimum set of sequences usd for the HMM profile construction (verion 1: 5 -> version 1.1: 3) as well as length for split of ambigous regions (verion 1: 1 -> verion 1.1: 5) were changed.

## NeoRdRp-HMM.vX.hmm.gz
This file is the HMM profile created from the RdRp collection collected by NeoRdRp. Specify the unzipped NeoRdRp-HMM file as the hmmsearch option. An example of the command is shown as follows:

```hmmsearch --domtblout OUTPUT.domtblout --cpu 20 -E 1e-10 NeoRdRp-HMM.vX.hmm AMINO_ACID.FASTA```

## NeoRdRp-seq.vX.fasta.gz
This file is the amino acid sequence of the RdRp collection of NeoRdRp. By adding any RdRp sequence to this FASTA file, NeoRdRp can be trained.

## NeoRdRp_annotation.vX.txt
This file shows the results of annotation and score of the HMM profiles used in the UniProtKB search with NeoRdRp. Detected RdRp of interest can be verified by referring to the annotation information of the HMM profile.

## NeoRdRp-score.vX.xlsx
This file shows the following information for each HMM profile: the number and length of amino acid sequences used to build the HMM profile, the number of times they were used in the UniProtKB search, the mean score, and the standard deviation.

## Making RdRp HMM profiles by yourself
The pipeline script and related scripts are located in the script directory; see the README documentation in the script directory.

## Reference
The details of this dataset were described at the following paper available at the Microbes and Environments website: <br>
Sakaguchi, S., Urayama, S.-I., Takaki, Y., Hirosuna, K., Wu, H., Suzuki, Y., Nunoura, T. and Nakagawa, S. (2022). NeoRdRp: A Comprehensive Dataset for Identifying RNA-dependent RNA Polymerases of Various RNA Viruses from Metatranscriptomic Data. <i>Microbes Environ</i>. 37, ME22001. doi: 10.1264/jsme2.ME22001. <br>
https://www.jstage.jst.go.jp/article/jsme2/37/3/37_ME22001/_article
