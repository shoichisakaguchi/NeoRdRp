## NeoRdRp

## Introduction
NeoRdRp is a dataset of Hidden Markov Model (HMM) profiles of RNA-dependent RNA polymerase (RdRp) domains. The RdRp gene is essential in almost all RNA viruses, excluding retroviruses. Therefore, The NeoRdRp can be used to search for various RNA viruses from metatranscriptome sequencing data.

## What's New
### Version 2.1
[Download NeoRdRp 2.1](https://doi.org/10.5281/zenodo.10851672)
- Threshold for splitting a gapped region was changed from 6 to 8, resulting in better performance (fewer non-specific hits) of UniProtKB.
- HMM Database: Now includes 19,394 HMM profiles from 557,197 RdRp-containing sequences.

### Version 2.0
[Download NeoRdRp 2.0](https://doi.org/10.5281/zenodo.10807561)
- Expanded Database: Now includes 21,025 HMM profiles from 557,197 RdRp-containing sequences.
- Improved Accuracy: Recall and specificity rates improved to 99.4% and 81.6% using UniProtKB datasets.
- Enhanced Detection: Outperforms other RdRp search tools in balanced detection.

### Version 1.1
[Download NeoRdRp 1.1](https://github.com/shoichisakaguchi/NeoRdRp/blob/fe1070534804d9c9f9725421538cd92390a78ea8/NeoRdRp-HMM.v1.1.hmm.gz)
- Initial Release: 1,182 HMM profiles created from 12,502 RdRp-containing  sequences.
- Refinement: 426 initial HMM profiles refined for better accuracy.
- Performance: Successfully identified the majority of RNA viruses in the UniProtKB dataset.

### Version 1.0
[Download NeoRdRp 1.0](https://github.com/shoichisakaguchi/NeoRdRp/blob/9b7aa7e258b91866cb531d3184acdc79f3b9d6dc/archive/v1.0/NeoRdRp-HMM.v1.hmm.gz)

*Note that there was a small bug in the construction of the HMM profiles, which has been fixed in version 1.1. In addition, some parameters, including the minimum set of sequences used for the HMM profile construction (5 -> 3) as well as length for a split of ambiguous regions (1 -> 3), were changed.

## Future Updates
As NeoRdRp continues to evolve, new versions will be documented here, detailing enhancements and new features.

## HMMER Installation and Usage

### Installing HMMER

HMMER is required to utilize the HMM profiles in NeoRdRp. Follow these steps to install HMMER:

1. Visit the [HMMER website](http://hmmer.org/) and download the appropriate version for your system.
2. Extract the downloaded file.
3. Navigate to the extracted directory and compile HMMER by running:
   ```
   ./configure
   make
   make install
   ```
   (Note: You might need administrative privileges to install.)

### Using `hmmsearch` with NeoRdRp
After installing HMMER, you can use `hmmsearch` to identify RNA viruses using NeoRdRp HMM profiles:

1. Run `hmmsearch` with a NeoRdRp HMM profile against your sequence data. For example:
   ```
   hmmsearch [options] neordrp_profile.hmm your_sequence_data.fasta
   ```
2. Interpret the output based on `hmmsearch` documentation and your research requirements.

For more detailed instructions and options, refer to the [HMMER user guide](http://hmmer.org/documentation.html).

## Features
- Comprehensive RdRp dataset.
- High accuracy and specificity in RNA virus identification.
- User-friendly for diverse metatranscriptome data analysis.

## Contributing to NeoRdRp

We welcome contributions to the NeoRdRp dataset. If you have compiled a dataset of RdRp sequences, please consider sharing it with us. We will evaluate its integration into NeoRdRp, and if it enhances the performance, we will release a new version (incremented by 0.1) incorporating your data.

### Contribution Process:
1. **Send Your Dataset**: Share your RdRp dataset with us.
2. **Evaluation**: We will review and merge your dataset after thorough evaluation.
3. **Version Update**: If your contribution improves NeoRdRp, we will release an updated version.

For detailed instructions on how to contribute and the dataset requirements, please contact us. The contact details can be found in the [Contact](#contact) section.

## License
Released under the MIT License.

## Contact

For inquiries, support, or further information about NeoRdRp, please contact us:

- **Email**: shoichi.sakaguchi at ompu.ac.jp
- **X (Formerly Twitter)**: [@sho1sakaguchi](https://x.com/sho1sakaguchi?t=FLQLjUY28GH40V2ckTwy3g&s=09)

## Acknowledgments

We extend our heartfelt thanks to the RdRp summit, an initiative focused on building consensus in RdRp studies. Their efforts and collaboration have been invaluable to our work. 

We also want to express our gratitude to all those who have supported us throughout this project. Their assistance, insights, and encouragement have played a crucial role in the development and success of NeoRdRp.

## Reference
### NeoRdRp2
Sakaguchi S, Nakano T and Nakagawa S (2024) NeoRdRp2 with improved seed data, annotations, and scoring. <i>Front. Virol.</i> 4:1378695. doi: 10.3389/fviro.2024.1378695
https://www.frontiersin.org/articles/10.3389/fviro.2024.1378695/full

### NeoRdRp1
Sakaguchi, S., Urayama, S.-I., Takaki, Y., Hirosuna, K., Wu, H., Suzuki, Y., Nunoura, T. and Nakagawa, S. (2022). NeoRdRp: A Comprehensive Dataset for Identifying RNA-dependent RNA Polymerases of Various RNA Viruses from Metatranscriptomic Data. <i>Microbes Environ</i>. 37, ME22001. doi: 10.1264/jsme2.ME22001.
https://www.jstage.jst.go.jp/article/jsme2/37/3/37_ME22001/_article
