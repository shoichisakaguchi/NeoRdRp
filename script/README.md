## Creating RdRp HMM profiles by yourself

You can create your original RdRp HMM profiles by runnning the following script.

`python3 pipeline.py INPUT.FASTA`

## Recuirements
- [CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [Divvier](https://github.com/simonwhelan/Divvier)
- [HMMER](http://hmmer.org)
- [Pandas](https://pandas.pydata.org)
- [Biopython](https://biopython.org)

## Changing parameters

Open pipeline.py and change the parameters in the parameter section: you can change CD-HIT identity (default 0.6), wordsize (default 4), and minimum sequences in the cluster (default 3).

```
# parameter
param_cdhit_threshould = str(0.6)
param_cdhit_wordsize = str(4)
param_cdhit_cluster = str(3)
```

If you need to change the definition of boundary, open cutgap.py and change the parameter: you can change threshold (defaul 0.25) and minimum length (defaul 9).

```
threshold=0.25
minlength = 9
```
