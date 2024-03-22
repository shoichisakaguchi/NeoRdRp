import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_regions(records, threshold=0.25, min_skips=8):
    
    df_aas = pd.DataFrame( [ list( f.seq ) for f in records] ) 
    df_aas.columns += 1 
    df_aas["z"] = "-"

    # sr_gaprate = df_aas.apply( lambda x:x.value_counts()['-'] , axis=0) / len( df_aas)
    sr_gaprate = df_aas.apply( lambda x:x.value_counts() ).fillna(0).loc['-'] / len( df_aas)
    targets = sr_gaprate.loc[ sr_gaprate <= threshold ].index.to_frame(name='pos')
    targets['diff'] = targets.diff()

    if len( targets) == 0:
        print('warning : target is not found !')
        return pd.DataFrame({'start':[], 'end':[]})

    return pd.DataFrame({ 
        'start' : [ int(targets.iloc[0]['pos'])]  
                + list( targets.loc[ targets['diff'] >= min_skips,'pos']),
        'end'   : list( targets.loc[ (targets['diff'] >= min_skips).shift(-1, fill_value=False),'pos' ])
                + [ int(targets.iloc[-1]['pos'])] 
    })

def get_seqs(records, start, end):
    return [
        SeqRecord( 
            i.seq[start:end], 
            id='%s_%d-%d'%( i.id, start, end), 
            description='dropped_gaps_by_hiroumauma' )
        for i in records
    ]

import glob
import os
import sys

threshold=0.25
datadir = sys.argv[1] # './data/'
files = datadir + "*.aln.fasta"
minlength = 9

for file in glob.glob(files):
    print('processing ... %s' % file)
    records  = list( SeqIO.parse(file, "fasta") )

    for i,(start,end) in get_regions(records, threshold=threshold).iterrows():
        if end - start >= minlength - 1:
            SeqIO.write(
                get_seqs( records, start, end),
                "%s%s_RDRP_%.2f_%d-%d.fasta" % (
                    datadir,
                    os.path.splitext( os.path.basename(file) )[0],
                    threshold, start, end
                ),
                'fasta'
            )