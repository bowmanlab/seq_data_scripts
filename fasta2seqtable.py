# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 16:23:22 2025

@author: jeff
"""

from Bio import SeqIO
import os
import pandas as pd

seqtable = pd.DataFrame()

for f in os.listdir('.'):
    if f.endswith('exp.fasta'):
        print(f)
        seqseries = pd.Series()
        name = f[:-10]
        seqseries.name = name
        for record in SeqIO.parse(f, 'fasta'):
            try:
                seqseries[str(record.seq)] = seqseries[str(record.seq)] + 1
            except KeyError:
                seqseries[str(record.seq)] = 1
                
        seqtable = pd.concat([seqtable, seqseries], axis = 1)

seqtable.fillna(0, inplace = True)

seqtable = seqtable.transpose()            

seqtable.to_csv('18S_R1_seqtable.csv')
