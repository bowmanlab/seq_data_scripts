#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 07:56:57 2023

@author: jeff
"""


import pandas as pd
import sys

csv_in = sys.argv[1]

seqtable = pd.read_csv(csv_in, index_col = 0)

for col in seqtable.columns:
    i = 0
    with open(col + '.exp.fasta', 'w') as fasta_out:
        for seq in seqtable.index:
            i = i + 1
            n = seqtable.loc[seq, col]
            
            for j in range(0,n):
                print('>' + str(i) + '_' + str(j) + '\n' + str(seq), file = fasta_out)
                
    print(col)