# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 07:32:12 2025

@author: jeff
"""

import pandas as pd
import os

n_reads_df = pd.read_csv('sample_read_counts.csv', index_col = 0)
bin_info = pd.read_csv('hiqh_qual_draft_info.csv', index_col = 0)

col_names = bin_info.index
        
data_rpm = pd.DataFrame(columns = col_names)

columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'MRNM', 'MPOS', 'ISIZE', 'SEQ', 'QUAL'] 
        
for f in os.listdir('mag_map'):
    if f.endswith('.sam.gz'):
        print(f, 'calculating RPM')
        name = f.split('_map.q20.sam.gz')[0]
        n_reads = n_reads_df.loc[name, 'n_reads']
        sam = pd.read_csv('mag_map/' + f, sep='\s+', usecols = list(range(0,11)), names = columns,  on_bad_lines = 'skip')
        sam.drop_duplicates(subset = 'QNAME', inplace = True)
        bin_hits = sam.RNAME.str.split('|', expand = True)
        bin_hits.columns = ['bin', 'contig']
        pmf = n_reads/1000000
        bin_tally = bin_hits.bin.value_counts()/pmf
        bin_tally.name = name
        data_rpm = pd.concat([data_rpm, bin_tally.to_frame().T])

data_rpm.to_csv('WATL_mg_combined_RPM.csv')