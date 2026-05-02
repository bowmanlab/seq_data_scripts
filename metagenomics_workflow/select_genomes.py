# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 11:41:07 2025

@author: jeff
"""

import shutil
import os
import pandas as pd

if os.path.exists('high_qual_draft'):
    shutil.rmtree('high_qual_draft')
    
if os.path.exists('high_qual_draft_prodigal'):
    shutil.rmtree('high_qual_draft_prodigal')
    
os.mkdir('high_qual_draft')
os.mkdir('high_qual_draft_prodigal')

data = pd.read_csv('drep_output/data_tables/Widb.csv', index_col = 0)
data_hq = data.loc[data['completeness'] > 50]
data_hq = data_hq.loc[data_hq['contamination'] < 10]

for mag in data_hq.index:
    name = mag.split('.fa')[0]
    shutil.copy('drep_output/dereplicated_genomes/' + mag, 'high_qual_draft/' + mag)
    shutil.copy('drep_output/data/prodigal/' + mag + '.faa', 'high_qual_draft_prodigal/' + name + '.faa')
    shutil.copy('drep_output/data/prodigal/' + mag + '.fna', 'high_qual_draft_prodigal/' + name + '.fna')
    
data_hq.to_csv('hiqh_qual_draft_info.csv')
