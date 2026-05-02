# -*- coding: utf-8 -*-
"""
Created on Thu Dec  4 12:42:44 2025

@author: jeff
"""

import glob
import pandas as pd
from pathlib import Path

files = glob.glob('high_qual_draft_eggnog/*annotations')

list_out = []
names = []

table_out = pd.DataFrame(index = [Path(f).name.split('.emapper.annotations')[0] for f in files])

for f in files:
    name = Path(f).name.split('.emapper.annotations')[0]
    annotation = pd.read_csv(f, sep = '\t', skiprows = 4, index_col = 0, skipfooter = 3, engine = 'python')
    temp = list(annotation.KEGG_ko)
    temp = [item.replace(',', '|') for item in temp]
    temp = list(set(temp))
                
    temp_df = pd.DataFrame(data = 1, index = [name], columns = temp)
    list_out.append(temp_df)
    print(name)
    if name in names:
        break
    else:
        names.append(name)
    
table_out = pd.concat(list_out)

table_out.to_csv('WATL_mg_combined_KO.csv')

