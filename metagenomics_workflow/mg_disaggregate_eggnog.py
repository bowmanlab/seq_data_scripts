# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:03:23 2026

@author: jeff
"""

import pandas as pd

data = pd.read_csv('high_qual_draft_eggnog/combined.emapper.annotations', index_col = 0, sep =  '\t', comment = '#', header = 0)

data['bin'] = data.index.str.split('|').str[0]

# Group by bin and export each group as TSV
for bin_name, group in data.groupby('bin'):
    # Create filename (replace any problematic characters)
    filename = f"high_qual_draft_eggnog/{bin_name}.emapper.annotations"
    
    # Export group as TSV
    group.to_csv(filename, sep='\t', index=False)
    print(f"Exported {filename} with {len(group)} rows")
