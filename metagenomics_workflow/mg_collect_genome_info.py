# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 08:43:15 2026

@author: jeff
"""

import pandas as pd
import os

bins = os.listdir('high_qual_draft')
bins = [item for item in bins if not item.endswith('.faa')]

data = pd.read_csv('checkm2_output/drep_compatible_genomeInfo.csv', index_col = 0)
data = data.loc[data.index.isin(bins)]
data.index = data.index.str.split('.fa').str[0]

class_ar = pd.read_csv('high_qual_classify_output/gtdbtk.ar53.summary.tsv', index_col = 0, sep = '\t')

data.loc[class_ar.index, 'classification'] = class_ar['classification']
data.loc[class_ar.index, 'closest_genome_reference'] = class_ar['closest_genome_reference']
data.loc[class_ar.index, 'closest_genome_ani'] = class_ar['closest_genome_ani']

class_ba = pd.read_csv('high_qual_classify_output/gtdbtk.bac120.summary.tsv', index_col = 0, sep = '\t')

data.loc[class_ba.index, 'classification'] = class_ba['classification']
data.loc[class_ba.index, 'closest_genome_reference'] = class_ba['closest_genome_reference']
data.loc[class_ba.index, 'closest_genome_ani'] = class_ba['closest_genome_ani']

data.to_csv('high_qual_draft_info.csv')

