# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 16:18:15 2025

@author: jeff
"""

from Bio import SeqIO
import os

records = []

for f in os.listdir('high_qual_draft'):
    print(f)
    for record in SeqIO.parse('high_qual_draft/' + f, 'fasta'):
        record.id = f + '|' + record.id
        records.append(record)
        
with open('high_qual_draft.fasta', 'w') as fasta_out:
    SeqIO.write(records, fasta_out, 'fasta')
    