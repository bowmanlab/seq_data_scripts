# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:21:17 2026

@author: jeff
"""

import os
from pathlib import Path
import shutil

shutil.rmtree("high_qual_draft", ignore_errors=True); Path("high_qual_draft").mkdir()

for f in os.listdir('drep_output/dereplicated_genomes'):
    name = ".".join(f.split('.')[:-1])
    name_faa = name + '.faa'
    
    with open('drep_output/dereplicated_genomes/' + f, 'r') as fna_in, open('high_qual_draft/' + f, 'w') as fna_out:
        for line in fna_in:
            line = line.rstrip()
            if line.startswith('>'):
                line = line.strip('>')
                line = '>' + name + '|' + line
                print(line, file = fna_out)
            else:
                print(line, file = fna_out)
                
    with open('checkm2_output/protein_files/' + name_faa, 'r') as faa_in, open('high_qual_draft/' + name_faa, 'w') as faa_out:
        for line in faa_in:
            line = line.rstrip()
            if line.startswith('>'):
                line = line.strip('>')
                line = '>' + name + '|' + line
                print(line, file = faa_out)
            else:
                print(line, file = faa_out)
    
    
