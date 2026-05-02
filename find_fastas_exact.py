# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:56:26 2020

@author: jeff
"""
import subprocess
import shutil
import os

path = '/data_store/seq_data/export_16S_sccoos/'

found = set(os.listdir(path))

## Note that I've gotten a little inconsistent with naming.  If project used by 16S and 18S
## need to be specific about what this script is looking for and where to put.

with open('sccoos_combined_16S.txt', 'r') as file_in, open('filesnotfound.txt', 'w') as not_found:
    for line in file_in:
        line = line.rstrip()
        if line not in found:

            temp = subprocess.Popen('locate -i ' + line,
                                    shell = True,
                                    executable = '/bin/bash',
                                    stdout = subprocess.PIPE)
            
            locations = temp.communicate()[0]
            locations = locations.decode().split('\n')
            fastas = [i for i in locations]
                            
            fasta  = fastas[0]
            file_name = fasta.split('/')[-1]
            
            try:
                shutil.copyfile(fasta, path + file_name) 
                found.add(file_name)
                print(file_name)
            except shutil.SameFileError:
                print(file_name, 'already present')
                    
            if len(fastas) == 0:
                print(line, 'not found')
                print(line, file = not_found)


