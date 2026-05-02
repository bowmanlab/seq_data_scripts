# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:56:26 2020

@author: jeff
"""
import subprocess
import shutil
import os

path = '/data_store/seq_data/export_18S_sccoos/'

found = set(os.listdir(path))

## Note that I've gotten a little inconsistent with naming.  If project used by 16S and 18S
## need to be specific about what this script is looking for and where to put.

with open('files2find.txt', 'r') as file_in, open('filesnotfound.txt', 'w') as not_found:
    for line in file_in:
        line = line.rstrip()
        if line not in found:

            temp = subprocess.Popen('locate -i ' + line + '|grep -v \'filt\'|grep -v \'cutadapt\'|grep \'18S\'|grep \'fastq\'|grep -v \'v4\'',
                                    shell = True,
                                    executable = '/bin/bash',
                                    stdout = subprocess.PIPE)
            
            locations = temp.communicate()[0]
            locations = locations.decode().split('\n')
            fastas_R1 = [i for i in locations if 'R1' in i] 
            fastas_R2 = [i for i in locations if 'R2' in i]
            
            try:
                fastas = [fastas_R1[0]] + [fastas_R2[0]]
                                
                for fasta in fastas:
                    file_name = fasta.split('/')[-1]
                    
                    try:
                        shutil.copyfile(fasta, path + file_name) 
                        found.add(file_name)
                        print(fasta)
                    except shutil.SameFileError:
                        print(file_name, 'already present')
            except IndexError:
                print(line, 'not found')
