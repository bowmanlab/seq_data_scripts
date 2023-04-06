# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:56:26 2020

@author: jeff
"""
import subprocess
import shutil
import os

path = '/volumes/hd2/jeff/seq_data/export/'

found = set(os.listdir(path))

with open('files2find.txt') as file_in:
    for line in file_in:
        line = line.rstrip()

        temp = subprocess.Popen('locate ' + line,
                                shell = True,
                                executable = '/bin/bash',
                                stdout = subprocess.PIPE)
        
        locations = temp.communicate()[0]
        locations = locations.decode().split('\n')
        fastas = [i for i in locations if 'fastq' in i] 
                        
        for fasta in fastas:
            file_name = fasta.split('/')[-1]
            if file_name not in found:
                shutil.copyfile(fasta, path + file_name) 
                found.add(file_name)
                print(file_name)
                
        if len(fastas) == 0:
            print(line, 'not found')

