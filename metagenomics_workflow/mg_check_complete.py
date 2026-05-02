#!/usr/bin/env python3
"""
Script to check which samples from samples.txt need to be processed.
Creates samples_resume.txt with only samples that don't have corresponding directories.
"""

import os
import sys

input_file = sys.argv[1]

if input_file.endswith('_resume.txt'):
    base_name = input_file[:-11]
elif input_file.endswith('.txt'):
    base_name = input_file[:-4]
else:
    base_name = input_file
    
output_file = f"{base_name}_resume.txt"

to_run = []

with open(input_file, 'r') as file_in:
    for line in file_in.readlines():
        line = line.rstrip()
        line = line.split()[-1] # if more than one column, assume you want last column
        dir_name = 'binning_' + line
        if not os.path.exists(dir_name):
            to_run.append(line)
            print(line, 'has been added')
        else:
            print(line, 'is already complete')
            
with open(output_file, 'w') as file_out:
    file_out.write('\n'.join(to_run))