# -*- coding: utf-8 -*-
"""
Created on Tue May 21 19:54:18 2019

@author: jeff
"""
from Bio import Seq


with open('anl_191021/191021_Mapping_File_Bowman_18S.oligos', 'w') as file_out:
    with open('anl_191021/191217_Bowman_18S_AS_191219.txt', 'r') as file_in:
        gene = '18S'
        for line in file_in:
            if line.startswith('#') == False:
               line = line.rstrip()
               line = line.split('\t')
               if len(line) > 1:
                   print(line[0] + '_' + gene, line[1])
                   print(line[0] + '_' + gene + '\t' + line[1], file = file_out)

#with open('191007_Mapping_File_Bowman_18S.oligos', 'w') as file_out:                   
#    with open('anl_191007/191107_Bowman_18S_AS_191104.txt', 'r') as file_in:
#        gene = '18S'
#        for line in file_in:
#            if line.startswith('#') == False:
#               line = line.rstrip()
#               line = line.split('\t')               
#               if len(line) > 1:
#                   barcode = line[1]
#                   barcode = Seq.reverse_complement(barcode)
#                   print(line[0] + '_' + gene, barcode)
#                   print(line[0] + '_' + gene + '\t' + barcode, file = file_out)
           