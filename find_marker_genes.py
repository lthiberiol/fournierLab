#!/usr/bin/env python
#encoding: utf-8

############################################################
#                                                          #
# Script to assess taxonomy available in Silva and define  #
#     genomes to sampled...                                #
#                                                          #
#                                       L. Thibério Rangel #
#                                     lthiberiol@gmail.com #
#                                                          #
############################################################

from sys import exit

if __name__ != '__main__':
    exit()

import ete3
import re
import time
import xml.etree.ElementTree as et
import pandas as pd
import os
import pickle as pkl
from copy import deepcopy
import numpy as np
import multiprocessing
import random
import itertools
import ftplib as ftp
from Bio import SeqIO

genomes = pd.read_table( 'selected_genomes.tab', index_col=0 )

def download_genomes( ftp_path, target_sufix='protein.faa.gz' ):
    try: 
        ftp_path = ftp_path.replace('ftp://', '')
        ncbi = ftp.FTP('ftp.ncbi.nlm.nih.gov')
        ncbi.login()
        ncbi.cwd( ftp_path.replace('ftp.ncbi.nlm.nih.gov/', '') )
        for assembly_file in ncbi.nlst():
            if assembly_file.endswith( target_sufix ):
                with open( 'selected_genomes/'+assembly_file, 'wb') as file_handle:
                    ncbi.retrbinary( "RETR %s" %assembly_file, file_handle.write )
                    ncbi.quit()
                    return 0
        return 1
    except:
        return 1

pool              = multiprocessing.Pool( processes=20 )
retrieved_genomes = pool.map(download_genomes, genomes.ftp_path.values )
pool.close()
pool.join()

os.system( 'gunzip -d selected_genomes/*' )
print 'Genomes downloaded!'

os.mkdir( 'formated_faas' )
for faa_file_name in os.listdir( 'selected_genomes/' ):
    assembly_acc = re.match( '^GCA_\d+\.\d+', faa_file_name ).group()

    faa                = SeqIO.parse( 'selected_genomes/%s' % faa_file_name, 'fasta' )
    sequences_to_write = []
    for protein in faa:
        protein.id = '%s|%s' %(assembly_acc, protein.id)
        sequences_to_write.append( protein )

    SeqIO.write( sequences_to_write, 'formated_faas/%s.faa' %assembly_acc, 'fasta' )
print 'FAAs edited!'



def hmmsearch( (hmm, faa), output_folder='/mnt/work2/hgt/greg/fournierLab/ribo_searches', threads=2 ):
   subprocess.call( ['hmmsearch',
                     '-o', '%s/%s-search_for-%s' %(output_folder, faa.replace('.faa', ''), hmm.replace('.hmm', '')),
                     '--acc',
                     '--noali',
                     '--cpu', '%i' %threads,
                     '/mnt/work2/hgt/greg/ribo_db/%s' %hmm,
                     '/mnt/work2/hgt/greg/fournierLab/formated_faas/%s' %faa] )
