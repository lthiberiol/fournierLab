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

import sys

if __name__ != '__main__':
    sys.exit()

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

with open(sys.argv[1]) as configuration_file:

    #
    # try to find the block containing script's parameters
    parameter_block   = re.search( '^generate_reference_tree.py\s?\{([\s\S]*)\}', configuration_file.read(), re.M ) 
    if not parameter_block:
        sys.exit('Could not find parameters respective to "generate_reference_tree.py", sorry...')

    #
    # go through the block lines and atribute parameters to their respective variables
    for line in parameter_block.group(1).strip().split('\n'):
        if line.startswith('#') or not line.strip():
            continue

        key, value = line.split('=')
        key        = key.strip()
        value      = value.split('#')[0].strip() # ignore everything after a "#"

        if key == 'genomes_dataframe':
            genomes_dataframe = value
        elif key == 'folder_to_download_genomes':
            folder_to_download_genomes = value
        elif key == 'folder_to_save_formatted_faas':
            folder_to_download_genomes = value

    # check if all required variables were defined
    if not set( 'genomes_dataframe folder_to_download_genomes folder_to_save_formatted_faas'.split() ).issubset( locals() ):
        sys.exit( "Could not define values for the current parameters: %s" ', '.join( set( 'folder_to_save_formatted_faas'.split() ).difference( locals() ) ) )

genomes = pd.read_table( genomes_dataframe, index_col=0 )

def download_genomes( ftp_path, target_sufix='protein.faa.gz' ):
    try: 
        ftp_path = ftp_path.replace('ftp://', '')
        ncbi = ftp.FTP('ftp.ncbi.nlm.nih.gov')
        ncbi.login()
        ncbi.cwd( ftp_path.replace('ftp.ncbi.nlm.nih.gov/', '') )
        for assembly_file in ncbi.nlst():
            if assembly_file.endswith( target_sufix ):
                with open( folder_to_download_genomes+'/'+assembly_file, 'wb') as file_handle:
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

os.system( 'gunzip -d %s/*' %folder_to_download_genomes )
print 'Genomes downloaded!'

os.mkdir( folder_to_download_genomes )
for faa_file_name in os.listdir( folder_to_download_genomes+'/' ):
    assembly_acc = re.match( '^GCA_\d+\.\d+', faa_file_name ).group()

    faa                = SeqIO.parse( '%s/%s' %(folder_to_download_genomes, faa_file_name), 'fasta' )
    sequences_to_write = []
    for protein in faa:
        protein.id = '%s|%s' %(assembly_acc, protein.id)
        sequences_to_write.append( protein )

    SeqIO.write( sequences_to_write, '%s/%s.faa' %(folder_to_download_genomes, assembly_acc), 'fasta' )
print 'FAAs edited!'

##tree = ete3.Tree('RAxML_bestTree.ribosomal_concat')
##name_dict = { i:'%s %s' %(j, str(k)) for i,j,k in genomes['assembly_accession organism_name infraspecific_name'.split()].values}
##for leaf in tree.get_leaves():
##    leaf.name = name_dict[ leaf.name ]
##tree.write( outfile='RAxML_bestTree.ribosomal_concat-organism_names' )
