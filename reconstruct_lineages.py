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

def follow_lineage( genome_series ):
    lineage = []
    parent_id = genome_series['parent tax_id']
    while True:
        parent = taxonomy.loc[parent_id].squeeze()
        lineage.append( (parent_id, parent['rank']) )

        parent_id = parent['parent tax_id']
        if parent_id == '1':
            break
    return lineage

print "\t**Loadind main raw DataFrames ..."
header_names=['tax_id', 'parent tax_id', 'rank', 'embl code', 'division id', 'inherited div flag', 'genetic code id', 'inherited GC', 'mitochondrial genetic', 'inherited MGC flag', 'GenBank hidden flag', 'hidden subtree root', 'comments']
taxonomy                  = pd.read_table( '../ncbi_taxonomy/nodes.dmp', sep='\t\|\t', header=None, names=header_names, index_col=False, engine='python' )
taxonomy.tax_id           = taxonomy.tax_id.astype(str)
taxonomy['parent tax_id'] = taxonomy['parent tax_id'].astype(str)
taxonomy                  = taxonomy.set_index( 'tax_id' )

genbank_summary = pd.read_table( '/mnt/work2/hgt/greg/assembly_summary_genbank.txt', dtype={'taxid':str, 'infraspecific_name':str} )
print "\t**Main raw DataFrames loaded!!"

print '\t**Assembling lineages'
genomes_to_reconstruct = taxonomy[ taxonomy.index.isin( genbank_summary.taxid) ]
multithread_input = [ genome_series for index, genome_series in genomes_to_reconstruct.iterrows() ]

print '\t\t**starting threads!\n'
pool = multiprocessing.Pool( processes=30 )
lineages = pool.map(follow_lineage, multithread_input )
pool.close()
pool.join()

print '\t**Parsing lineages'
genome_lineages = pd.DataFrame()

#obligatory_taxonomic_levels = set(['phylum', 'family', 'order', 'class', 'genus'])
for taxid, lineage in zip(genomes_to_reconstruct.index, lineages):
    lineage.reverse()
    lineage = np.asarray( lineage )

#    if not obligatory_taxonomic_levels.issubset( lineage[:,1] ):
#        continue

    tmp_lineage = []
    for rank_id, rank_name in lineage:
        if rank_name == 'no rank':
            continue
        tmp_lineage.append( [rank_name, [rank_id]] )

    if not tmp_lineage:
        continue
    elif len(tmp_lineage) == 1:
        tmp_lineage = pd.DataFrame.from_items(tmp_lineage).iloc[0]
    else:
        tmp_lineage = pd.DataFrame.from_items(tmp_lineage).squeeze()
    tmp_lineage.name = taxid
    genome_lineages = genome_lineages.append( tmp_lineage )

print '\t**Saving lineages'
genome_lineages.to_csv( 'lineages_from_genbank_summary2.tab', sep='\t' )

print '\t**Parsing humand readable names...'
header_names      = ['tax_id', 'name_txt', 'unique name', 'name class']
taxa_names        = pd.read_table( '../ncbi_taxonomy/names.dmp', sep='\t\|\t', header=None, names=header_names, index_col=False, engine='python' )
taxa_names.tax_id = taxa_names.tax_id.astype(str)
taxa_names        = taxa_names.set_index( 'tax_id' )
taxa_names        = taxa_names[taxa_names['name class'] == 'scientific name\t|']

print '\t**Adding human readable names...'
for column in genome_lineages.columns:
    genome_lineages.loc[genome_lineages[column].notnull(), column] = taxa_names.loc[ genome_lineages.loc[genome_lineages[column].notnull(), column].astype(int).astype(str), 'name_txt' ].values

print '\t**Saving lineages, again!'
genome_lineages.to_csv( 'lineages_from_genbank_summary2-txt_names.tab', sep='\t' )
