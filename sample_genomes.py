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

def select_taxa_with_best_data( tmp_df, taxon_rank ):
    df = tmp_df.copy()

    tmp_grouped = df.groupby( by=taxon_rank )
    best_taxa = []
    for column_data in itertools.product( ['reference genome', 'representative genome', 'na'], ['complete genome', 'chromosome', 'scaffold', 'contig'], ['full', 'partial'] ):
        for taxon, taxon_df in tmp_grouped:
            if column_data in taxon_df['refseq_category assembly_level genome_rep'.split()].values:
                best_taxa.append( taxon_df )

            if len(best_taxa) == 2:
                break
        if len(best_taxa) == 2:
            break

    return best_taxa

def select_representative( tmp_df ):
    df = tmp_df.copy()
    df = tmp_df['refseq_category assembly_level genome_rep'.split()].copy()

    for column_data in itertools.product( ['reference genome', 'representative genome', 'na'], ['complete genome', 'chromosome', 'scaffold', 'contig'], ['full', 'partial'] ):
        column_data = pd.Series(index='refseq_category assembly_level genome_rep'.split(), data=column_data)
        df[ df == column_data ]
        
        if df.shape[0]:
            break

    df = tmp_df.loc[ df.index ].copy()
    if df.shape[0] == 1:
        return df
    else:
        return pd.DataFrame(df.iloc[0]).T

###    if df['refseq_category'].unique().shape[0] > 1:
###        for refseq_category in ['reference genome', 'representative genome', 'na']:
###            if refseq_category in df['refseq_category'].unique():
###                best_category = refseq_category
###                break
###        df.drop( df[df['refseq_category'] != best_category].index, axis='index', inplace=True )
###
###    if df['assembly_level'].unique().shape[0] > 1:
###        for assembly_level in ['complete genome', 'chromosome', 'scaffold', 'contig']:
###            if assembly_level in df['assembly_level'].unique():
###                best_assembly = assembly_level
###                break
###        df.drop( df[df['assembly_level'] != best_assembly].index, axis='index', inplace=True )
###
###    if df['genome_rep'].unique().shape[0] > 1:
###        for genome_rep in ['full', 'partial']:
###            if genome_rep in df['genome_rep'].unique():
###                best_representation = genome_rep
###                break
###        df.drop( df[df['genome_rep'] != best_representation].index, axis='index', inplace=True )
###
###    if df.shape[0] == 1:
###        return df
###    else:
###        return pd.DataFrame(df.iloc[0]).T

print "\t**Loadind main raw DataFrames ..."
genome_lineages = pd.read_table( 'lineages_from_genbank_summary2.tab', sep='\t', index_col=0, dtype=str )
genbank_summary = pd.read_table( '/mnt/work2/hgt/greg/assembly_summary_genbank.txt', dtype={'taxid':str, 'infraspecific_name':str} )

genome_lineages.index = genome_lineages.index.astype( str )
genome_lineages['taxid'] = genome_lineages.index
genbank_summary['refseq_category']  = genbank_summary['refseq_category'].str.lower()
genbank_summary['assembly_level']   = genbank_summary['assembly_level'].str.lower()
genbank_summary['genome_rep']       = genbank_summary['genome_rep'].str.lower()
print "\t**Main raw DataFrames loaded!!"

print "**\tFiltering down to genomes from target taxonomies.."
# Alphaproteobacteria = 28211
# Cyanobacteria       = 1117
target_taxonomies = ['28211', '1117']

genome_lineages  = genome_lineages[ (genome_lineages.phylum.isin( target_taxonomies )) | 
                                   (genome_lineages['class'].isin( target_taxonomies )) ]

merged           = genome_lineages.merge( genbank_summary, how='inner', on='taxid', suffixes=['_lineage','_summary'] )

merged.drop( merged.genus[merged.genus.isnull()].index, axis='index', inplace=True )
grouped_by_genus = merged.groupby( by='genus' )
selected_genomes = pd.DataFrame( columns=merged.columns )
single_representative = []
two_representative    = []
multiple_species      = []
multiple_strains      = []
for genus, df in grouped_by_genus:
    if df.shape[0] == 1:
        selected_genomes = selected_genomes.append( df )
        single_representative.append( genus )
    else:
        if df['species'].unique().shape[0] == 1 and df['infraspecific_name'].unique().shape[0] == 1:
            selected_genomes = selected_genomes.append( select_representative( df ) )
            two_representative.append( genus )
        elif df.species.unique().shape[0] > 1:
            for taxon in select_taxa_with_best_data( df, 'species' ):
                selected_genomes = selected_genomes.append( select_representative( taxon ) )
            multiple_species.append( genus )
        else:
            for taxon in select_taxa_with_best_data( df, 'infraspecific_name' ):
                selected_genomes = selected_genomes.append( select_representative( taxon ) )
            multiple_strains.append( genus )
print "**\tFiltering done!!!"
