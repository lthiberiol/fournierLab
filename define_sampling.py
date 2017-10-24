#!/usr/bin/env python
#enconding: utf-8

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
from Bio import Entrez, SeqIO
Entrez.email = 'lthiberiol@gmail.com'

target_taxonomies = 'alphaproteobacteria cyanobacteria'.split()

#
# parse information from xml to DataFrame, if it wasn't before
if not os.path.isfile( '../silva/silva.tab' ):
    print "\t*Loading Silva's full XML... (it might take a while, go grab something to eat)"
    root     = et.parse('../silva/silva-sub.xml').getroot()
    data     = {}
    xml_size = len(root)

    count      = 0
    percentage = 0
    for child in root:
        #
        # small useless thing to print progress
        count += 1
        if (count*100) / xml_size > percentage:
            percentage = (count*100) / xml_size
            print "%i%% complete         \r" %percentage,

        if child.tag != 'species':
            continue

        strain      = child.find( 'strain'        ).text
        taxid       = child.find( 'tax_xref_embl' ).text
        full_name   = child.find( 'full_name'     ).text
        tax_name    = child.find( 'tax_embl_name' ).text
        description = child.find( 'description'   ).text
        tax_embl    = child.find( 'tax_embl'      ).text
        if re.search( 'uncultured', full_name, re.I ) or re.search( 'uncultured', description, re.I ):
            print '** Fucked, still contains uncultivated taxa...'
            continue

        tax_class, tax_family = re.search( '.*;(\S+?ales)?;(\S+?aceae)?', tax_embl ).groups()

        data[ child.attrib['name'] ] = { 'full_name':full_name,
                                         'description':description,
                                         'tax_class':tax_class,
                                         'tax_family':tax_family,
                                         'infraspecific_name':'strain=%s' %strain,
                                         'taxid':taxid,
                                         'organism_name':tax_name,
                                         'taxonomy':tax_embl,
                                         'accession':child.find('acc').text,
                                       }

    leaf_info = pd.DataFrame.from_dict( data, orient='index' )
    leaf_info.to_csv( '../silva/silva-test.tab', sep='\t' )
    print '\t**Finished parsing XML...'
else:
    print "\t*Loading Silva's pre-parsed XML"
    leaf_info = pd.read_table( '../silva/silva-test.tab', dtype={'insdc':str, 'taxid':str, 'strain':str}, index_col=0 )
    print '\t**Finished loading XML...'

leaf_info_sub = leaf_info[ (leaf_info.taxonomy.str.contains('Cyanobacteria')) |
                           (leaf_info.taxonomy.str.contains('Alphaproteobacteria')) ].copy()

genbank_summary = pd.read_table( '/mnt/work2/hgt/greg/assembly_summary_genbank.txt', dtype={'taxid':str, 'infraspecific_name':str} )

samples = leaf_info_sub.merge( genbank_summary, on=['taxid', 'infraspecific_name', 'organism_name'], right_index=True )
samples['silva_name'] = samples.index
samples.drop_duplicates( subset=('silva_name', 'assembly_level', 'genome_rep', 'refseq_category'), inplace=True )
duplications = samples.index[samples.index.duplicated()]
for duplicate in duplications:
    print duplicate

    if samples.loc[duplicate, 'refseq_category'].unique().shape[0] > 1:
        for refseq_category in ['reference genome', 'representative genome', 'na']:
            if refseq_category in samples.loc[duplicate, 'refseq_category'].unique():
                best_category = refseq_category
                break
        to_keep = samples.loc[duplicate][samples.loc[duplicate, 'refseq_category'] == best_category].squeeze()
        samples.drop( labels=duplicate, inplace=True )
        samples.loc[duplicate] = to_keep.copy()

    elif samples.loc[duplicate, 'assembly_level'].unique().shape[0] > 1:
        for assembly_level in ['Complete genome', 'Chromosome', 'Scaffold', 'Contig']:
            if assembly_level in samples.loc[duplicate, 'assembly_level'].unique():
                best_assembly = assembly_level
                break
        to_keep = samples.loc[duplicate][samples.loc[duplicate, 'assembly_level'] == best_assembly].squeeze()
        samples.drop( labels=duplicate, inplace=True )
        samples.loc[duplicate] = to_keep.copy()

    elif samples.loc[duplicate, 'genome_rep'].unique().shape[0] > 1:
        for genome_rep in ['Full', 'Partial']:
            if genome_rep in samples.loc[duplicate, 'genome_rep'].unique():
                best_representation = genome_rep
                break
        to_keep = samples.loc[duplicate][samples.loc[duplicate, 'genome_rep'] == best_representation].squeeze()
        samples.drop( labels=duplicate, inplace=True )
        samples.loc[duplicate] = to_keep.copy()
print '\t**Merge done!'

#
# edit and load Silva's phylogeny
print "\n\t**Loading Silva's tree..."
reference_branch_labels      = {}
reference_raw_tree           = ''
branch_label_count           = 1
for line in open('../silva/silva-cultivated_only.tree').xreadlines():
    tmp_branch_name = re.search( "\)('.+?'):", line )
    if tmp_branch_name:
        reference_branch_labels[ "'branchLabel%i'" %branch_label_count ] = tmp_branch_name.group(1)
        line = line.replace( tmp_branch_name.group(1), "'branchLabel%i'" %branch_label_count )
        branch_label_count += 1
    reference_raw_tree += line
silva_tree = ete3.Tree( reference_raw_tree, format=1 )

print "**Pruning shitty taxa!!"
loop_count = 0
for node in silva_tree.traverse( is_leaf_fn=lambda node: True if target_leaves.issuperset(node.get_leaf_names()) or target_leaves.isdisjoint(node.get_leaf_names()) else False ):
    loop_count += 1

    if target_leaves.isdisjoint( node.get_leaf_names() ):
        node.detach()
print "**Pruning is done!!"
print '\t**Iterations required: %i' %loop_count

genus = {}
for window in range(0, samples.species_taxid.unique().shape[0], 200):
    query = Entrez.efetch( id=samples.species_taxid.unique()[window:window+200].astype(str).tolist(), db='taxonomy')
    result = et.parse(query).getroot()
    for entry in result:
        lineage = entry.find('LineageEx')
        for rank in lineage.getchildren():
            if rank.find( 'Rank' ).text == 'genus':
                genus[entry.find( 'TaxId' ).text] = rank.find( 'ScientificName' ).text
#                entry_genus = rank.find( 'ScientificName' ).text
                break
    time.sleep(3)
print '** Done retrieving genus lineage data!'

for leaf in silva_tree.get_leaves():
    if not samples.loc[leaf.name, 'species_taxid'] or pd.isnull( samples.loc[leaf.name, 'species_taxid'] ):
        print leaf.name
        break
    else:
        try:
            leaf.add_feature( 'genus', genus[samples.loc[leaf.name, 'species_taxid'].astype(str)] )
        except KeyError:
            leaf.add_feature( 'genus', 'NA' )

monophyletic_genera = []
paraphyletic_genera = []
for genus in set( genus.values() ):
    isMonophyletic, clade_type, fucking_up = silva_tree.check_monophyly( values=[genus], target_attr='genus' )
    if isMonophyletic:
        monophyletic_genera.append( genus )
    else:
        paraphyletic_genera.append( genus )
            

#
# traverse tree and select representatives
representatives = {}
for family in leaf_info.tax_family.unique():
    taxa = leaf_info.index[ leaf_info.tax_family == family ]
    isMonophyletic, cladeType, fuckingThingsUp = silva_tree.check_monophyly( taxa, target_attr='name' )
    if not isMonophyletic:
        continue
    
    representatives[family] = []
    for child in silva_tree.get_common_ancestor( taxa ):
        representatives[family].append( child.get_topology_only() )
