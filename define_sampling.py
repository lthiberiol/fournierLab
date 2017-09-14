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
import xml.etree.ElementTree as et
import pandas as pd
import os

target_taxonomies = 'alphaproteobacteria cyanobacteria'.split()

#
# parse information from xml to DataFrame, if it wasn't before
if not os.path.isfile( '../silva/silva.tab' ):
    data = {}
#    root = et.parse('../silva/silva.xml').getroot()
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

        full_name   = child.find( 'full_name'   ).text
        description = child.find( 'description' ).text
        if re.search( 'uncultured', full_name, re.I ) or re.search( 'uncultured', description, re.I ):
            continue

        tax_silva = child.find( 'tax_slv' ).text
        tax_flag = False
        for target in target_taxonomies:
            if re.search( target, tax_silva, re.I ):
                tax_flag = True
                break
        if not tax_flag:
            continue

        try:
            tax_class, tax_family = re.search( '.*;(\S+?ales);(\S+?aceae)', tax_silva ).groups()
        except AttributeError:
            continue

        data[ child.attrib['name'] ] = { 'full_name':full_name,
                                         'description':description,
                                         'tax_class':tax_class,
                                         'tax_family':tax_family,
                                         'accession':child.find('acc').text,
                                       }

    leaf_info = pd.DataFrame.from_dict( data, orient='index' )
    leaf_info.to_csv( '../silva/silva.tab', sep='\t' )
    print '\t**Finished parsing XML...'
else:
    leaf_info = pd.read_table( '../silva/silva.tab', index_col=1 )
    print '\t**Finished loading XML...'

#
# edit and load Silva's phylogeny
branch_labels      = {}
branch_label_count = 1
raw_tree           = ''
for line in open('../silva/SSURefNR99_1200_slv_128.tree').xreadlines():
    tmp_branch_name = re.search( "\)('.+?'):", line )
    if tmp_branch_name:
        branch_labels[ "'branchLabel%i'" %count ] = tmp_branch_name.group(1)
        line = line.replace( tmp_branch_name.group(1), "'branchLabel%i'" %count )
        count += 1
    raw_tree += line

silva_tree = ete3.Tree( raw_tree, format=1 )
#silva_tree.prune( leaf_info.index, preserve_branch_length=True )
silva_tree.prune( leaf_info.index )
print '\t**Finished loading tree...'
