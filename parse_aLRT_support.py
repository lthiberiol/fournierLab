#!/usr/bin/env python2.7
#coding: utf-8

############################################################
#                                                          #
# Script to change RAxML aLRT support trees into FigTree   #
#   happy style (move support values to its regular place) #
#                                                          #
#                                       L. Thibério Rangel #
#                                     lthiberiol@gmail.com #
#                                                          #
############################################################

from Bio import Phylo
import sys
import re

input_tree = sys.argv[-1]

tree = Phylo.read(input_tree, 'newick', comments_are_confidence=True)
for node in tree.get_nonterminals():
    if not node.comment:
        continue
    node.confidence = float(node.comment)
    node.comment = None

newick = tree.format('newick')
print re.sub('\):(\d+\.\d+:\d+\.\d+)', ')\\1', newick).strip()
