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

input_tree = sys.argv[-1]

tree = Phylo.read(input_tree, 'newick')
for node in tree.get_nonterminals():
    if not node.comment:
        continue
    node.confidence = float(node.comment)
    node.comment = None
Phylo.write( tree, '%s.parsed' %input_tree, 'newick', format_branch_length='%1.20f')
