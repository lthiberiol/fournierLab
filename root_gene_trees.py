#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
import operator
import multiprocessing
import pickle as pkl
import pandas as pd
import numpy as np
import random
from commands import getoutput
import itertools
import seaborn as sns
from matplotlib import pyplot as plt
from time import time
import random
import plotly
import plotly.plotly as ptl
from plotly import graph_objs as go
ptl.sign_in('lthiberiol', 'm15ikp59lt')

class cd:
    """
    Context manager for changing the current working directory
    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

def jaccard(array1, array2):
    return 1 - len(array1.intersection(array2))/float(len(array1.union(array2)))

def root_like( ref_tree, tree2 ):
    tree_to_root = tree2.copy()
    for node in sorted( ref_tree.children, key=len ):
        if node.is_leaf():
            leaf = tree_to_root.get_leaves_by_name(node.name)[0]
            tree_to_root.set_outgroup(leaf)
            break
        else:
            is_it_monophyletic, clade_type, fucking_up = tree_to_root.check_monophyly(node.get_leaf_names(), 'name', unrooted=False)
            if is_it_monophyletic:
                equivalent = tree_to_root.get_common_ancestor(node.get_leaf_names())
                tree_to_root.set_outgroup(equivalent)
            else:
                gene_tree.set_outgroup(fucking_up.pop())
                equivalent = tree_to_root.get_common_ancestor(node.get_leaf_names())
                tree_to_root.set_outgroup(equivalent)
            break

    return tree_to_root

ncbi = ete3.NCBITaxa()

os.chdir('/work/Alphas_and_Cyanos')
species_tree = ete3.Tree('rooted_partitions-with_BB_support.treefile', format=0)
gene_tree    = ete3.Tree('ranger_input_trees/000070.tree',             format=0)

genomes = set()
for leaf in gene_tree.get_leaves():
    genome = re.match('GC[AF]\d+', leaf.name).group()
    leaf.add_feature('genome', genome)
    genomes.add(genome)

pruned_species_tree = species_tree.copy()
pruned_species_tree.prune(genomes)

child1           = set(pruned_species_tree.children[0].get_leaf_names())
child2           = set(pruned_species_tree.children[1].get_leaf_names())
gene_tree_leaves = set(gene_tree.get_leaves())

loop_count        = 0
possible_rootings = {}
for node in gene_tree.traverse():
    if node.is_leaf():
        continue

    bipartition1 = set([leaf.genome for leaf in node.get_leaves()])
    bipartition2 = genomes.difference(bipartition1)

    if not bipartition2:
        continue

    loop_count += 1
    node.add_feature('name', 'node_%i' %loop_count)

    jaccard_distance1 = jaccard(child1, bipartition1) + jaccard(child2, bipartition2)
    jaccard_distance2 = jaccard(child1, bipartition2) + jaccard(child2, bipartition1)
    if jaccard_distance1 < jaccard_distance2:
        possible_rootings[node.name] = jaccard_distance1
    else:
        possible_rootings[node.name] = jaccard_distance2

sorted_by_jaccard = sorted(possible_rootings.items(), key=operator.itemgetter(1))
optimal_root      = gene_tree.search_nodes(name=sorted_by_jaccard[0][0])[0]
gene_tree.set_outgroup(optimal_root)
gene_tree.write(outfile='test_optimal_rootings/000070.tax_optimal_root.tree', format=0)

########################################################################################################################
# Compare rooting with RANGER's                                                                                        #
########################################################################################################################

gene_tree_optRooted = root_like(ete3.Tree('test_optimal_rootings/000070.ranger_rooting.tree'), gene_tree)

random_node_sample        = random.sample(list(gene_tree_optRooted.traverse()), 100)
optimal_root_in_optRooted = gene_tree_optRooted.search_nodes(name=sorted_by_jaccard[0][0])[0]
ancestors                 = optimal_root_in_optRooted.get_ancestors()
ancestors.insert(0, optimal_root_in_optRooted)

with cd('test_optimal_rootings/random_rootings'):
    count = 0
    for node in random_node_sample:
        gene_tree_optRooted.set_outgroup(node)
        out    = open('random_gene_tree_rootings-%i.ranger_input' %count, 'wb')
        out.write(species_tree.write(format=9))
        out.write('\n')
        out.write(gene_tree_optRooted.write(format=9))
        out.close()
        count += 1

with cd('test_optimal_rootings'):
    count = 0
    for ancestor in ancestors:
        gene_tree_optRooted.set_outgroup(ancestor)
        print ancestor.name
        print gene_tree_optRooted.get_topology_id()
        out    = open('steps_to_optRoot-%i.ranger_input' %count, 'wb')
        out.write(species_tree.write(format=9))
        out.write('\n')
        out.write(gene_tree_optRooted.write(format=9))
        out.close()
        count += 1

reconciliation_costs = [int(cost) for cost in re.findall(
    'The minimum reconciliation cost is: (\d+)',
    getoutput('grep "The minimum reconciliation cost is" test_optimal_rootings/steps_to_optRoot-{0..13}.ranger_input.out'
))]
























