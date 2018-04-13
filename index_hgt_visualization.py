#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import itertools
import re
from commands import getoutput
import seaborn as sns
from matplotlib import pyplot as plt
import multiprocessing

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

def name_matching_branches( named, unamed ):
    if named.get_topology_id() == unamed.get_topology_id():
        unamed.name = named.name

    for node1 in unamed.children:
        if node1.is_leaf():
            continue

        for node2 in named.children:
            if node2.is_leaf():
                continue

            if node1.get_topology_id() == node2.get_topology_id():
                node1.name = node2.name
                name_matching_branches( node2, node1 )

def root_like( ref_tree, tree_to_root ):
    if len(ref_tree) < len(tree_to_root):
        tree_to_root.prune(ref_tree.get_leaf_names())
    elif len(ref_tree) > len(tree_to_root):
        print 'Different tree sizes'
        return
    outgroup = ''
    for node in sorted( ref_tree.children, key=len ):
        if node.is_leaf():
            putative_outgroups = tree_to_root.get_leaves_by_name(node.name)
            for leaf in putative_outgroups:
                tree_to_root.set_outgroup( leaf )
                if tree_to_root.get_topology_id() == ref_tree.get_topology_id():
                    outgroup = leaf
                    break
        else:
            outgroup_members = []
            for leaf in node:
                outgroup_members.append( tree_to_root.get_leaves_by_name(leaf.name) )

            for outgroup_combination in itertools.product( *outgroup_members ):
                if tree_to_root.check_monophyly( [x.name for x in outgroup_combination], 'name' )[0]:
                    putative_outgroup = tree_to_root.get_common_ancestor( outgroup_combination )
                    tree_to_root.set_outgroup( putative_outgroup )
                    if tree_to_root.get_topology_id() == ref_tree.get_topology_id():
                        outgroup = putative_outgroup
                        break
        if outgroup:
            return tree_to_root

    if not outgroup:
        print 'fuck'
        return None

os.chdir('/work/Alphas_and_Cyanos/groups_with_nested_transfers')

homologue_group = '004119'
named_gene_tree = ete3.Tree('%s.optimal_root.tree'        %homologue_group, format=1)
for node in named_gene_tree.traverse():
    if node.name == 'm1' or node.is_leaf() or not node.name:
        continue

    node_name, node_support = re.search('^(m\d+?)(100|\d{2})?$', node.name).groups()
    node.support            = int(node_support if node_support else 0)
    node.name               =     node_name

gene_tree       = ete3.Tree('RAxML_result.%s' %homologue_group)

for named_node in named_gene_tree.traverse():
    if named_node.is_leaf():
        continue

    equivalent_node = gene_tree.get_common_ancestor(named_node.get_leaf_names())
    if len(named_node) != len(equivalent_node):
        print named_node.name
        break




name_matching_branches(named_gene_tree, gene_tree)
named_gene_tree.get_topology_id()
gene_tree.get_topology_id()