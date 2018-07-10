#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
import multiprocessing
import pickle as pkl
import linecache

os.chdir('/work/Alphas_and_Cyanos')

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

def rename_branches(reconciliation_file, tree):
    branches = re.findall('^(m\d+) = LCA\[(\S+), (\S+)\]:', reconciliation_file, re.M)
    for name, leaf1, leaf2 in branches:
        node = tree.get_common_ancestor(leaf1, leaf2)
        if node.name:
            print node.get_leaves_by_name()
            break
        node.name = name
    return tree

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
                tree_to_root.set_outgroup(fucking_up.pop())
                equivalent = tree_to_root.get_common_ancestor(node.get_leaf_names())
                tree_to_root.set_outgroup(equivalent)
            break

    return tree_to_root

def parse_aggregated(folder, threshold=0.9, leaves_allowed=False):
    with cd(folder):
        aggregated    = open('aggregated').read()
        gene_tree     = {'named':ete3.Tree(linecache.getline('%s-MAD.ranger_out1' %folder, 8), format=1)}

    gene_tree['support'] = root_like(gene_tree['named'], ete3.Tree('/work/Alphas_and_Cyanos/ranger_input_trees-no_long_branches/%s.tree' %folder))
    gene_tree            = rename_branches(aggregated, gene_tree['support'])

    num_replicates = float(re.match('Processed (\d+) files', aggregated).group(1))

    if not leaves_allowed:
        transfers = re.findall('^(m\d+) = .*, Transfers = [^0]\d+?\], \[Most Frequent mapping --> (n\d+), (\d+) times\], \[Most Frequent recipient --> (n\d+), (\d+) times\].', aggregated, re.M)
    else:
        transfers = re.findall('^(m\d+) = .*, Transfers = [^0]\d+?\], \[Most Frequent mapping --> (\S+), (\d+) times\], \[Most Frequent recipient --> (\S+), (\d+) times\].',   aggregated, re.M)

    supported_transfers = []
    for donor_map, donor, ranger_confidence_donor, recipient, ranger_confidence_recipient in transfers:
        if int(ranger_confidence_donor) < threshold*num_replicates or int(ranger_confidence_recipient) < threshold*num_replicates:
            continue
        supported_transfers.append((donor_map, donor, recipient))

    selected_transfers = []
    for donor_map_name, donor_name, recipient_name in supported_transfers:
        donor_map = gene_tree.search_nodes(name=donor_map_name)[0]
        if donor_map.support < 95:
            continue

        recipient_map_search = re.search('^({children[0]}|{children[1]}).*Most Frequent mapping --> {recipient}'.format(recipient=recipient_name, children=[child.name for child in donor_map.children]), aggregated, re.M)
        if recipient_map_search:
            recipient_map_name = recipient_map_search.group(1)
            selected_transfers.append({'donor':donor_name, 'recipient':recipient_name, 'donor_map':donor_map_name, 'recipient_map':recipient_map_name})

    return {folder:selected_transfers}

with cd('reconciliations/mad_roots'):
    pool = multiprocessing.Pool(processes=4)
    transfers = pool.map(parse_aggregated, os.listdir('.'))
    pool.close()
    pool.join()

    yeah = {}
    for filtered in transfers:
        if filtered.values() != [[]]:
            yeah.update(filtered)

out = open('aggregated/mad-90_threshold.pkl', 'wb')
pkl.dump(yeah, out)
out.close()