#!/usr/bin/env python
#====================

import ete3
import os
import numpy as np

def remove_long_branches(tree_file):
    tree = ete3.Tree('/home/thiberio/Alphas_and_Cyanos/ranger_input_trees/%s' %tree_file)
    branch_lengths = []
    for node in tree.traverse():
        branch_lengths.append(node.dist)

    quartile_3 = np.percentile(branch_lengths, 75)

    taxa_to_remove = set()
    for node in tree.traverse():
        if node.dist > quartile_3*100:
            taxa_to_remove.update(node.get_leaf_names())
    tree.prune(taxa_to_remove.symmetric_difference(tree.get_leaf_names()), preserve_branch_length=True)

    tree.write(outfile='/home/thiberio/Alphas_and_Cyanos/ranger_input_trees-no_long_branches/%s' %tree_file, dist_formatter='%.20f')

    return 1
