#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
from commands import getoutput
import multiprocessing
import itertools
import pickle as pkl
import seaborn as sns
from matplotlib import pyplot as plt

os.chdir('/work/Alphas_and_Cyanos')
named_reference_tree = ete3.Tree('rooted_partitions-with_named_branches.treefile', format=1)

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

def assess_compatibilities(donor, recipient, reticulation):
    #
    # ignore reticulation if there is a duplication in it
    if len(set([leaf.genome_name for leaf in reticulation.get_leaves()])) < len(reticulation):
        return False

    #
    # assess donor branch compatibility
    rf, rf_max, names, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = reticulation.robinson_foulds(donor, attr_t1='genome_name')
    if rf_max:
        rf_normalized = rf/float(rf_max)
    else:
        rf_normalized = 0.0

    if rf_normalized > 0.2:
        return False

    #
    # assess recipient branch compatibility
    rf, rf_max, names, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = reticulation.robinson_foulds(recipient, attr_t1='genome_name')
    if rf_max:
        rf_normalized = rf/float(rf_max)
    else:
        rf_normalized = 0.0

    if rf_normalized > 0.2:
        return False

    return True

def assess_transfers(transfer_file, reference_tree=named_reference_tree, debug=False):
    if not transfer_file.endswith('.pkl'):
        return None

    selected_transfers = []
    group              = transfer_file.replace('.pkl', '')
    gene_tree          = ete3.Tree('/work/Alphas_and_Cyanos/ranger_input_trees/%s.tree' %group)
    for leaf in gene_tree.get_leaves():
        leaf.name = re.sub('\.\d+', '', leaf.name.replace('GCA_', 'GCA'))
        leaf.add_feature('genome_name', leaf.name.split('_')[0])

    transfer_data = pkl.load(open(transfer_file))
    for donor, recipient, reticulation, topology_id in transfer_data:
        is_it_monophyletic, clade_type, fucking_up = gene_tree.check_monophyly(reticulation.get_leaf_names(), 'name', unrooted=False)
        if is_it_monophyletic:
            reticulation_with_support = gene_tree.get_common_ancestor(reticulation.get_leaf_names())
        else:
            gene_tree.set_outgroup(fucking_up.pop())

            if debug:
                is_it_monophyletic, clade_type, fucking_up = gene_tree.check_monophyly(reticulation.get_leaf_names(), 'name', unrooted=False)
                if not is_it_monophyletic:
                    raise KeyError('the reticulation is not monophyletic within the gene tree!')

            reticulation_with_support = gene_tree.get_common_ancestor(reticulation.get_leaf_names())

        if reticulation_with_support.support < 95:
            continue

        recipient_branch = reference_tree.search_nodes(name=recipient)[0].copy(method='deepcopy')
        donor_branch     = reference_tree.search_nodes(name=donor    )[0].copy(method='deepcopy')

##        flag = assess_compatibilities(donor_branch, recipient_branch, reticulation_with_support)
##        if not flag:
##            continue

        donor_leaves     = [leaf.name for leaf in reticulation_with_support.get_leaves() if leaf.genome_name in donor_branch.get_leaf_names()    ]
        recipient_leaves = [leaf.name for leaf in reticulation_with_support.get_leaves() if leaf.genome_name in recipient_branch.get_leaf_names()]

        if debug:
            is_it_monophyletic, clade_type, fucking_up = reticulation_with_support.check_monophyly(donor_leaves.get_leaf_names(), 'name', unrooted=False)
            if not is_it_monophyletic:
                raise KeyError('the donor leaves are not monophyletic within the reticulation!')

            is_it_monophyletic, clade_type, fucking_up = reticulation_with_support.check_monophyly(recipient_leaves.get_leaf_names(), 'name', unrooted=False)
            if not is_it_monophyletic:
                raise KeyError('the recipient leaves are not monophyletic within the reticulation!')

        donor_branch_in_reticulation     = reticulation_with_support.get_common_ancestor(    donor_leaves)
        recipient_branch_in_reticulation = reticulation_with_support.get_common_ancestor(recipient_leaves)

        if donor_branch_in_reticulation in recipient_branch_in_reticulation.get_ancestors():
            print '%s: recipient nested into donor' %group
            selected_transfers.append((donor, recipient, reticulation, topology_id))

    return {group:selected_transfers}


with cd('aggregated'):
    pool = multiprocessing.Pool(processes=6)
    results = pool.map(assess_transfers, os.listdir('.'))
    pool.close()
    pool.join()

usefull_results = {}
for entry in results:
    if entry and entry.values() != [[]]:
        for key, value in entry.items():
            if key not in usefull_results:
                usefull_results[key] = []
            usefull_results[key].extend(value)

### fig, axs = plt.subplots(nrows=2)
### axs[0].set_title('Donor branch Robinson-Foulds distances')
### sns.distplot(rf_values['donor'],     ax=axs[0])
### axs[1].set_title('Recipient branch Robinson-Foulds distances')
### sns.distplot(rf_values['recipient'], ax=axs[1])
### fig.tight_layout()
### fig.savefig('rf_distances.pdf', dpi=600)
###
### fig, axs = plt.subplots(nrows=2)
### axs[0].set_title('Donor branch normalized Robinson-Foulds distances')
### sns.distplot(rf_values['donor'],     ax=axs[0])
### axs[1].set_title('Recipient branch normalized Robinson-Foulds distances')
### sns.distplot(rf_values['recipient'], ax=axs[1])
### fig.tight_layout()
### fig.savefig('rf_distances_norm.pdf', dpi=600)
###
### transfer_supports     = []
### mapping_consistencies = []
### for element in results1:
###     if element == None:
###         continue
###     transfer_supports.extend(element[0])
###     mapping_consistencies.extend(element[1])
### fig, axs = plt.subplots(nrows=2)
### axs[0].set_title('Ranger HGT confidences')
### sns.distplot(transfer_supports,     ax=axs[0])
### axs[1].set_title('Ranger mapping consistencies')
### sns.distplot(mapping_consistencies, ax=axs[1])
### fig.tight_layout()
### fig.savefig('ranger_support_metrics.pdf', dpi=600)
###
### farthest_leaf_from_root   = named_reference_tree.get_farthest_leaf(   topology_only=True)[0]
### longest_possible_distance = farthest_leaf_from_root.get_farthest_node(topology_only=True)[1]
### bipartition_distances        = []
### bipartition_distances_normed = []
### donor_distance_to_root       = []
### for transfers in final_transfers.values():
###     for transfer in transfers:
###         donor, recipient = re.search('Mapping --> (\S+), Recipient --> (\S+)$', transfer, re.M).groups()
###         recipient_branch = named_reference_tree.search_nodes(name=recipient)[0]
###         donor_branch     = named_reference_tree.search_nodes(name=donor    )[0]
###         bipartition_distances.append(donor_branch.get_distance(recipient_branch, topology_only=True))
###         bipartition_distances_normed.append(donor_branch.get_distance(recipient_branch, topology_only=True)/longest_possible_distance)
###         donor_distance_to_root.append(named_reference_tree.get_distance(donor_branch, topology_only=True))
###
### fig, axs = plt.subplots(nrows=2)
### axs[0].set_title('# of bipartitions from root to donor')
### sns.distplot(donor_distance_to_root, ax=axs[0])
### axs[1].set_title('# of bipartitions between donor/recipient (normalized)')
### sns.distplot(bipartition_distances_normed, ax=axs[1])
### fig.tight_layout()
### fig.savefig('hgt_distances.pdf', dpi=600)
