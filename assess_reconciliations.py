#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
import multiprocessing
import pickle as pkl
import pandas as pd
import random
from commands import getoutput
import itertools
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

def visualize_tree(group, transfers, reference_tree, taxonomy_df, output_folder='.'):
    tmp_names = {}
    gene_tree = ete3.Tree('/work/Alphas_and_Cyanos/ranger_input_trees/%s.tree' %group)
    for leaf in gene_tree.get_leaves():
        leaf.name = re.sub('\.\d+', '', leaf.name.replace('GCA_', 'GCA'))
        leaf.add_feature('genome_name', leaf.name.split('_')[0])

    count = 0
    for [donor, recipient, reticulation, topology_id] in transfers:
        count += 1

        is_it_monophyletic, clade_type, fucking_up = gene_tree.check_monophyly(reticulation.get_leaf_names(), 'name', unrooted=False)
        if is_it_monophyletic:
            reticulation_with_support = gene_tree.get_common_ancestor(reticulation.get_leaf_names())
        else:
            gene_tree.set_outgroup(fucking_up.pop())
            is_it_monophyletic, clade_type, fucking_up = gene_tree.check_monophyly(reticulation.get_leaf_names(), 'name', unrooted=False)
            if not is_it_monophyletic:
                raise KeyError("Reticulation is not monphyletic, something is really broken...")
            reticulation_with_support = gene_tree.get_common_ancestor(reticulation.get_leaf_names())

        #
        # set reticulation name
        tmp_id = str(random.random())
        tmp_names[tmp_id] = '[&support=%.2f,hg_role="reticulation"]' %reticulation_with_support.support
        reticulation_with_support.name  = tmp_id

        donor_branch      = reference_tree.search_nodes(name=donor    )[0]
        if 'annotation' not in donor_branch.features:
            donor_branch.add_feature('annotation', [])
        recipient_branch      = reference_tree.search_nodes(name=recipient)[0]
        if 'annotation' not in recipient_branch.features:
            recipient_branch.add_feature('annotation', [])
        if count == 1:
            donor_branch.annotation.append('%s_role="donor"' %group)
            recipient_branch.annotation.append('%s_role="recipient"' %group)
        else:
            donor_branch.annotation.append('%s_role="donor%i"' %(group, count))
            recipient_branch.annotation.append('%s_role="recipient%i"' %(group, count))

        donor_leaves     = [leaf.name for leaf in reticulation_with_support.get_leaves() if leaf.genome_name in donor_branch.get_leaf_names()    ]
        recipient_leaves = [leaf.name for leaf in reticulation_with_support.get_leaves() if leaf.genome_name in recipient_branch.get_leaf_names()]

        #
        # set donor name
        donor_branch_in_reticulation      = reticulation_with_support.get_common_ancestor(    donor_leaves)
        tmp_id                            = str(random.random())
        tmp_names[tmp_id]                 = '[&support=%.2f,hgt_role="donor",event_count=%i]' %(donor_branch_in_reticulation.support, count)
        donor_branch_in_reticulation.name = tmp_id

        #
        # set recipient name
        recipient_branch_in_reticulation       = reticulation_with_support.get_common_ancestor(recipient_leaves)
        tmp_id                                 = str(random.random())
        tmp_names[tmp_id]                      = '[&support=%.2f,hgt_role="recipient",event_count=%i]' %(recipient_branch_in_reticulation.support, count)
        recipient_branch_in_reticulation.name  = tmp_id

    out  = open('%s/%s.Figtree.tree' %(output_folder, group), 'wb')
    out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(gene_tree))
    for node in gene_tree.traverse():
        if node.is_leaf():
            out.write('\t%s ' %(node.name))
            comment = []
            for rank in ['organism_name', 'class', 'phylum', 'order']:
                comment.append('tax_%s="%s"' %(rank, taxonomy_df.loc[node.genome_name, rank]))
            out.write('[&%s]\n' %' '.join(comment))

        else:
            if node.support and not node.name:
                node.name = '[&support=%.2f]' %node.support

    newick_text = gene_tree.write(format=1)
    newick_text = re.sub('_&support_(\d+\.\d\d)_', '[&support=\\1]', newick_text)
    for key, value in tmp_names.items():
        newick_text = newick_text.replace(key, value)
    out.write(';\nend;\n')
    out.write('begin trees;\n\ttree tree_1 = [&R] %s\nend;' %newick_text)
    out.close()

    return reference_tree

def assess_compatibilities(donor, recipient, reticulation):
    #
    # ignore reticulation if there is a duplication in it
    if len(set([leaf.genome_name for leaf in reticulation.get_leaves()])) < len(reticulation):
        genomes_in_reticulation = [leaf.genome_name for leaf in reticulation.get_leaves()]
        non_duplicated_leaves   = [leaf.name for leaf in reticulation.get_leaves() if genomes_in_reticulation.count(leaf.genome_name) == 1]
        reticulation.prune(non_duplicated_leaves)

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

        if not assess_compatibilities(donor_branch, recipient_branch, reticulation_with_support.copy(method='deepcopy')):
            continue

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

aggregated_folder = 'aggregated/80_threshold'
with cd(aggregated_folder):
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

with cd(aggregated_folder):
    out = open('index_transfers.pkl', 'wb')
    pkl.dump(usefull_results, out)
    out.close()

########################################################################################################################
#                                                                                                                      #
# Plot HGT data                                                                                                        #
#                                                                                                                      #
########################################################################################################################
most_distant_internal_node = ''
distance_from_root         = 0
for node in named_reference_tree.traverse():
    if node.is_leaf():
        continue

    if named_reference_tree.get_distance(node, topology_only=True) > distance_from_root:
        most_distant_internal_node = node
        distance_from_root         = named_reference_tree.get_distance(most_distant_internal_node, topology_only=True)

most_distant_internal_node2 = ''
longest_distance            = 0
for node in named_reference_tree.traverse():
    if node.is_leaf():
        continue

    if most_distant_internal_node.get_distance(node, topology_only=True) > longest_distance:
        most_distant_internal_node2 = node
        longest_distance            = most_distant_internal_node.get_distance(most_distant_internal_node2, topology_only=True)

aggregated_folder = 'aggregated/90_threshold'
usefull_results   = pkl.load(open('%s/index_transfers.pkl' %aggregated_folder))

transfer_distances       = []
donor2root_distances     = []
recipient2root_distances = []
for group, transfers in usefull_results.items():
    for [donor, recipient, reticulation, topology_id] in transfers:

        donor_branch     = named_reference_tree.search_nodes(name=donor    )[0]
        recipient_branch = named_reference_tree.search_nodes(name=recipient)[0]

        transfer_distances.append(              donor_branch.get_distance(recipient_branch, topology_only=True))
        donor2root_distances.append(    named_reference_tree.get_distance(donor_branch,     topology_only=True))
        recipient2root_distances.append(named_reference_tree.get_distance(recipient_branch, topology_only=True))

fig, axs = plt.subplots(nrows=2)
axs[0].set_title('Donor-Recipient bipartition distances')
sns.distplot(transfer_distances, ax=axs[0])
axs[1].set_title('Donor-Recipient bipartition distances (normalized)')
sns.distplot([distance/longest_distance for distance in transfer_distances], ax=axs[1])
fig.tight_layout()
fig.savefig('%s/transfer_distance.pdf' %aggregated_folder, dpi=600)

plot = sns.jointplot(x=pd.Series(name='Donor distances to root', data=[distance/distance_from_root for distance in donor2root_distances]), y=pd.Series(name='Recipient distances to root', data=[distance/distance_from_root for distance in recipient2root_distances]))
plot.savefig('%s/distances_to_root.pdf' %aggregated_folder, dpi=600)

########################################################################################################################
#                                                                                                                      #
# HGT vizualitation with FigTree                                                                                       #
#                                                                                                                      #
########################################################################################################################
genome_lineages = pd.read_table( '/work/lbi_backup/fournierLab/lineages_from_genbank_summary.tab', sep='\t', index_col=0, dtype=str )
genome_lineages.index = genome_lineages.index.astype( str )
genome_lineages['taxid'] = genome_lineages.index

genbank_summary                     = pd.read_table( '/work/lbi_backup/assembly_summary_genbank.txt', dtype={'taxid':str, 'infraspecific_name':str} )
genbank_summary['refseq_category']  = genbank_summary['refseq_category'].str.lower()
genbank_summary['assembly_level']   = genbank_summary['assembly_level'].str.lower()
genbank_summary['genome_rep']       = genbank_summary['genome_rep'].str.lower()
genbank_summary                     = genbank_summary.merge( genome_lineages[['taxid', 'class', 'phylum', 'order']], how='left', on='taxid' )
genbank_summary.set_index('assembly_accession', inplace=True)
genbank_summary.index               = [re.sub('\.\d+$', '', index).replace('_', '') for index in genbank_summary.index]
genbank_summary                     = genbank_summary.reindex(named_reference_tree.get_leaf_names())

tree_visualization_folder = '/work/Alphas_and_Cyanos/index_transfer_trees/80_threshold'
species_tree_to_color = named_reference_tree.copy(method='deepcopy')
for group, transfers in usefull_results.items():
    species_tree_to_color = visualize_tree(group, transfers, species_tree_to_color, genbank_summary, output_folder=tree_visualization_folder)

out  = open('%s/species_tree.Figtree.tree' %tree_visualization_folder, 'wb')
out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(species_tree_to_color))
for node in species_tree_to_color.traverse():
    if node.is_leaf():
        out.write('\t%s ' %(node.name))
        comment = []
        for rank in ['organism_name', 'class', 'phylum', 'order']:
            comment.append('tax_%s="%s"' %(rank, genbank_summary.loc[node.name, rank]))
        out.write('[&%s]\n' %' '.join(comment))

    else:
        if 'annotation' in node.features:
            node.name = '[&%s]' %','.join(node.annotation)

newick_text = species_tree_to_color.write(format=1)
newick_text = newick_text.replace('role_"', 'role="')
newick_text = newick_text.replace('"_:', '"]:')
newick_text = newick_text.replace(')_&', ')[&')
newick_text = re.sub('_(\d{6}_role="(?:donor\d?|recipient\d?)")', ',\\1', newick_text)
out.write(';\nend;\n')
out.write('begin trees;\n\ttree tree_1 = [&R] %s\nend;' %newick_text)
out.close()

########################################################################################################################
### fig, axs = plt.subplots(nrows=2)
### axs[0].set_title('Donor branch Robinson-Foulds distances')
### sns.distplot(rf_values['donor'],     ax=axs[0])
### axs[1].set_title('Recipient branch Robinson-Foulds distances')
### sns.distplot(rf_values['recipient'], ax=axs[1])
### fig.tight_layout()
### fig.savefig('rf_distances-index_transfer.pdf', dpi=600)
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
