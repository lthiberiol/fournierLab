#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
import multiprocessing
import pickle as pkl
import pandas as pd
import numpy as np
import random
from commands import getoutput
import itertools
import seaborn as sns
from matplotlib import pyplot as plt

ncbi = ete3.NCBITaxa()

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

rf_values = {'donor':[], 'recipient':[]}
rf_values_max = {'donor':[], 'recipient':[]}
def visualize_tree(group, transfers, reference_tree, taxonomy_df, output_folder='.'):
    tmp_names = {}

    gene_tree = transfers['gene_tree'].copy('deepcopy')
    count = -1
    for donor, recipient in transfers['donor/recipient']:
        count += 1

        reticulation = gene_tree.search_nodes(name='reticulation', hgt_event=count)
        reticulation = reticulation[0] if reticulation else None

        try:
            donor_in_gene_tree     = gene_tree.search_nodes(name='donor'    , hgt_event=count)[0]
            recipient_in_gene_tree = gene_tree.search_nodes(name='recipient', hgt_event=count)[0]
        except IndexError:
            continue

        #
        # set reticulation name
        if reticulation:
            tmp_id = str(random.random())
            tmp_names[tmp_id] = '[&support=%.2f,hgt_role="reticulation"]' %reticulation.support
            reticulation.name  = tmp_id

        donor_branch = reference_tree.search_nodes(name=donor)[0]
        if 'annotation' not in donor_branch.features:
            donor_branch.add_feature('annotation', [])
        recipient_branch = reference_tree.search_nodes(name=recipient)[0]
        if 'annotation' not in recipient_branch.features:
            recipient_branch.add_feature('annotation', [])

        donor_branch.annotation.append('%s_role="donor%i"' %(group, count))
        recipient_branch.annotation.append('%s_role="recipient%i"' %(group, count))

        #
        # set donor name
        tmp_id                            = str(random.random())
        tmp_names[tmp_id]                 = '[&support=%.2f,hgt_role="donor",event_count=%i]' %(donor_in_gene_tree.support, count)
        donor_in_gene_tree.name = tmp_id
        try:
            rf, rf_max, names, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = donor_in_gene_tree.robinson_foulds(donor_branch, attr_t1='genome_name')
            if rf_max:
                rf_normalized = rf/float(rf_max)
            else:
                rf_normalized = 0.0
        except ete3.coretype.tree.TreeError:
            rf_normalized = np.nan
            rf_max        = np.nan
        rf_values['donor'].append(rf_normalized)
        rf_values_max['donor'].append(rf_max)

        #
        # set recipient name
        tmp_id                                 = str(random.random())
        tmp_names[tmp_id]                      = '[&support=%.2f,hgt_role="recipient",event_count=%i]' %(recipient_in_gene_tree.support, count)
        recipient_in_gene_tree.name  = tmp_id
        try:
            rf, rf_max, names, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = recipient_in_gene_tree.robinson_foulds(recipient_branch, attr_t1='genome_name')
            if rf_max:
                rf_normalized = rf/float(rf_max)
            else:
                rf_normalized = 0.0
        except ete3.coretype.tree.TreeError:
            rf_normalized = np.nan
            rf_max        = np.nan
        rf_values['recipient'].append(rf_normalized)
        rf_values_max['recipient'].append(rf_max)

    out  = open('%s/%s.Figtree.tree' %(output_folder, group), 'wb')
    out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(gene_tree))
    for node in gene_tree.traverse():
        if node.is_leaf():
            taxid         = taxonomy_df.loc[node.genome_name, 'taxid']
            lineage       = {j: i for i, j in ncbi.get_rank(ncbi.get_lineage(taxid)).items()}
            lineage_names = ncbi.get_taxid_translator(lineage.values())

            out.write('\t%s ' %(node.name))
            comment = ['source_name="%s"' %ncbi.get_taxid_translator([taxid]).itervalues().next()]
            for rank in ['class', 'phylum', 'order']:
                if rank in lineage:
                    comment.append('tax_%s="%s"' %(rank, lineage_names[lineage[rank]]))
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

def jaccard(i,j):
    j_index = len(i.intersection(j))/float(len(i.union(j)))
    return 1-j_index

def find_most_similar_branch(branch, tree, jaccard_threshold=0.25, exceptions=None):
    branch_leaves     = set(branch.get_leaf_names())
    tree_genomes      = set([leaf.genome_name for leaf in tree.get_leaves()])
    intersect_genomes = branch_leaves.intersection(tree_genomes)

    most_similar_branch = tree
    minimal_rf_distance = 100000000 #any absurdly large number would work...
    for node in tree.traverse():
        if node.is_leaf():
            continue

        if not intersect_genomes.issubset([tmp.genome_name for tmp in node.get_leaves()]):
            continue

        try:
            rf, rf_max, names, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = node.robinson_foulds(branch, attr_t1='genome_name')
        except ete3.coretype.tree.TreeError, error_message:
            if str(error_message) != "'Duplicated items found in source tree'":
                print '**ERROR! different error message than expected, should be assessed: %s' %error_message
            continue

        if rf < minimal_rf_distance:
            minimal_rf_distance = rf
            most_similar_branch = node

    most_similar_branch_genomes = set([leaf.genome_name for leaf in most_similar_branch.get_leaves()])

    if exceptions is not None:
        exceptions_genomes = set([leaf.genome_name for leaf in exceptions.get_leaves()])
        most_similar_branch_genomes.difference_update(exceptions_genomes)

    j_distance = jaccard(most_similar_branch_genomes, branch_leaves)
    if j_distance > jaccard_threshold:
        return None,None
    return (most_similar_branch, j_distance)

def assess_transfers(transfer_file, reference_tree=named_reference_tree, debug=False):
    if not transfer_file.endswith('.pkl'):
        return None

    selected_transfers = []
    jaccard = {'donor':[], 'recipient':[]}
    group              = transfer_file.replace('.pkl', '')
    try:
        gene_tree          = ete3.Tree('/work/Alphas_and_Cyanos/ranger_input_trees/%s.tree' %group)
    except ete3.parser.newick.NewickError:
        return None

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

        if reticulation_with_support.support < 95 and reticulation_with_support.get_ancestors()[0].support < 95:
            continue

        donor_branch     = reference_tree.search_nodes(name=donor    )[0].copy(method='deepcopy')
        recipient_branch = reference_tree.search_nodes(name=recipient)[0].copy(method='deepcopy')

        recipient_branch_in_reticulation, recipient_j_distance = find_most_similar_branch(recipient_branch, reticulation_with_support, jaccard_threshold=0.5)
        donor_branch_in_reticulation, donor_j_distance     = find_most_similar_branch(donor_branch, reticulation_with_support, jaccard_threshold=0.5, exceptions=recipient_branch_in_reticulation)

        if not recipient_branch_in_reticulation or not donor_branch_in_reticulation:
            continue

        if donor_branch_in_reticulation in recipient_branch_in_reticulation.get_ancestors():
            print '%s: recipient nested into donor' %group
            reticulation_with_support.name = 'reticulation'
            reticulation_with_support.add_feature('hgt_event', len(selected_transfers))

            donor_branch_in_reticulation.name = 'donor'
            donor_branch_in_reticulation.add_feature('hgt_event', len(selected_transfers))

            recipient_branch_in_reticulation.name = 'recipient'
            recipient_branch_in_reticulation.add_feature('hgt_event', len(selected_transfers))

            selected_transfers.append((donor, recipient))
            jaccard['donor'].append(donor_j_distance)
            jaccard['recipient'].append(recipient_j_distance)

    return {group:{'donor/recipient':selected_transfers,'gene_tree':gene_tree, 'jaccard':jaccard}}

aggregated_folder = 'aggregated/90_threshold'
with cd(aggregated_folder):
    pool = multiprocessing.Pool(processes=6)
    results = pool.map(assess_transfers, os.listdir('.'))
    pool.close()
    pool.join()

jaccard_dists = {'donor':[], 'recipient':[]}
usefull_results = {}
for entry in results:
    if entry and entry.values()[0]['donor/recipient']:
        usefull_results.update(entry)
        for yeah in 'donor recipient'.split():
            jaccard_dists[yeah].extend(entry.values()[0]['jaccard'][yeah])

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
header = 'assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'.split()
genbank_summary                     = pd.read_table( 'assembly_summary_genbank.txt', comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str} )
genbank_summary['refseq_category']  = genbank_summary['refseq_category'].str.lower()
genbank_summary['assembly_level']   = genbank_summary['assembly_level'].str.lower()
genbank_summary['genome_rep']       = genbank_summary['genome_rep'].str.lower()
genbank_summary.set_index('assembly_accession', inplace=True)
genbank_summary.index               = [re.sub('\.\d+$', '', index).replace('_', '') for index in genbank_summary.index]
genbank_summary                     = genbank_summary.reindex(named_reference_tree.get_leaf_names())

tree_visualization_folder = '/work/Alphas_and_Cyanos/index_transfer_trees/90_threshold'
species_tree_to_color = named_reference_tree.copy(method='deepcopy')
for group, transfers in usefull_results.items():
    species_tree_to_color = visualize_tree(group, transfers, species_tree_to_color, genbank_summary, output_folder=tree_visualization_folder)

out  = open('%s/species_tree.Figtree.tree' %tree_visualization_folder, 'wb')
out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(species_tree_to_color))
for node in species_tree_to_color.traverse():
    if node.is_leaf():
        taxid = genbank_summary.loc[node.name, 'taxid']
        lineage = {j: i for i, j in ncbi.get_rank(ncbi.get_lineage(taxid)).items()}
        lineage_names = ncbi.get_taxid_translator(lineage.values())

        out.write('\t%s ' % (node.name))
        comment = ['source_name="%s"' % ncbi.get_taxid_translator([taxid]).itervalues().next()]
        for rank in ['class', 'phylum', 'order']:
            if rank in lineage:
                comment.append('tax_%s="%s"' % (rank, lineage_names[lineage[rank]]))
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
fig, axs = plt.subplots(nrows=2)
axs[0].set_title('Donor branch normalized Robinson-Foulds distances')
sns.distplot([rf for rf in rf_values['donor'] if pd.notnull(rf)],     ax=axs[0])
axs[1].set_title('Recipient branch normalized Robinson-Foulds distances')
sns.distplot([rf for rf in rf_values['recipient'] if pd.notnull(rf)], ax=axs[1])
fig.tight_layout()
fig.savefig('rf_distances_norm.pdf', dpi=600)
###
fig, axs = plt.subplots(nrows=2)
axs[0].set_title('Donor branch size variation')
sns.distplot(jaccard_dists['donor'],     ax=axs[0])
axs[1].set_title('Recipient branch size variation')
sns.distplot(jaccard_dists['recipient'], ax=axs[1])
fig.tight_layout()
fig.savefig('jaccard_branch_distances.pdf', dpi=1200)
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
