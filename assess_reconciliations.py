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
from time import time
import plotly
import plotly.plotly as ptl
from plotly import graph_objs as go
ptl.sign_in('lthiberiol', 'm15ikp59lt')

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

def visualize_tree((group, transfers), taxonomy_df, reference_tree, output_folder='/work/Alphas_and_Cyanos/index_transfer_trees/90_threshold'):
    tmp_names = {}

    gene_tree = ete3.Tree('ranger_input_trees/{group}.tree'.format(group=group), format=0)
    for leaf in gene_tree.get_leaves():
        genome_name, gene_name = leaf.name.split('_')
        leaf.add_feature('genome', genome_name)

    count = 0
    for mappings in transfers:
        is_it_monophyletic, clade_type, fucking_up = gene_tree.check_monophyly(mappings['reticulation'].get_leaf_names(), 'name', unrooted=False)
        if is_it_monophyletic:
            donor_map = gene_tree.get_common_ancestor(mappings['reticulation'].get_leaf_names())
        else:
            gene_tree.set_outgroup(fucking_up.pop())
            donor_map = gene_tree.get_common_ancestor(mappings['reticulation'].get_leaf_names())

        if donor_map.get_topology_id() != mappings['donor_map']:
            raise ValueError('Topology ids not matching, {group} between {donor} and {recipient}' .format(group=group, donor=mappings['donor'], recipient=mappings['recipient']))

        for node in donor_map.traverse():
            if not node.is_leaf() and node.get_topology_id() == mappings['recipient_map']:
                recipient_map = node
                break

        donor_branch = reference_tree.search_nodes(name=mappings[ 'donor'])[0]
        if 'annotation' not in donor_branch.features:
            donor_branch.add_feature('annotation', {'donor':[], 'recipient':[]})
        recipient_branch = reference_tree.search_nodes(name=mappings['recipient'])[0]
        if 'annotation' not in recipient_branch.features:
            recipient_branch.add_feature('annotation', {'donor':[], 'recipient':[]})

        donor_branch.annotation['donor'].append([group, count])
        recipient_branch.annotation['recipient'].append([group, count])

        #
        # set donor name
        tmp_id                            = str(random.random())
        tmp_names[tmp_id]                 = '[&support=%.2f,hgt_role="donor%i",hgt_count=%i]' %(donor_map.support, count, count)
        donor_map.name = tmp_id

        #
        # set recipient name
        tmp_id                                 = str(random.random())
        tmp_names[tmp_id]                      = '[&support=%.2f,hgt_role="recipient%i",hgt_count=%i]' %(recipient_map.support, count, count)
        recipient_map.name  = tmp_id

        count += 1
        del(donor_map, recipient_map)

    out  = open('%s/%s.Figtree.tree' %(output_folder, group), 'wb')
    out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(gene_tree))
    for node in gene_tree.traverse():
        if node.is_leaf():
            taxid         = taxonomy_df.loc[node.genome, 'taxid']
            lineage       = {j: i for i, j in ncbi.get_rank(ncbi.get_lineage(taxid)).items()}
            lineage_names = ncbi.get_taxid_translator(lineage.values())

            out.write('\t%s ' %(node.name))
            comment = ['source_name="%s"' %ncbi.get_taxid_translator([taxid]).itervalues().next()]
            for rank in ['phylum', 'class', 'order', 'family', 'genus']:
                if rank in lineage:
                    comment.append('tax_%s="%s"' %(rank, lineage_names[lineage[rank]]))
            out.write('[&%s]\n' %' '.join(comment))

        else:
            if node.support and not node.name:
                node.name = '[&support=%.2f]' %node.support

    rooted_tree = ete3.Tree(open('reconciliations/{group}/{group}.reconciliation1'.format(group=group)).readlines()[7], format=1)
    if rooted_tree.children[0].is_leaf():
        gene_tree.set_outgroup(rooted_tree.children[0].name)
    elif rooted_tree.children[1].is_leaf():
        gene_tree.set_outgroup(rooted_tree.children[1].name)
    else:
        is_it_monophyletic, clade_type, fucking_up = gene_tree.check_monophyly(rooted_tree.children[0].get_leaf_names(), 'name')
        if is_it_monophyletic:
            gene_tree.set_outgroup(gene_tree.get_common_ancestor(rooted_tree.children[0].get_leaf_names()))
        else:
            gene_tree.set_outgroup(gene_tree.get_common_ancestor(rooted_tree.children[1].get_leaf_names()))

    newick_text = gene_tree.write(format=1)
    newick_text = re.sub('_&support_(\d+\.\d\d)_', '[&support=\\1]', newick_text)
    for key, value in tmp_names.items():
        newick_text = newick_text.replace(key, value)
    out.write(';\nend;\n')
    out.write('begin trees;\n\ttree tree_1 = [&R] %s\nend;' %newick_text)
    out.close()

    return reference_tree

def rename_branches(reconciliation_file):
    reconciliation = open(reconciliation_file).read()
    tree = ete3.Tree(reconciliation.split('\n')[7], format=1)
    branches = re.findall('^(m\d+) = LCA\[(\S+), (\S+)\]', reconciliation, re.M)
    for node in tree.traverse():
        if node.is_leaf():
            continue
        else:
            node.name = ''
    for name, leaf1, leaf2 in branches:
        node = tree.get_common_ancestor(leaf1, leaf2)
        if node.name:
            print node.get_leaves_by_name()
            break
        node.name = name
    return tree

def assess_transfers(transfer_file, species_tree=named_reference_tree, working_dir='/work/Alphas_and_Cyanos'):
    if not transfer_file.endswith('.pkl'):
        return None

    group              = transfer_file.replace('.pkl', '')
    transfer_data = pkl.load(open(transfer_file))
    mapped_nodes  = {}

    support_tree = ete3.Tree('%s/ranger_input_trees/%s.tree' %(working_dir, group))

    with cd('%s/reconciliations/%s' %(working_dir, group)):
        reconciliation_count = 1
        descendants = {}
        for donor, recipient, reticulation, topology_id in transfer_data:
            recipient_branch       = species_tree.search_nodes(name=recipient)[0]
            descendants[recipient] = [descendant.name for descendant in recipient_branch.get_descendants() if not descendant.is_leaf()]
            descendants[recipient].append(recipient)

        while True:
            tmp_tree = rename_branches('%s.reconciliation%i' %(group, reconciliation_count))
            if tmp_tree.children[0].is_leaf():
                support_tree.set_outgroup(tmp_tree.children[0].name)
            elif tmp_tree.children[1].is_leaf():
                support_tree.set_outgroup(tmp_tree.children[1].name)
            else:
                is_it_monophyletic, clade_type, fucking_up = support_tree.check_monophyly(tmp_tree.children[0].get_leaf_names(), 'name')
                if is_it_monophyletic:
                    support_tree.set_outgroup(support_tree.get_common_ancestor(tmp_tree.children[0].get_leaf_names()))
                else:
                    support_tree.set_outgroup(support_tree.get_common_ancestor(tmp_tree.children[1].get_leaf_names()))

            if support_tree.get_topology_id() != tmp_tree.get_topology_id():
                raise ValueError('%s support and named trees not matching: reconciliation_count=%i' %(group, reconciliation_count))

            mapped_donors       = {topology_id:'' for donor, recipient, reticulation, topology_id in transfer_data}
            map_count           = 0
            for node in tmp_tree.traverse():
                if not node.is_leaf() and node.get_topology_id() in mapped_donors:
                    if support_tree.get_common_ancestor(node.get_leaf_names()).support >= 95:
                        mapped_donors[node.get_topology_id()] = node.name
                    map_count += 1
                if map_count == len(mapped_donors):
                    break

            for donor, recipient, reticulation, topology_id in transfer_data:
                index = tuple([donor, recipient, reticulation, topology_id])
                if not mapped_donors[topology_id]:
                    continue

                mapped_donor = tmp_tree.search_nodes(name=mapped_donors[topology_id])[0]

                mapping_child1 = getoutput('grep -E "^%s = LCA\[\S+, \S+\]: \S+, Mapping --> (%s)" %s.reconciliation{%i..%i}' %(mapped_donor.children[0].name, '|'.join(descendants[recipient]), group, reconciliation_count, reconciliation_count+19))
                mapping_child2 = getoutput('grep -E "^%s = LCA\[\S+, \S+\]: \S+, Mapping --> (%s)" %s.reconciliation{%i..%i}' %(mapped_donor.children[1].name, '|'.join(descendants[recipient]), group, reconciliation_count, reconciliation_count+19))

                mapping = mapping_child1
                child   = mapped_donor.children[0]
                if len(mapping_child1.split('\n')) < len(mapping_child2.split('\n')):
                    mapping = mapping_child2
                    child   = mapped_donor.children[1]

                if len(mapping.split('\n')) > 20:
                    raise ValueError('something is wrong with group %s, between %s and %s' %(group, donor, recipient))

                if index not in mapped_nodes:
                    mapped_nodes[index] = []
                mapped_nodes[index].append((topology_id, child.get_topology_id()))

            reconciliation_count += 20
            if not os.path.isfile('%s.reconciliation%i' %(group, reconciliation_count)):
                break

    selected_transfers = []
    for index, mappings in mapped_nodes.items():
        if len(set(mappings)) != 1:
            break
            raise ValueError('Multiple reconciliations mappings for group %s, between %s and %s' %(group, index[0], index[1]))
        else:
            selected_transfers.append({'donor':index[0], 'recipient':index[1], 'reticulation':index[2], 'donor_map':mappings[0][0], 'recipient_map':mappings[0][1]})

    return {group:selected_transfers}

aggregated_folder = 'aggregated/90_threshold'
with cd(aggregated_folder):
    pool = multiprocessing.Pool(processes=4)
    results = pool.map(assess_transfers, os.listdir('.'))
    pool.close()
    pool.join()

    out = open('index_transfers_candidates.pkl', 'wb')
    pkl.dump(usefull_results, out)
    out.close()

usefull_results = {}
for entry in results:
    key, value = entry.items()[0]
    if value:
        usefull_results[key] = value

out = open('index_transfers_dtl_distances.pkl', 'wb')
pkl.dump(usefull_results, out)
out.close()

def assess_dtl_dist((group, transfer_data), working_dir='/home/thiberio/Alphas_and_Cyanos'):
    dtl_distances = {'donor':[], 'recipient':[]}
    donor_trees     = []
    recipient_trees = []
    for transfer in transfer_data:
        if transfer['reticulation'].children[1].get_topology_id() == transfer['recipient_map']:
            donor_trees.append(    transfer['reticulation'].children[0].write(format=9))
            recipient_trees.append(transfer['reticulation'].children[1].write(format=9))
        elif transfer['reticulation'].children[0].get_topology_id() == transfer['recipient_map']:
            donor_trees.append(    transfer['reticulation'].children[1].write(format=9))
            recipient_trees.append(transfer['reticulation'].children[0].write(format=9))
        else:
            raise ValueError('wrong tree mapping on group %s!' %group)

    #
    # donor compatibility assessment
    os.system( 'cp %s/ranger_input_trees/species_tree.template tmp_ranger-%s.input' %(working_dir, multiprocessing.current_process().name))
    out = open('tmp_ranger-%s.input' %multiprocessing.current_process().name, 'a')
    out.write('\n'.join(donor_trees))
    out.close()
    os.system('/home/thiberio/ranger/CorePrograms/Ranger-DTL.linux -q -i tmp_ranger-%s.input -o tmp_ranger-%s.output' %(multiprocessing.current_process().name, multiprocessing.current_process().name))
    dtl_distances['donor'].extend([int(reconciliation_cost) for reconciliation_cost in re.findall('^The minimum reconciliation cost is: (\d+)', open('tmp_ranger-%s.output' %multiprocessing.current_process().name).read(), re.M)])

    #
    # recipient compatibility assessment
    os.system( 'cp %s/ranger_input_trees/species_tree.template tmp_ranger-%s.input' %(working_dir, multiprocessing.current_process().name))
    out = open('tmp_ranger-%s.input' %multiprocessing.current_process().name, 'a')
    out.write('\n'.join(recipient_trees))
    out.close()
    os.system('/home/thiberio/ranger/CorePrograms/Ranger-DTL.linux -q -i tmp_ranger-%s.input -o tmp_ranger-%s.output' %(multiprocessing.current_process().name, multiprocessing.current_process().name))
    dtl_distances['recipient'].extend([int(reconciliation_cost) for reconciliation_cost in re.findall('^The minimum reconciliation cost is: (\d+)', open('tmp_ranger-%s.output' %multiprocessing.current_process().name).read(), re.M)])

    return {group:dtl_distances}

pool      = multiprocessing.Pool(processes=8)
distances = pool.map(assess_dtl_dist, usefull_results.items())
pool.close()
pool.join()

dtl_distances = {}
for entry in distances:
    dtl_distances.update(entry)

transfer_distances       = {group:[] for group in usefull_results.keys()}
donor2root_distances     = {group:[] for group in usefull_results.keys()}
recipient2root_distances = {group:[] for group in usefull_results.keys()}
for group, transfers in usefull_results.items():
    for transfer in transfers:

        donor_branch     = named_reference_tree.search_nodes(name=transfer['donor']    )[0]
        recipient_branch = named_reference_tree.search_nodes(name=transfer['recipient'])[0]

        transfer_distances[group].append(              donor_branch.get_distance(recipient_branch))
        donor2root_distances[group].append(    named_reference_tree.get_distance(donor_branch,     topology_only=True))
        recipient2root_distances[group].append(named_reference_tree.get_distance(recipient_branch, topology_only=True))

out = open('index_transfers_dtl_distances.pkl', 'wb')
pkl.dump(distances, out)
out.close()
########################################################################################################################
#                                                                                                                      #
# Plot HGT data                                                                                                        #
#                                                                                                                      #
########################################################################################################################
x, y = [],[]
for group in dtl_distances.keys():
    for i,j in zip(dtl_distances[group]['donor'], donor2root_distances[group]):
        x.append(i)
        y.append(j)
plot = sns.jointplot(x=pd.Series(name='Donor subtrees DTL distances', data=x), y=pd.Series(name='Donor node "bipartition distance" to root', data=y), color='black', joint_kws={'s':3, 'alpha':0.3})
plot.savefig('%s/dtl_VS_rootDistance-donors.pdf' %aggregated_folder, dpi=600)

x, y = [],[]
for group in dtl_distances.keys():
    for i,j in zip(dtl_distances[group]['recipient'], recipient2root_distances[group]):
        x.append(i)
        y.append(j)
plot = sns.jointplot(x=pd.Series(name='Recipient subtrees DTL distances', data=x), y=pd.Series(name='Recipient node "bipartition distance" to root', data=y), color='black', joint_kws={'s':3, 'alpha':0.3})
plot.savefig('%s/dtl_VS_rootDistance-recipients.pdf' %aggregated_folder, dpi=600)

tracer = {'color':'black', 'x':[], 'y':[], 'text':[], 'linecolor':'transparent', 'marker':'circle'}
for group in dtl_distances.keys():
    for position in range(len(dtl_distances[group]['recipient'])):
        tracer['x'].append(dtl_distances[group]['recipient'][position])
        tracer['y'].append(transfer_distances[group][position])
        tracer['text'].append('%s-#%i' %(group, position))

plot_data = go.Data([go.Scatter(x=tracer['x'], y=tracer['y'], mode='markers', text=tracer['text'], name='Subtree differences VS transfer distances', hoverinfo='text',
                        marker=dict(size=10, color=tracer['color'], symbol=tracer['marker'], opacity=0.75))])
layout    = go.Layout(title='Subtree differences VS transfer distances', hovermode='closest', width=1200, height=1000, xaxis=dict(title='Recipient subtrees DTL distances'), yaxis=dict(title='Bipartition distance between donor/recipient'))
fig       = go.Figure(data=plot_data, layout=layout)
plot      = plotly.offline.plot(fig, filename='subtreeDiff_VS_transferDist-patristic.html', auto_open=False)
########################################################################################################################
#                                                                                                                      #
# HGT vizualitation with FigTree                                                                                       #
#                                                                                                                      #
########################################################################################################################
header = 'assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'.split()
genbank_summary                     = pd.read_table('/work/Alphas_and_Cyanos/assembly_summary_genbank.txt', comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str} )
genbank_summary['refseq_category']  = genbank_summary['refseq_category'].str.lower()
genbank_summary['assembly_level']   = genbank_summary['assembly_level'].str.lower()
genbank_summary['genome_rep']       = genbank_summary['genome_rep'].str.lower()
genbank_summary.set_index('assembly_accession', inplace=True)
genbank_summary.index               = [re.sub('\.\d+$', '', index).replace('_', '') for index in genbank_summary.index]
genbank_summary                     = genbank_summary.reindex(named_reference_tree.get_leaf_names())

tree_visualization_folder = '/work/Alphas_and_Cyanos/index_transfer_trees/90_threshold'
species_tree_to_color = named_reference_tree.copy(method='deepcopy')
for entry in usefull_results.items():
    species_tree_to_color = visualize_tree(entry, genbank_summary, species_tree_to_color, output_folder=tree_visualization_folder)

out  = open('%s/species_tree.Figtree.tree' %tree_visualization_folder, 'wb')
out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(species_tree_to_color))
tmp_names = {}
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
            tmp_id            = str(random.random())
            tmp_names[tmp_id] = {'ranger_name':node.name}
            node.name         = tmp_id
            for group, count in node.annotation['donor']:
                if group not in tmp_names[tmp_id]:
                    tmp_names[tmp_id][group] = {'donor':[], 'recipient':[]}
                tmp_names[tmp_id][group]['donor'].append(str(count))
            for group, count in node.annotation['recipient']:
                if group not in tmp_names[tmp_id]:
                    tmp_names[tmp_id][group] = {'donor':[], 'recipient':[]}
                tmp_names[tmp_id][group]['recipient'].append(str(count))

newick_text = species_tree_to_color.write(format=1)
for key, value in tmp_names.items():
    comment = []
    for group, roles in value.items():
        if group == 'ranger_name':
            continue
        donor_fragment = ''
        if roles['donor']:
            donor_fragment = 'donor%s' %','.join(roles['donor'])

        recipient_fragment = ''
        if roles['recipient']:
            recipient_fragment = 'recipient%s' %','.join(roles['recipient'])

        comment.append('%s_role="%s %s"' %(group, donor_fragment, recipient_fragment) )

    newick_text = newick_text.replace(key, '[&ranger_name="%s",%s]' %(value['ranger_name'], ','.join(comment)))
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
