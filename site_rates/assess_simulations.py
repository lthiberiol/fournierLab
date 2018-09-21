#!/usr/bin/env python2
#coding: utf-8

import ete3
import os
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.colors as colors
import matplotlib.cm as cmx
import itertools
import numpy as np
import subprocess
from Bio import SeqIO, SearchIO, AlignIO, Align, Alphabet
import re
import multiprocessing
from copy import deepcopy
from scipy.stats import spearmanr

os.chdir('/work/site_rate')
ncbi = ete3.NCBITaxa()

#
# already done!
#
header = '''assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name 
            infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name 
            submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'''.split()
assembly_summary_outgroup                    = pd.read_table('outgroups.tab', comment='#', header=None, names=header,
                                                             dtype={'taxid':str, 'infraspecific_name':str})
assembly_summary_outgroup['refseq_category'] = assembly_summary_outgroup['refseq_category'].str.lower()
assembly_summary_outgroup['assembly_level']  = assembly_summary_outgroup['assembly_level'].str.lower()
assembly_summary_outgroup['genome_rep']      = assembly_summary_outgroup['genome_rep'].str.lower()
assembly_summary_outgroup.set_index('assembly_accession', inplace=True)

sampled_genomes = pd.read_table('sampled_genomes.tab', index_col=0)

full_tree = ete3.Tree('engaging.treefile', format=1)
if full_tree.check_monophyly(assembly_summary_outgroup.index, 'name'):
    outgroup = full_tree.get_common_ancestor(assembly_summary_outgroup.index.tolist())
full_tree.set_outgroup(outgroup)

sampled_genomes = sampled_genomes.reindex(full_tree.get_leaf_names())

lineages = {}
for index, row in sampled_genomes.iterrows():
    tmp_lineage = {j: int(i) for i, j in ncbi.get_rank(ncbi.get_lineage(row.taxid)).items()}
    lineages[index] = tmp_lineage

lineages = pd.DataFrame.from_dict(lineages, dtype=int).T
lineages.drop('no rank', axis=1, inplace=True)
lineages.dropna(how='any', subset=['genus'], inplace=True)
########################################################################################################################
os.chdir('/work/site_rate/simulated_alignmentsG/')

full_tree.write(outfile='species.tree', dist_formatter='%.20f', format=5)

for category in range(1,9):
    subprocess.call(['/Users/thiberio/anaconda2/bin/raxml', '-m', 'PROTGAMMALG', '-p', '12345', '-f', 'b',
                     '-t', 'species.tree', '-z', '%i.boottrees' %category, '-n', '%i_supported' %category])

categories_supported_trees = {}
for category in range(1,9):
    tmp_tree = ete3.Tree('RAxML_bipartitions.%i_supported' %category)

    outgroup = tmp_tree.get_common_ancestor(assembly_summary_outgroup.index.tolist())
    tmp_tree.set_outgroup(outgroup)

    for node in tmp_tree.traverse():
        if not node.is_leaf():
            node.name = node.get_topology_id()

    categories_supported_trees[category] = tmp_tree.copy()

for node in full_tree.traverse():
    if node.is_leaf():
        continue

    node_topology_id = node.get_topology_id()
    for category in range(1,9):
        equivalent_node = categories_supported_trees[category].search_nodes(name=node_topology_id)[0]
        node.add_feature('category_support_%i' %category, equivalent_node.support)

#
# assess correlations for each branch
correlation_values = {category:[] for category in range(1,9)}
slope_values = {category:[] for category in range(1,9)}
for node in full_tree.traverse():
    pretip_node = False
    for child in node.children:
        if child.is_leaf():
            pretip_node = True
            break

    if not pretip_node or len(node.get_ancestors()) <= 10:
        continue

    for category in correlation_values.keys():
        tmp_supports                 = []
        tmp_num_bipartitions_to_root = []
        for ancestor in node.get_ancestors()[:-1]:
            if not ancestor.name:
                continue
            tmp_alrt, tmp_bb = [float(value) for value in ancestor.name.split('/')]
            if tmp_bb < 80 or tmp_alrt < 80:
                continue
            tmp_supports.append(eval('ancestor.category_support_%i' %category))
            tmp_num_bipartitions_to_root.append(full_tree.get_distance(ancestor, topology_only=True))

        if len(tmp_supports) < 5:
            continue
        tmp_spearman, tmp_p = spearmanr(tmp_supports, tmp_num_bipartitions_to_root)
        tmp_regression      = linregress(tmp_supports, tmp_num_bipartitions_to_root)
        correlation_values[category].append(tmp_spearman)
        slope_values[      category].append(tmp_regression.slope)
        print '%i: %.2f' %(category, tmp_regression.rvalue)

    print ''

fig, axs = plt.subplots(nrows=8, sharex=True)
for category, values in correlation_values.items():
    axs[category-1].set_title('Site-rate category %i support correlation' %category)
    sns.kdeplot(values, shade=True, ax=axs[category-1])
fig.set_size_inches(15,12)
fig.tight_layout()
fig.savefig('test.pdf', dpi=300)
plt.close()

#
# assess correlations for each branch
support_values = {category:[] for category in range(1,9)}
branch_lengths = []
for node in full_tree.traverse():
    if node.is_leaf()::row
        continue

    branch_lengths.append(node.dist)
    for category in correlation_values.keys():
        support_values[category].append(eval('node.category_support_%i' %category))

#branch_length_bins = np.linspace(0, np.max(branch_lengths), 50)
branch_length_bins = [np.percentile(branch_lengths, decile) for decile in range(20, 81, 20)]
binning            = np.digitize(branch_lengths, branch_length_bins)

fig, axs = plt.subplots(nrows=8, sharex=True)
yeah = pd.DataFrame(columns='bin category support'.split())
branch_lengths = np.asarray(branch_lengths)
for category in [1,3,4,5,6,8]:
    tmp_supports = np.asarray(support_values[category])
    hell = []
    for bin in set(binning):
        binned_support = tmp_supports[binning==bin]
        hell.append(binned_support)
        tmp_df = pd.DataFrame(zip([bin]*binned_support.shape[0], [category]*binned_support.shape[0], binned_support),
                              columns='bin category support'.split())
        yeah = yeah.append(tmp_df)
    sns.kdeplot(hell, ax=category - 1)

fig, ax = plt.subplots()
sns.boxplot(x='bin', y='support', hue='category', data=yeah, ax=ax)
fig.set_size_inches(15,6)
fig.tight_layout()
fig.savefig('test1.pdf', dpi=300)
plt.close()

for category, values in support_values.items():
    chart = sns.jointplot(pd.Series(data=branch_lengths, name='branch length'),
                          pd.Series(data=support_values[category], name='Site-rate category %i support' %category))
    chart.fig.set_size_inches(10,10)
    chart.fig.tight_layout()
    chart.fig.savefig('test-%i.png' %category, dpi=300)
    plt.close()

support_values = {}
for category in range(1,9):
    support_values[category] = {'distance_to_root':[], 'bootstrap':[]}

for node in full_tree.traverse():
    if node.is_leaf():
        continue

    distance_to_root = full_tree.get_distance(node, topology_only=True)

    for category in range(1,9):
        support_values[category]['distance_to_root'].append(distance_to_root)
        support_values[category]['bootstrap'].append(eval('node.category_support_%i' %category))

cNorm      = colors.Normalize(vmin=0, vmax=8)
scalarMap  = cmx.ScalarMappable(norm=cNorm, cmap=plt.get_cmap('Set1'))

fig, ax = plt.subplots()
for num in range(18):
    ax.add_patch(plt.Rectangle((num-0.4, -5), 0.8, 110, fill=True, color='grey', alpha=0.1))

#for category in [1,2,4,6,8]:
for category in range(1,9):
    tmp_df             = pd.DataFrame.from_dict(support_values[category])
    tmp_df['category'] = category

    color = '#%02x%02x%02x' %scalarMap.to_rgba(category, bytes=True)[:3]

    sns.stripplot(x='distance_to_root', y='bootstrap', data=tmp_df, color=color, edgecolor='white', linewidth=1, size=10, jitter=True, alpha=0.6, ax=ax)
    sns.lineplot( x='distance_to_root', y='bootstrap', data=tmp_df, color=color, ci=None, ax=ax, label=str(category))

fig.set_size_inches(15,8)
ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
fig.tight_layout()
fig.savefig('yeah.pdf', dpi=300)
plt.close()

########################################################################################################################
os.chdir('/work/site_rate/fixed_tree/simulated_alignmentsG')

support_values = {}
for category in range(1,9):
    tree = ete3.Tree('%i.aln.treefile' %category, format=1)
    if tree.check_monophyly(assembly_summary_outgroup.index, 'name'):
        outgroup = tree.get_common_ancestor(assembly_summary_outgroup.index.tolist())
        tree.set_outgroup(outgroup)
    else:
        print 'not monophyletic outgroup'
        continue

    support_values[category] = {'distance_to_root':[], 'alrt':[], 'branch_length':[]}

    for node in tree.traverse():
        if node.is_leaf() or not node.name:
            continue
        original_alrt, original_bb, alrt = [float(support) for support in node.name.split('/')]
        support_values[category]['alrt'            ].append(alrt)
        support_values[category]['branch_length'   ].append(node.dist)
        support_values[category]['distance_to_root'].append(tree.get_distance(node, topology_only=True))

colors = 'red blue green magenta black'.split()
fig, ax = plt.subplots()
for category, color in zip([1,4,8], colors):
    tmp_df             = pd.DataFrame.from_dict(support_values[category])
    tmp_df['category'] = category

    sns.stripplot(x='distance_to_root', y='alrt', data=tmp_df, color=color, edgecolor='white', linewidth=1, size=10, jitter=True, alpha=0.6, ax=ax)
    sns.lineplot( x='distance_to_root', y='alrt', data=tmp_df, color=color, ci=None, ax=ax, label=str(category))
    #sns.lineplot( x='distance_to_root', y='alrt', style='category', markers=True, dashes=True, data=tmp_df, color=color, ci=None, ax=ax, label=str(category))

fig.set_size_inches(15,8)
ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
fig.tight_layout()
fig.savefig('../yeah2.pdf', dpi=300)
plt.close()

























