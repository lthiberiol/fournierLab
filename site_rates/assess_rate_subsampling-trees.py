#!/usr/bin/env python2
#coding: utf-8

import ete3
import os
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import subprocess
from scipy.stats import spearmanr, linregress
import itertools
from Bio import SeqIO, SearchIO, AlignIO, Align, Alphabet
import re
import multiprocessing
from copy import deepcopy

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
association_data   = {category:dict(x=[], y=[]) for category in range(1,9)}
correlation_values = {category:[] for category in range(1,9)}
slope_values       = {category:[] for category in range(1,9)}
for leaf in full_tree.get_leaves():
    if len(leaf.get_ancestors()) < 10:
        continue

    for category in correlation_values.keys():
        tmp_supports                 = []
        tmp_num_bipartitions_to_root = []
        for ancestor in leaf.get_ancestors()[:-1]:
            if not ancestor.name:
                continue
            tmp_alrt, tmp_bb = [float(value) for value in ancestor.name.split('/')]
            if tmp_bb < 80 or tmp_alrt < 80:
                continue
            tmp_supports.append(eval('ancestor.category_support_%i' %category))
            tmp_num_bipartitions_to_root.append(full_tree.get_distance(ancestor, topology_only=True))

        if len(tmp_supports) < 7:
            continue
        association_data[category]['x'].append(tmp_num_bipartitions_to_root)
        association_data[category]['y'].append(tmp_supports)

        tmp_spearman, tmp_p = spearmanr(tmp_supports, tmp_num_bipartitions_to_root)
        tmp_regression      = linregress(tmp_supports, tmp_num_bipartitions_to_root)

        correlation_values[category].append(tmp_spearman)
        slope_values[      category].append(tmp_regression.slope)
#        print '%i: %.2f' %(category, tmp_regression.rvalue)

#    print ''

fig, axs = plt.subplots(nrows=8, sharex=True)
for category, values in correlation_values.items():
    axs[category-1].set_title('Site-rate category %i support correlation' %category)
    sns.kdeplot(values, shade=True, ax=axs[category-1])
fig.set_size_inches(15,12)
fig.tight_layout()
fig.savefig('site_rate_support_correlations.pdf', dpi=300)
plt.close()

fig, axs = plt.subplots(nrows=8, sharex=True)
for category, values in slope_values.items():
    axs[category-1].set_title('Site-rate category %i support regression slope distribution' %category)
    sns.kdeplot(values, shade=True, ax=axs[category-1])
fig.set_size_inches(15,12)
fig.tight_layout()
fig.savefig('site_rate_support_regression_slopes.pdf', dpi=300)
plt.close()

fig, ax = plt.subplots()
for x, y in zip(association_data[1]['x'],association_data[1]['y']):
    sns.regplot(x, y, robust=True, truncate=True, ci=False, color='grey', scatter_kws=dict(alpha=0.0), line_kws=dict(alpha=0.6), ax=ax)
fig.tight_layout()
fig.savefig('yeah.pdf', dpi=300)
plt.close()

#
# assess correlations for branch length
support_values = {category:[] for category in range(1,9)}
branch_lengths = []
for node in full_tree.traverse():
    if node.is_leaf():
        continue

    branch_lengths.append(node.dist)
    for category in correlation_values.keys():
        support_values[category].append(eval('node.category_support_%i' %category))

branch_length_bins = [np.percentile(branch_lengths, decile) for decile in range(20, 81, 20)]
binning            = np.digitize(branch_lengths, branch_length_bins)

support_df        = pd.DataFrame(columns='branch length bin\tcategory\tsupport'.split('\t'))
branch_lengths    = np.asarray(branch_lengths)
#for category in [1,3,5,8]:
for category in range(1,9):
    tmp_supports = np.asarray(support_values[category])
    for bin in set(binning):
        binned_support = tmp_supports[binning==bin]
        tmp_df = pd.DataFrame(zip([bin+1]*binned_support.shape[0], [category]*binned_support.shape[0], binned_support),
                              columns='branch length bin\tcategory\tsupport'.split('\t'))
        support_df = support_df.append(tmp_df)

fig, ax = plt.subplots()
sns.boxplot(x='branch length bin', y='support', hue='category', data=support_df, ax=ax)
fig.set_size_inches(15,6)
ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.015), title='Site-rate category',frameon=False)
fig.tight_layout()
fig.savefig('support_binned_by_branch_length.pdf', dpi=300)
plt.close()
