import ete3
import os
import re
import subprocess
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import pickle as pkl
from itertools import combinations
import multiprocessing
import pandas as pd

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

########################################################################################################################
#                                                                                                                      #
# ASSESS BALANCE                                                                                                       #
########################################################################################################################
def assess_tree_balance_ranger(folder):
    ratios = []
    trees = []
    count = 1
    while os.path.isfile('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)):
        trees.append(ete3.Tree(open('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)).readlines()[1]))
        count += 1

    for tree in trees:
        children = sorted(tree.children, key=len)
        ratios.append(len(children[0])/float(len(children[1])))

    return ratios

def assess_tree_balance_mad(folder):
    if not os.path.isfile('%s/%s-MAD.ranger_input' %(folder, folder)):
        return []

    tree = ete3.Tree(open('%s/%s-MAD.ranger_input' %(folder, folder)).readlines()[1])
    children = sorted(tree.children, key=len)

    return len(children[0])/float(len(children[1]))

pool = multiprocessing.Pool(processes=6)
tmp_balance = pool.map(assess_tree_balance_mad, os.listdir('.'))
mad_balance = []
for n in tmp_balance:
    if n:
        mad_balance.append(n)

pool = multiprocessing.Pool(processes=8)
tmp_balance = pool.map(assess_tree_balance_ranger, os.listdir('.'))
ranger_balance = []
for n in tmp_balance:
    ranger_balance.extend(n)

fig, axs = plt.subplots(nrows=2, sharex=True)
sns.kdeplot(ranger_balance, shade=True, color='blue', ax=axs[0])
axs[0].set_title('RANGER topologies tree balance')
axs[0].set_ylabel('Frequency')
sns.kdeplot(mad_balance, shade=True, color='green', ax=axs[1])
axs[1].set_title('MAD topologies tree balance')
axs[1].set_xlabel('root children balance distribution')
axs[1].set_ylabel('Frequency')
fig.tight_layout()
fig.savefig('balance_distribution.pdf', dpi=300)

########################################################################################################################
#                                                                                                                      #
# ASSESS ALPHA AND CYANO SEPARATION                                                                                    #
########################################################################################################################
os.chdir('/work/Alphas_and_Cyanos/reconciliations')

reference_tree = ete3.Tree('/work/Alphas_and_Cyanos/rooted_partitions-with_named_branches.treefile', format=1)

ncbi = ete3.NCBITaxa()
header = 'assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'.split()
genbank_summary                     = pd.read_table('/work/assembly_summary_genbank.txt', comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str})
genbank_summary['refseq_category']  = genbank_summary['refseq_category'].str.lower()
genbank_summary['assembly_level']   = genbank_summary['assembly_level'].str.lower()
genbank_summary['genome_rep']       = genbank_summary['genome_rep'].str.lower()
genbank_summary.set_index('assembly_accession', inplace=True)
genbank_summary.index               = [re.sub('\.\d+$', '', index).replace('_', '') for index in genbank_summary.index]
genbank_summary                     = genbank_summary.reindex(labels=reference_tree.get_leaf_names(), axis='index')

lineages = pd.DataFrame.from_dict({taxid:{j: i for i, j in ncbi.get_rank(ncbi.get_lineage(taxid)).items()} for taxid in genbank_summary.taxid}).T

def assess_separation(folder):
    if not os.path.isfile('%s/%s-MAD.ranger_input' %(folder, folder)):
        return 0

    tree   = ete3.Tree(open('%s/%s-MAD.ranger_input' %(folder, folder)).readlines()[1])
    taxids = genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.get_leaf_names()]), 'taxid']
    phyla  = lineages.loc[taxids, 'phylum'].values.astype('int').astype('str').tolist()

    if not len(set(phyla)) == 2 or not set(['1117', '1224']).issubset(phyla):
        return 0

    child1_phyla = set(lineages.loc[genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.children[0].get_leaf_names()]), 'taxid'],
                                    'phylum'].values.astype('int').astype('str').tolist())
    child2_phyla = set(lineages.loc[genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.children[1].get_leaf_names()]), 'taxid'],
                                    'phylum'].values.astype('int').astype('str').tolist())

    if len(child1_phyla) == len(child2_phyla) == 1:
        return 1
    else:
        return 0

pool = multiprocessing.Pool(processes=8)
tmp_separation = pool.map(assess_separation, os.listdir('.'))

def assess_separation_ranger(folder):
    trees = []
    count = 1
    while os.path.isfile('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)):
        trees.append(ete3.Tree(open('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)).readlines()[1]))
        count += 1

    if not trees:
        return 0

    for tree in trees:
        taxids = genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.get_leaf_names()]), 'taxid']
        phyla  = lineages.loc[taxids, 'phylum'].values.astype('int').astype('str').tolist()

        if not len(set(phyla)) == 2 or not set(['1117', '1224']).issubset(phyla):
            continue

        child1_phyla = set(lineages.loc[genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.children[0].get_leaf_names()]), 'taxid'],
                                    'phylum'].values.astype('int').astype('str').tolist())
        child2_phyla = set(lineages.loc[genbank_summary.loc[set([leaf.split('_')[0] for leaf in tree.children[1].get_leaf_names()]), 'taxid'],
                                    'phylum'].values.astype('int').astype('str').tolist())

        if len(child1_phyla) == len(child2_phyla) == 1:
            return 1

    return 0

pool = multiprocessing.Pool(processes=8)
ranger_separation = pool.map(assess_separation_ranger, os.listdir('.'))
