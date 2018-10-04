import os
import ete3
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import re
import random


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
os.chdir('/work/site_rate/sequence_simulation/')

num_replicates = 100
trees          = {partition_name:ete3.Tree('%s/reference.tre' % partition_name, format=5) for partition_name in 'base_tree short2long long2short'.split()}
branch_lengths = np.geomspace(0.01, 1, int(trees['base_tree'].get_farthest_leaf(topology_only=True)[1]+1))

support_values = {}
for partition_name, tree in trees.items():
    support_values[partition_name] = {category:[] for category in range(1, 9)}

    for node in tree.traverse():
        if node.is_leaf():
            continue

        for category in range(1,9):
            node.add_feature('category_support_%i' % category, [])

    with cd('%s/categories' % partition_name):
        for replicate in range(1, num_replicates+1):
            tmp_trees = [ete3.Tree('RAxML_bipartitions.%i.%i' % (replicate, category)) for category in range(1,9)]
#            tree_txts = [open('%i.%i.boottrees.suptree' % (replicate, category)).readline().replace('((a', '(a') for category in range(1,9)]
#            tmp_trees = [ete3.Tree(re.sub('\d+:1.0000000000\);', ';', tmp_txt)) for tmp_txt in tree_txts]
#
#            [tmp_tree.set_outgroup(tmp_tree.get_common_ancestor('a b'.split())) for tmp_tree in tmp_trees]
#            [tmp_tree.ladderize                                                 for tmp_tree in tmp_trees]
#
#            traversing_trees = [tmp_tree.traverse() for tmp_tree in tmp_trees]

            for node, replicated_nodes in zip(trees[partition_name].traverse(), zip(*tmp_trees)):
                if node.is_leaf():
                    continue

                ########################################################################################################
                #
                # reality check if the tree and replicates are traversed equally!!!
                topology_ids = set([replicated_node.get_topology_id() for replicated_node in replicated_nodes])
                topology_ids.add(node.get_topology_id())

                if len(topology_ids) != 1:
                    print 'SOMETHING IS FUCKING WRONG!!!!!'
                    break
                ########################################################################################################

                for category, replicated_node in enumerate(replicated_nodes):
                    eval('node.category_support_%i.append(replicated_node.support)' % (category+1))

for partition_name, tree in trees.items():
    sampled_branch_lengths = [node.dist for node in tree.traverse() if not node.is_leaf()]

    branch_length_bins     = [np.percentile(sampled_branch_lengths, decile) for decile in range(20, 81, 20)]
    binning                = np.digitize(sampled_branch_lengths, branch_length_bins)

    support_df             = pd.DataFrame(columns='branch length bin\tcategory\tsupport'.split('\t'))
    sampled_branch_lengths = np.asarray(sampled_branch_lengths)

    for category in range(1,9):
        tmp_supports = np.asarray([eval('node.category_support_%i' %category) for node in tree.traverse() if not node.is_leaf()])

        for bin in set(binning):
            binned_support = tmp_supports[binning==bin, :].flatten()
            tmp_df = pd.DataFrame(zip([bin+1]*binned_support.shape[0], [category]*binned_support.shape[0], binned_support),
                                  columns='branch length bin\tcategory\tsupport'.split('\t'))
            support_df = support_df.append(tmp_df)

    fig, ax = plt.subplots()
    sns.boxplot(x='branch length bin', y='support', hue='category', data=support_df, ax=ax)
    ax.set_title(partition_name)
    fig.set_size_inches(15,6)
    ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.015), title='Site-rate category',frameon=False)
    fig.tight_layout()
    fig.savefig('support_binned_by_branch_length-%s.pdf' %partition_name, dpi=300)
    plt.close()

    break
