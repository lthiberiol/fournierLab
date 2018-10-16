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

num_replicates      = 100
num_rate_categories = 12
trees          = {partition_name:ete3.Tree('%s/reference.tre' % partition_name, format=5) for partition_name in 'base_tree short2long long2short'.split()}

support_values = {}
for partition_name, tree in trees.items():
    outgroup_taxa                  = sorted(tree.children, key=len)[0].get_leaf_names()
    support_values[partition_name] = {category:[] for category in range(1, num_rate_categories+1)}

    for node in tree.traverse():
        if node.is_leaf():
            continue

        node.add_feature('topology_id', node.get_topology_id())

        for category in range(1, num_rate_categories+1):
            node.add_feature('category_support_%i' % category, [])

    with cd('%s/categories' % partition_name):
        for replicate in range(1, num_replicates+1):
            tmp_trees = [ete3.Tree('%i.%i.treefile' % (replicate, category)) for category in range(1, num_rate_categories+1)]
            [tmp_tree.set_outgroup(tmp_tree.get_common_ancestor(outgroup_taxa)) for tmp_tree in tmp_trees]
            [tmp_tree.ladderize()                                               for tmp_tree in tmp_trees]

            traversing_trees = [tmp_tree.traverse() for tmp_tree in tmp_trees]

            for replicated_nodes in zip(*traversing_trees):
                if replicated_nodes[0].is_leaf():
                    continue

                topology_ids = set([replicated_node.get_topology_id() for replicated_node in replicated_nodes])
                if len(topology_ids) != 1:
                    print 'SOMETHING IS FUCKING WRONG!!!!!'
                    break

                reference_topology_id = topology_ids.pop()
                reference_node        = tree.search_nodes(topology_id=reference_topology_id)[0]

                for category, replicated_node in enumerate(replicated_nodes):
                    eval('reference_node.category_support_%i.append(replicated_node.support)' % (category+1))

for partition_name, tree in trees.items():
    sampled_branch_lengths = np.asarray([node.dist for node in tree.traverse() if not node.is_leaf() and not node.is_root()])

    branch_length_bins     = [np.percentile(sampled_branch_lengths, decile) for decile in range(10, 91, 10)]
#    branch_length_bins     = np.linspace(sampled_branch_lengths.min(),
#                                         sampled_branch_lengths.max(),
#                                         10)
    binning                = np.digitize(sampled_branch_lengths, branch_length_bins)

    support_df             = pd.DataFrame(columns='branch length bin\tcategory\tsupport'.split('\t'))

    for category in range(1,num_rate_categories+1):
        tmp_supports = np.asarray([eval('node.category_support_%i' %category) for node in tree.traverse() if not node.is_leaf() and not node.is_root()])

        for bin in set(binning):
            binned_support = tmp_supports[binning==bin, :].flatten()
            min_binned_branch_len = sampled_branch_lengths[binning==bin].min()
            max_binned_branch_len = sampled_branch_lengths[binning==bin].max()

            tmp_df = pd.DataFrame(
                zip(['%.4f - %.4f' % (min_binned_branch_len, max_binned_branch_len)]*binned_support.shape[0],
                    [category]*binned_support.shape[0],
                    binned_support
                    ),
                columns='branch length bin\tcategory\tsupport'.split('\t'))
            support_df = support_df.append(tmp_df)

    fig, ax = plt.subplots()
    sns.boxplot(x='branch length bin', y='support', hue='category', data=support_df, ax=ax)
    ax.set_title(partition_name)
    fig.set_size_inches(25,6)
    ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1.015), title='Site-rate category',frameon=False)
    fig.tight_layout()
    fig.savefig('support_binned_by_branch_length-decile_spaced-%s.png' %partition_name, dpi=300)
    plt.close()
