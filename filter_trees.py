import ete3
import os
import numpy as np
import pandas as pd
import re
from matplotlib import pyplot as plt
import seaborn as sns
import pickle as pkl

os.chdir('/work/Alphas_and_Cyanos')

mad_stats = pd.DataFrame(columns='deviation ratio'.split())
multiple_roots = []
for result in os.listdir('ranger_input_trees-pruned_long_branches'):
    if not result.endswith('.tree.rooted'):
        continue
    tmp_stats = re.findall('^>> \[MAD=(\d+\.\d+)_AI=(\d+\.\d+)_CCV=\S+?%_N=\d/\d\]$', open('ranger_input_trees-pruned_long_branches/%s' %result).read(), re.M)
    group     = result.replace('.tree.rooted', '')
    if len(tmp_stats) > 1:
        multiple_roots.append(group)
    else:
        mad_stats.loc[result.replace('.tree.rooted', '')] = [float(n) for n in tmp_stats[0]]

fig, axs = plt.subplots(nrows=2)
sns.distplot(mad_stats.deviation, ax=axs[0])
sns.distplot(mad_stats.ratio,     ax=axs[1])
axs[0].set_title('Minimal ancestral deviation')
axs[1].set_title('Ratio with 2nd best root position')
fig.tight_layout()
fig.savefig('mad_rooting_stats.pdf', dpi=300)

tree_support_desc = pd.DataFrame(columns='mean median 1st_quartile 3rd_quartile'.split())
for result in os.listdir('ranger_input_trees-pruned_long_branches'):
    if not result.endswith('.tree'):
        continue

    group = result.replace('.tree', '')
    if group in multiple_roots:
        continue

    tmp_tree = ete3.Tree('ranger_input_trees-pruned_long_branches/%s' %result)
    branch_support = np.asarray([node.support for node in tmp_tree.traverse() if not node.is_leaf()])
    tree_support_desc.loc[group] = [branch_support.mean(),
                                    np.median(branch_support),
                                    np.percentile(branch_support, 25),
                                    np.percentile(branch_support, 75)]

filtered_groups = tree_support_desc.index[tree_support_desc['1st_quartile'] >= 80]
out = open('families_filtered_by_multiRoots_and_support.pkl', 'wb')
pkl.dump(filtered_groups.tolist(), out)
out.close()