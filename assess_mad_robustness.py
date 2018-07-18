import ete3
import os
import re
import subprocess
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

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

os.chdir('/work/Alphas_and_Cyanos/test_mad_rooting_consistency')

def compare_root_position(subtree, newly_rooted):
    for child in subtree.children:
        is_it_monophyletic, clade_type, fucking_up = newly_rooted.check_monophyly(child.get_leaf_names(), 'name',
                                                                                  unrooted=False)
        if is_it_monophyletic or child.is_leaf():
            break

    if not is_it_monophyletic and not child.is_leaf():
        return None

    if child.is_leaf():
        new_root = newly_rooted.get_leaves_by_name(name=child.name)[0]
    else:
        new_root = newly_rooted.get_common_ancestor(child.get_leaf_names())
    return newly_rooted.get_distance(new_root, topology_only=True)


def assess_mad_consistency(tree_file):
    group     = tree_file.replace('.tree.rooted', '')
    full_tree = ete3.Tree('../ranger_input_trees-no_long_branches/%s' %tree_file)
    if os.path.isdir(group):
        os.system('rm -rf %s' %group)
    os.mkdir(group)
    root_distances = []
    with cd(group):
        for count, child in enumerate(full_tree.children):
            txt = child.write(format=5)
            out = open('child_%i.tree' %count, 'wb')
            out.write(re.sub('\):\d+\.\d+;$', ');', txt, flags=re.M))
            out.close()

            subprocess.call(['/Users/thiberio/anaconda2/envs/py37/bin/python', '/work/mad.py', 'child_%i.tree' %count])

            if not os.path.getsize('child_%i.tree.rooted' %count):
                continue

            minimal_value = 1000000000
            for line in open('child_%i.tree.rooted' %count).readlines():
                if not line.strip():
                    continue
                newly_rooted = ete3.Tree(line)
                root_distance = compare_root_position(child, newly_rooted)
                if root_distance < minimal_value:
                    minimal_value = root_distance
            root_distances.append(minimal_value)

    return root_distances

os.chdir('/work/Alphas_and_Cyanos')

########################################################################################################################
less_stringent = []
for tmp in pkl.load(open('ranger_input_trees-no_long_branches/subtrees_root_distances.pkl')):
    less_stringent.extend(tmp)

more_stringent = []
similar = []
dissimilar = []
for tmp in pkl.load(open('ranger_input_trees-no_long_branches2/subtrees_root_distances.pkl')):
    more_stringent.extend(tmp)
    close_flag = True
    for n in tmp:
        if n > 2:
            close_flag = False
            break

    if close_flag:
        similar.append(tmp)
    else:
        dissimilar.append(tmp)

fig, ax = plt.subplots()
sns.kdeplot(less_stringent, ax=ax, label='Less stringent branch length cut-off')
sns.kdeplot(more_stringent, ax=ax, label='More stringent branch length cut-off')
ax.set_xlabel('Distance between root positions')
ax.set_ylabel('Frequency')
fig.tight_layout()
fig.savefig('yeah.pdf')
