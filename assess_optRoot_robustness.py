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

os.chdir('/work/Alphas_and_Cyanos/test_ranger_rooting_consistency')

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

def assess_mad_consistency(group):
    root_distances = []
    for tree_file in os.listdir('../reconciliations/ranger_roots/%s' %group):
        if not tree_file.endswith('.ranger_input'):
            continue

        print tree_file
        full_tree = ete3.Tree(open('../reconciliations/ranger_roots/%s/%s' %(group, tree_file)).readlines()[1])
        if os.path.isdir(group):
            os.system('rm -rf %s' %group)
        os.mkdir(group)
        with cd(group):
            for count, child in enumerate(full_tree.children):
                if child.is_leaf() or len(child) <= 3:
                    continue

                subprocess.call(['cp', '../../rooted_partitions-with_BB_support.treefile', 'child_%i.tree' %count])
                txt = child.write(format=5)
                out = open('child_%i.tree' %count, 'a')
                out.write(re.sub('\):\d+\.\d+;$', ');', txt, flags=re.M))
                out.close()

                subprocess.call(['/work/ranger/CorePrograms/OptRoot.mac', '-i', 'child_%i.tree' %count, '-o', 'child_%i.tree.rooted' %count])

                if not os.path.getsize('child_%i.tree.rooted' %count):
                    continue

                minimal_value = 1000000000
                for line in open('child_%i.tree.rooted' %count).readlines():
                    line = line.strip()
                    if line.startswith('(') and line.endswith(';'):
                        newly_rooted = ete3.Tree(line)
                        root_distance = compare_root_position(child, newly_rooted)
                        if root_distance < minimal_value:
                            minimal_value = root_distance
                    root_distances.append(minimal_value)

    return root_distances

os.chdir('/work/Alphas_and_Cyanos')

########################################################################################################################
def distances_between_roots(folder):
    trees = []
    count = 1
    while os.path.isfile('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)):
        trees.append(ete3.Tree(open('%s/%s.optResolution%i.ranger_input' %(folder, folder, count)).readlines()[1]))
        count += 1

    root_distances = []
    for tree1, tree2 in combinations(trees, 2):
        for child in tree2.children:
            if child.is_leaf():
                new_root = tree1.get_leaves_by_name(name=child.name)[0]
                break
            else:
                is_it_monophyletic, clade_type, fucking_up = tree1.check_monophyly(child.get_leaf_names(), 'name',
                                                                                          unrooted=False)
                if is_it_monophyletic or child.is_leaf():
                    new_root = tree1.get_common_ancestor(child.get_leaf_names())
                    break
        root_distances.append(tree1.get_distance(new_root, topology_only=True))

    return root_distances

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

root_distances = []
for tmp in pkl.load(open('subtrees_root_distances.pkl')):
    root_distances.extend(tmp)

root_distances = np.asarray(root_distances)
root_distances = root_distances[root_distances < 1000000000]

fig, ax = plt.subplots()
sns.kdeplot(root_distances, ax=ax)
ax.set_xlabel('Distance between root positions')
ax.set_ylabel('Frequency')
fig.tight_layout()
fig.savefig('yeah.pdf')
