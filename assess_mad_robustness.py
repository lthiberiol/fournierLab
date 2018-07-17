import ete3
import os
import re
import subprocess

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
