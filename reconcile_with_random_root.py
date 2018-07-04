#/usr/bin/env python
#coding: utf-8

import os
import multiprocessing
import ete3
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

def root_randomly(tree):
    random_node = random.choice(list(tree.traverse()))
    tree.set_outgroup(random_node)
    return tree.write(format=9)

def run_ranger(tree_file):
    os.mkdir('random_root_reconciliations/%s' %tree_file)
    unrooted_tree = ete3.Tree('ranger_input_trees/%s.tree' %tree_file)
    rooted_tree   = root_randomly(unrooted_tree)
    os.system('cat rooted_partitions-with_BB_support.treefile > random_root_reconciliations/%s/%s-random_root.ranger_input' %(tree_file, tree_file))
    with cd('random_root_reconciliations/%s' %tree_file):
        out = open('%s-random_root.ranger_input' %tree_file, 'a')
        out.write(rooted_tree)
        out.close()
        for count in range(20):
            os.system('/home/thiberio/ranger/CorePrograms/Ranger-DTL.linux -i {tree_file}-random_root.ranger_input -o {tree_file}-random_root.ranger_out{count}'.format(tree_file=tree_file, count=count+1))

tree_files = open('groups_with_single_optimal_DTL_root.list').read().split()

pool = multiprocessing.Pool(processes=20)
pool.map_async(run_mad, tree_files)
pool.close()
pool.join()
#run_mad(tree_files[-1])
