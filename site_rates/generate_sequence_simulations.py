import subprocess
import os
import ete3
import random
import numpy as np

#
# initial definitions
os.chdir('/work/site_rate/sequence_simulation')
random.seed(5)

indelible_conf = '''\
/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
//  INDELible V1.03 control file - site-rate project                               //
//                                                                                 //
//      Automaticaly generated, more information:                                  //
//          github.com/lthiberiol/fournierLab/tree/master/site_rates               //
//                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////

[TYPE]  AMINOACID 1

[MODEL] model1  [submodel]  LG         //   LG
                [rates]     0 1.3 8    //   pinv=0.1, alpha=1.3, discrete gamma categories=8

[TREE] T1  {t1}
[TREE] T2  {t2}
[TREE] T3  {t3}

[PARTITIONS] partition1 [T1 model1 1000]   // tree T1, model model1, root length 1000
[PARTITIONS] partition2 [T2 model1 1000]   // tree T2, model model1, root length 1000
[PARTITIONS] partition3 [T3 model1 1000]   // tree T3, model model1, root length 1000

[EVOLVE]    partition1  10   {t1_name}
            partition2  10   {t2_name}
            partition3  10   {t3_name}
'''

#
# generate base tree, if necessary...
#base_tree = ete3.Tree()
#base_tree.populate(10)
#for leaf in base_tree.get_leaves():
#    leaf.name = leaf.name.replace('aaaaaaaaa', '')
base_tree = ete3.Tree('(((c,(d,e)),((f,(g,h)),(i,j))),(a,b));')
for node in base_tree.traverse():
    if node.is_root():
        continue
    node.dist = random.random()

branch_lengths = np.geomspace(0.01, 1, int(base_tree.get_farthest_leaf(topology_only=True)[1]+1))
short2long = base_tree.copy()
for node in short2long.traverse():
    if node.is_root():
        continue
    distance_from_root = int(short2long.get_distance(node, topology_only=True))
    node.dist          = branch_lengths[distance_from_root]

branch_lengths = sorted(branch_lengths, reverse=True)
long2short = base_tree.copy()
for node in long2short.traverse():
    if node.is_root():
        continue
    distance_from_root = int(long2short.get_distance(node, topology_only=True))
    node.dist          = branch_lengths[distance_from_root]

out = open('control.txt', 'w')
out.write(indelible_conf.format(t1=base_tree.write( format=5), t1_name='random_branch_length',
                                t2=short2long.write(format=5), t2_name='short_to_long_branches',
                                t3=long2short.write(format=5), t3_name='long_to_short_branches'))
out.close()

subprocess.call(['/work/site_rate/indelible/INDELibleV1.03/bin/indelible_1.03_OSX_intel'])

for partition_name in 'random_branch_length short_to_long_branches long_to_short_branches'.split():
    fasta = open('%s.fas' %partition_name).read().strip()
    for count, block in enumerate(fasta.split('\n     \n')):
        out = open('%s/%i.fas' %(partition_name, count+1), 'w')
        out.write(block)
        out.close()
