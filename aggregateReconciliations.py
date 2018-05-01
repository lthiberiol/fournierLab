#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
import itertools
from sys import argv, exit
import pickle as pkl
import pandas as pd

output_folder = '/work/Alphas_and_Cyanos/aggregated'
folder = argv[1]
group  =  folder.strip('/').split('/')[-1]
if not os.path.isdir(folder) or not os.path.isfile('%s/%s.reconciliation1' %(folder, group)):
    exit('\n\t**ERROR, no valid reconciliation!\n')

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

def get_branches_with_transfers(reconciliation_file, leaves_allowed=False):
    reconciliation = open(reconciliation_file).read()
    gene_tree = ete3.Tree(reconciliation.split('\n')[7], format=1)
    if leaves_allowed:
        for leaf1, leaf2, donor, recipient in re.findall('^m\d+\s=\sLCA\[(\S+),\s(\S+)\]:\sTransfer,\sMapping\s-->\s(\S+),\sRecipient\s-->\s(\S+)$', reconciliation, re.M):
            reticulation = gene_tree.get_common_ancestor(leaf1, leaf2)
            yield {'reticulation':reticulation, 'topology_id':reticulation.get_topology_id(), 'donor':donor, 'recipient':recipient}
    else:
        for leaf1, leaf2, donor, recipient in re.findall('^m\d+\s=\sLCA\[(\S+),\s(\S+)\]:\sTransfer,\sMapping\s-->\s(n\d+),\sRecipient\s-->\s(n\d+)$', reconciliation, re.M):
            reticulation = gene_tree.get_common_ancestor(leaf1, leaf2)
            yield {'reticulation':reticulation.copy(method='deepcopy'), 'topology_id':reticulation.get_topology_id(), 'donor':donor, 'recipient':recipient}

def get_shared_transfers(transfer_descriptions, num_replicates, threshold=0.9):
    topology_ids        = [[transfer['topology_id'] for transfer in transfers] for transfers in transfer_descriptions]
    shared_topology_ids = set.intersection(*map(set, topology_ids))

    #
    # if no reticulation is shared by all reconciliations, ignore it
    if not shared_topology_ids:
        return None

    shared_reticulations = [tmp for tmp in itertools.chain.from_iterable(transfer_descriptions) if tmp['topology_id'] in shared_topology_ids]
    shared_reticulations = pd.DataFrame(shared_reticulations)

    equivalent_transfers = shared_reticulations.groupby(by='topology_id donor recipient'.split(), sort=False)
    shared_transfers     = [equivalent_transfers.get_group(keys).iloc[0].values.tolist() for keys, tmp_df in equivalent_transfers.groups.items() if tmp_df.shape[0] >= num_replicates*threshold]

    return shared_transfers

def traverse_reconciliations(folder):
    transfer_descriptions = []
    counter              = 0
    while True:
        counter += 1
        if not os.path.isfile('%s.reconciliation%i' %(folder, counter)):
            break
        transfer_descriptions.append(list(get_branches_with_transfers('%s.reconciliation%i' %(folder, counter))))

    shared_transfers = get_shared_transfers(transfer_descriptions, num_replicates=counter-1, threshold=0.9)

    return shared_transfers

with cd(folder):
    result = traverse_reconciliations(group)

out = open('%s/%s.pkl' %(output_folder, group), 'wb')
pkl.dump(result, out)
out.close()