#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
from commands import getoutput
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import pearsonr
import pandas as pd
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

def name_matching_branches( named, unamed ):
    if named.get_topology_id() == unamed.get_topology_id():
        unamed.name = named.name

    for node1 in unamed.children:
        if node1.is_leaf():
            continue

        for node2 in named.children:
            if node2.is_leaf():
                continue

            if node1.get_topology_id() == node2.get_topology_id():
                node1.name = node2.name
                name_matching_branches( node2, node1 )

os.chdir('/work/Alphas_and_Cyanos')
named_reference_tree = ete3.Tree('rooted_partitions-with_named_branches.treefile', format=1)

def assess_branch_compatibility(folder, transfers, gene_tree, named_reference_tree):
    final_transfers = []
    rf_values = {'donor':[], 'recipient':[], 'donor_normed':[], 'recipient_normed':[]}
    for transfer in transfers:
        grep_query  = re.match('^(m\d+ = LCA\[\S+, \S+\]:)', transfer, re.M).group(1)
        grep_result = getoutput( 'grep "%s" %s.reconciliation*' %(re.escape(grep_query), folder))
        grep_result = set(re.sub('^%s.reconciliation\d+:' %folder, '', grep_result, flags=re.M).split('\n'))
        #
        # assess consistency of donor/recipient pair
        if len(grep_result) > 1:
            #
            # inconsistent... sorry
            continue

        transfer         = grep_result.pop()
        donor, recipient = re.search('Mapping --> (\S+), Recipient --> (\S+)$', transfer, re.M).groups()
        recipient_branch = named_reference_tree.search_nodes(name=recipient)[0].copy(method='deepcopy')
        donor_branch     = named_reference_tree.search_nodes(name=donor    )[0].copy(method='deepcopy')
        if recipient_branch.is_leaf() or donor_branch.is_leaf():
            continue

        reticulation = gene_tree.search_nodes(name=transfer.split()[0])[0]

        #
        # assess donor branch compatibility
        rf, rf_max, names, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = reticulation.robinson_foulds(donor_branch, attr_t1='genome_name')
        if rf:
            return None

        #
        # assess recipient branch compatibility
        rf, rf_max, names, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = reticulation.robinson_foulds(recipient_branch, attr_t1='genome_name')
        if rf:
            return None

        donor_intersection     = [leaf.name for leaf in reticulation.get_leaves() if leaf.genome_name in donor_branch.get_leaf_names()    ]
        recipient_intersection = [leaf.name for leaf in reticulation.get_leaves() if leaf.genome_name in recipient_branch.get_leaf_names()]

        if not donor_intersection or not recipient_intersection:
            continue

        donor_branch_in_reticulation     = reticulation.get_common_ancestor(donor_intersection    )
        recipient_branch_in_reticulation = reticulation.get_common_ancestor(recipient_intersection)

        if donor_branch_in_reticulation.name == 'm1' or recipient_branch_in_reticulation.name == 'm1':
            continue

        if donor_branch_in_reticulation.name in [node.name for node in recipient_branch_in_reticulation.get_ancestors()]:
            print 'yeah'
        else:
            continue
        print reticulation.name
        print reticulation
        print donor_branch_in_reticulation.name
        print donor_branch_in_reticulation
        print recipient_branch_in_reticulation.get_ancestors()
        print recipient_branch_in_reticulation
        return
        final_transfers.append(transfer)
    return {folder:final_transfers}

def assess_reconciliation(folder):
    reconciliation_data = {'transfer_supports':[], 'mapping_consistencies':[], 'homologue_group':folder}

    if not os.path.isdir('reconciliations/%s' %folder) or not os.path.isfile('reconciliations/%s/aggregate.reconciliation' %folder):
        return

    with cd('reconciliations/%s' %folder):
        print folder

        #
        # load trees with named branches
        gene_tree = ete3.Tree(open('%s.reconciliation1' %folder).readlines()[7], format=1 )
        for node in gene_tree.traverse():
            if node.name == 'm1':
                continue

            if node.is_leaf():
                node.add_feature('genome_name', node.name.split('_')[0])
                continue

            node_name, node_support = re.search('^(m\d+?)(100|\d{2})?$', node.name).groups()
            node.support            = int(node_support if node_support else 0)
            node.name               =     node_name

        reconciliation = open( 'aggregate.reconciliation' ).read()

        number_of_replications = float(re.match( '^Processed (\d+) files', reconciliation).group(1))
        events                 = re.findall('^(m\d+\s=.*)$', reconciliation, re.M)

        flag = False
        for event in events:
            transfer_support, mapping_node, mapping_consistency = re.search( ',\sTransfers\s=\s(\d+)], \[Most Frequent mapping --> (\S+), (\d+) times\]', event ).groups()
            transfer_support    = float(transfer_support)
            mapping_consistency = float(mapping_consistency)

            reconciliation_data['transfer_supports'].append(transfer_support/number_of_replications)
            if transfer_support:
                reconciliation_data['mapping_consistencies'].append(mapping_consistency/transfer_support)

            if transfer_support/number_of_replications >= 1 and transfer_support == mapping_consistency:
                node_name    = event.split()[0]
                try:
                    reticulation = gene_tree.search_nodes(name=node_name)[0]
                except:
                    continue

                if reticulation.support < 95:
                    continue

                if len(set([leaf.genome_name for leaf in reticulation.get_leaves()])) < len(reticulation):
                    continue

                if not flag:
                    flag = True
                    reconciliation_data['tree']      = gene_tree.copy(method='deepcopy')
                    reconciliation_data['transfers'] = []
                reconciliation_data['transfers'].append( event )

        if not 'tree' in reconciliation_data:
            return None
        return assess_branch_compatibility(folder, reconciliation_data['transfers'], gene_tree, named_reference_tree)

#with multiprocessing.Pool(processes=6) as pool:
pool = multiprocessing.Pool(processes=10)
results = pool.map(assess_reconciliation, os.listdir('reconciliations/')[100:1000])
pool.close()
pool.join()
print 'yeah'

final_transfers = {}
for element in results:
    if type(element) is not dict or element.values() == [[]]:
        continue

    final_transfers.update(element)

### fig, axs = plt.subplots(nrows=2)
### axs[0].set_title('Donor branch Robinson-Foulds distances')
### sns.distplot(rf_values['donor'],     ax=axs[0])
### axs[1].set_title('Recipient branch Robinson-Foulds distances')
### sns.distplot(rf_values['recipient'], ax=axs[1])
### fig.tight_layout()
### fig.savefig('rf_distances.pdf', dpi=600)
###
### fig, axs = plt.subplots(nrows=2)
### axs[0].set_title('Donor branch normalized Robinson-Foulds distances')
### sns.distplot(rf_values['donor_normed'],     ax=axs[0])
### axs[1].set_title('Recipient branch normalized Robinson-Foulds distances')
### sns.distplot(rf_values['recipient_normed'], ax=axs[1])
### fig.tight_layout()
### fig.savefig('rf_distances_norm.pdf', dpi=600)

fig, axs = plt.subplots(nrows=2)
axs[0].set_title('Ranger HGT confidences')
sns.distplot(transfer_supports,     ax=axs[0])
axs[1].set_title('Ranger mapping consistencies')
sns.distplot(mapping_consistencies, ax=axs[1])
fig.tight_layout()
fig.savefig('ranger_support_metrics.pdf', dpi=600)


# test distance to root and branch support correlation #####
distances_from_root = []
supports            = []
for node in reference_tree.traverse():
    if node.is_leaf() or node.is_root():
        continue
    distances_from_root.append(reference_tree.get_distance(node, topology_only=True))
    supports.append(node.support)

supports            = pd.Series(supports, name='Branch support')
distances_from_root = pd.Series(distances_from_root, name='# of bipartitions from root')
pearson_value = pearsonr(distances_from_root, supports)[0]
fig, ax = plt.subplots()
ax.set_title('Pearson correlation: %f' % round(pearson_value, 4))
sns.regplot(distances_from_root, supports, ax=ax)
fig.tight_layout()
fig.savefig('root_distanceVSsupport-correlation.pdf', dpi=300)
