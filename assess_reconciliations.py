#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
from commands import getoutput
import multiprocessing
import itertools
import pickle as pkl
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

os.chdir('/work/Alphas_and_Cyanos')
named_reference_tree = ete3.Tree('rooted_partitions-with_named_branches.treefile', format=1)

def root_like( ref_tree, tree_to_root ):
    outgroup = ''
    for node in sorted( ref_tree.children, key=len ):
        if node.is_leaf():
            putative_outgroups = tree_to_root.get_leaves_by_name(node.name)
            for leaf in putative_outgroups:
                tree_to_root.set_outgroup( leaf )
                if tree_to_root.get_topology_id() == ref_tree.get_topology_id():
                    outgroup = leaf
                    break
        else:
            outgroup_members = []
            for leaf in node:
                outgroup_members.append( tree_to_root.get_leaves_by_name(leaf.name) )

            for outgroup_combination in itertools.product( *outgroup_members ):
                if tree_to_root.check_monophyly( [x.name for x in outgroup_combination], 'name' )[0]:
                    putative_outgroup = tree_to_root.get_common_ancestor( outgroup_combination )
                    tree_to_root.set_outgroup( putative_outgroup )
                    if tree_to_root.get_topology_id() == ref_tree.get_topology_id():
                        outgroup = putative_outgroup
                        break
        if outgroup:
            return tree_to_root

    if not outgroup:
        return None

def assess_branch_compatibility(folder, transfers, gene_tree, named_reference_tree):
    selected_transfers = []
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
        if rf_max and rf/rf_max > 0.2:
            continue

        #
        # assess recipient branch compatibility
        rf, rf_max, names, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = reticulation.robinson_foulds(recipient_branch, attr_t1='genome_name')
        if rf_max and rf/rf_max > 0.2:
            continue

        donor_intersection     = [leaf.name for leaf in reticulation.get_leaves() if leaf.genome_name in donor_branch.get_leaf_names()    ]
        recipient_intersection = [leaf.name for leaf in reticulation.get_leaves() if leaf.genome_name in recipient_branch.get_leaf_names()]

        if not donor_intersection or not recipient_intersection:
            continue

        donor_branch_in_reticulation     = reticulation.get_common_ancestor(donor_intersection    )
        recipient_branch_in_reticulation = reticulation.get_common_ancestor(recipient_intersection)

        if donor_branch_in_reticulation in recipient_branch_in_reticulation.get_ancestors():
            print '%s: recipient nested into donor' %folder
            selected_transfers.append(transfer)

    return {folder:selected_transfers}

def assess_reconciliation(folder):
    reconciliation_data = {'transfer_supports':[], 'mapping_consistencies':[], 'homologue_group':folder}

    if not os.path.isdir('reconciliations/%s' %folder) or not os.path.isfile('reconciliations/%s/aggregate.reconciliation' %folder):
        return

    with cd('reconciliations/%s' %folder):
        #
        # load trees with named branches
        named_gene_tree = ete3.Tree(open('%s.reconciliation1' %folder).readlines()[7], format=1 )

        gene_tree       = ete3.Tree('../../iqtrees/%s.aln.treefile'      %folder, format=0)
        for leaf in gene_tree.get_leaves():
            leaf.name = re.sub('\.\d+', '', leaf.name.replace('GCA_', 'GCA'))
            leaf.add_feature('genome_name', leaf.name.split('_')[0])

        if set(named_gene_tree.get_leaf_names()).issubset(gene_tree.get_leaf_names()):
            gene_tree.prune(named_gene_tree.get_leaf_names())
        else:
            print 'FUUUUUCK'
        gene_tree = root_like(named_gene_tree, gene_tree)

        reconciliation = open( 'aggregate.reconciliation' ).read()
        for node_name, descendant1, descendant2 in re.findall('^(m\d+)\s=\sLCA\[(\S+),\s(\S+)\]', reconciliation, re.M):
            branch      = gene_tree.get_common_ancestor([descendant1, descendant2])
            branch.name = node_name

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
results = pool.map(assess_reconciliation, os.listdir('reconciliations/'))
pool.close()
pool.join()

final_transfers = {}
for element in results:
    if type(element) is not dict or element.values() == [[]]:
        continue
    final_transfers.update(element)

out = open('final_transfers.pkl', 'wb')
pkl.dump(final_transfers, out)
out.close()
print 'yeah'

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
###
### transfer_supports     = []
### mapping_consistencies = []
### for element in results1:
###     if element == None:
###         continue
###     transfer_supports.extend(element[0])
###     mapping_consistencies.extend(element[1])
### fig, axs = plt.subplots(nrows=2)
### axs[0].set_title('Ranger HGT confidences')
### sns.distplot(transfer_supports,     ax=axs[0])
### axs[1].set_title('Ranger mapping consistencies')
### sns.distplot(mapping_consistencies, ax=axs[1])
### fig.tight_layout()
### fig.savefig('ranger_support_metrics.pdf', dpi=600)
###
### farthest_leaf_from_root   = named_reference_tree.get_farthest_leaf(   topology_only=True)[0]
### longest_possible_distance = farthest_leaf_from_root.get_farthest_node(topology_only=True)[1]
### bipartition_distances        = []
### bipartition_distances_normed = []
### donor_distance_to_root       = []
### for transfers in final_transfers.values():
###     for transfer in transfers:
###         donor, recipient = re.search('Mapping --> (\S+), Recipient --> (\S+)$', transfer, re.M).groups()
###         recipient_branch = named_reference_tree.search_nodes(name=recipient)[0]
###         donor_branch     = named_reference_tree.search_nodes(name=donor    )[0]
###         bipartition_distances.append(donor_branch.get_distance(recipient_branch, topology_only=True))
###         bipartition_distances_normed.append(donor_branch.get_distance(recipient_branch, topology_only=True)/longest_possible_distance)
###         donor_distance_to_root.append(named_reference_tree.get_distance(donor_branch, topology_only=True))
###
### fig, axs = plt.subplots(nrows=2)
### axs[0].set_title('# of bipartitions from root to donor')
### sns.distplot(donor_distance_to_root, ax=axs[0])
### axs[1].set_title('# of bipartitions between donor/recipient (normalized)')
### sns.distplot(bipartition_distances_normed, ax=axs[1])
### fig.tight_layout()
### fig.savefig('hgt_distances.pdf', dpi=600)
