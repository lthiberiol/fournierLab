import ete3
import os
import re
from commands import getoutput

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

reference_tree       = ete3.Tree( 'rooted_partitions-with_BB_support.treefile' )
named_reference_tree = ete3.Tree( 'rooted_partitions-with_named_branches.treefile', format=1 )
counter              = 0
supported_transfers  = {}
rf_values            = []
rf_max_values        = []
for folder in os.listdir( 'reconciliations/' ):
    if not os.path.isdir( 'reconciliations/%s' %folder ):
        continue

    with cd( 'reconciliations/%s' %folder ):
        reconciliation_files = os.listdir('.')

        if 'aggregate.reconciliation' not in reconciliation_files:
            continue

        counter += 1

        #
        # load trees with named branches
        gene_tree = ete3.Tree( open('%s.reconciliation1' %folder).readlines()[7], format=1 )
        for node in gene_tree.traverse():
            if node.is_leaf():
                continue

            if node.name == 'm1':
                continue

            node_name, node_support = re.search('^(m\d+?)(100|\d{2})?$', node.name).groups()
            node.support            = int( node_support if node_support else 0 )
            node.name               =      node_name

        reconciliation = open( 'aggregate.reconciliation' ).read()

        number_of_solutions    = int( re.search( '^Total number of optimal solutions: (\d+)$', reconciliation, re.M ).group(1) )
        if number_of_solutions > 100000:
            continue

        number_of_replications = float( re.match( '^Processed (\d+) files', reconciliation ).group(1) )
        events                 = re.findall( '^(m\d+\s=.*)$', reconciliation, re.M )

        flag                = False
        for event in events:
            transfer_support, mapping_node, mapping_consistency = re.search( ',\sTransfers\s=\s(\d+)], \[Most Frequent mapping --> (\S+), (\d+) times\]', event ).groups()
            transfer_support    = float(transfer_support)
            mapping_consistency = float(mapping_consistency)
            mapping_node        = named_reference_tree.search_nodes( name=mapping_node)[0] 

            if mapping_node.is_leaf():
                continue

            if transfer_support/number_of_replications >= 1 and transfer_support == mapping_consistency:
                node_name    = event.split()[0]
                try:
                    reticulation = gene_tree.search_nodes( name=node_name )[0]
                except:
                    continue
                if reticulation.support < 95:
                    continue

                for leaf in reticulation.get_leaves():
                    leaf.add_feature( 'genome_name', leaf.name.split('_')[0] )

                if len(set([leaf.genome_name for leaf in reticulation.get_leaves()])) < len(reticulation):
                    continue

                rf_flag = False
                for child in reticulation.children:
                    rf, rf_max, names, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = reference_tree.robinson_foulds(child, attr_t2='genome_name')
                    if rf:
                        rf_values.append( rf )
                        rf_max_values.append( rf_max )
#                        rf_flag = True
#                        break

                if rf_flag:
                    continue

                if not flag:
                    flag = True
                    supported_transfers[folder] = []
                supported_transfers[folder].append( event )
print 'yeah'

final_transfers = {}
for gene_family, transfers in supported_transfers.items():
    with cd( 'reconciliations/%s' %gene_family ):
        for transfer in transfers:
            grep_query  = re.match( '^(m\d+ = LCA\[\S+, \S+\]:)', transfer, re.M ).group(1)
            grep_result = getoutput( 'grep "%s" %s.reconciliation*' %(re.escape(grep_query), gene_family ) )
            grep_result = set(re.sub( '^%s.reconciliation\d+:' %gene_family, '', grep_result, flags=re.M ).split('\n'))
            if len(grep_result) > 1:
                continue

            transfer         = grep_result.pop()
            donor, recipient = re.search( 'Mapping --> (\S+), Recipient --> (\S+)$', transfer, re.M ).groups()
            recipient_branch = named_reference_tree.search_nodes( name=recipient )[0]
            donor_branch     = named_reference_tree.search_nodes( name=donor )[0]
            if recipient_branch.is_leaf() or donor_branch.is_leaf():
                continue

            if gene_family not in final_transfers:
                final_transfers[gene_family] = []
            final_transfers[gene_family].append( transfer )
print 'yeah'

distances_from_root2 = []
transfer_distance2   = []
for gene_family, transfers in final_transfers.items():
    for transfer in transfers:
        donor, recipient = re.search( 'Mapping --> (\S+), Recipient --> (\S+)$', transfer, re.M ).groups()
        recipient_branch = named_reference_tree.search_nodes( name=recipient )[0]
        donor_branch     = named_reference_tree.search_nodes( name=donor )[0]
        if recipient_branch.is_leaf():
            continue

        distances_from_root2.append( named_reference_tree.get_distance( donor_branch,     topology_only=True ) )
        distances_from_root2.append( named_reference_tree.get_distance( recipient_branch, topology_only=True ) )
        transfer_distance2.append( donor_branch.get_distance( recipient_branch, topology_only=True ) )
print 'yeah'

from scipy.stats import pearsonr
from scipy.stats import spearmanr
supports = []
distances_from_root3 = []
for node in reference_tree.traverse():
	if not node.is_leaf() and not node.is_root():
		supports.append( node.support )
		distances_from_root3.append( reference_tree.get_distance( node, topology_only=True) )
spearmanr( supports, distances_from_root3)
            
for gene_family, transfers in final_transfers.items():
    print gene_family
    print '\t%s' %'\n\t'.join(transfers)
