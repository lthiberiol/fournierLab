#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
import multiprocessing
import pickle as pkl
import linecache
import pandas as pd

os.chdir('/work/Alphas_and_Cyanos')
ncbi = ete3.NCBITaxa()
header = 'assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'.split()
genbank_summary                     = pd.read_table('/work/assembly_summary_genbank.txt', comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str})
genbank_summary['refseq_category']  = genbank_summary['refseq_category'].str.lower()
genbank_summary['assembly_level']   = genbank_summary['assembly_level'].str.lower()
genbank_summary['genome_rep']       = genbank_summary['genome_rep'].str.lower()
genbank_summary.set_index('assembly_accession', inplace=True)
genbank_summary.index               = [re.sub('\.\d+$', '', index).replace('_', '') for index in genbank_summary.index]

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

def rename_branches(reconciliation_file, tree):
    branches = re.findall('^(m\d+) = LCA\[(\S+), (\S+)\]:', reconciliation_file, re.M)
    for name, leaf1, leaf2 in branches:
        node = tree.get_common_ancestor(leaf1, leaf2)
        if node.name:
            print node.get_leaves_by_name()
            break
        node.name = name
    return tree

def root_like( ref_tree, tree2 ):
    tree_to_root = tree2.copy()
    for node in sorted( ref_tree.children, key=len ):
        if node.is_leaf():
            leaf = tree_to_root.get_leaves_by_name(node.name)[0]
            tree_to_root.set_outgroup(leaf)
            break
        else:
            is_it_monophyletic, clade_type, fucking_up = tree_to_root.check_monophyly(node.get_leaf_names(), 'name', unrooted=False)
            if is_it_monophyletic:
                equivalent = tree_to_root.get_common_ancestor(node.get_leaf_names())
                tree_to_root.set_outgroup(equivalent)
            else:
                tree_to_root.set_outgroup(fucking_up.pop())
                equivalent = tree_to_root.get_common_ancestor(node.get_leaf_names())
                tree_to_root.set_outgroup(equivalent)
            break

    return tree_to_root

def parse_aggregated(folder, threshold=0.9, leaves_allowed=False):
    with cd(folder):
        aggregated    = open('aggregated').read()
        gene_tree     = {'named':ete3.Tree(linecache.getline('%s-MAD.ranger_out1' %folder, 8), format=1)}

    gene_tree['support'] = root_like(gene_tree['named'], ete3.Tree('/work/Alphas_and_Cyanos/ranger_input_trees-no_long_branches/%s.tree' %folder))
    gene_tree            = rename_branches(aggregated, gene_tree['support'])

    num_replicates = float(re.match('Processed (\d+) files', aggregated).group(1))

    if not leaves_allowed:
        transfers = re.findall('^(m\d+) = .*, Transfers = [^0]\d+?\], \[Most Frequent mapping --> (n\d+), (\d+) times\], \[Most Frequent recipient --> (n\d+), (\d+) times\].', aggregated, re.M)
    else:
        transfers = re.findall('^(m\d+) = .*, Transfers = [^0]\d+?\], \[Most Frequent mapping --> (\S+), (\d+) times\], \[Most Frequent recipient --> (\S+), (\d+) times\].',   aggregated, re.M)

    supported_transfers = []
    for donor_map, donor, ranger_confidence_donor, recipient, ranger_confidence_recipient in transfers:
        if int(ranger_confidence_donor) < threshold*num_replicates or int(ranger_confidence_recipient) < threshold*num_replicates:
            continue
        supported_transfers.append((donor_map, donor, recipient))

    selected_transfers = []
    for donor_map_name, donor_name, recipient_name in supported_transfers:
        donor_map = gene_tree.search_nodes(name=donor_map_name)[0]
        if donor_map.support < 95:
            continue

        recipient_map_search = re.search('^({children[0]}|{children[1]}).*Most Frequent mapping --> {recipient}'.format(recipient=recipient_name, children=[child.name for child in donor_map.children]), aggregated, re.M)
        if recipient_map_search:
            recipient_map_name = recipient_map_search.group(1)
            selected_transfers.append({'donor':donor_name, 'recipient':recipient_name, 'donor_map':donor_map_name, 'recipient_map':recipient_map_name})

    return {folder:[selected_transfers, gene_tree]}

with cd('reconciliations/mad_roots'):
    yeah = map(parse_aggregated, '000049 000050'.split())
    pool = multiprocessing.Pool(processes=4)
    transfers = pool.map(parse_aggregated, os.listdir('.'))
    pool.close()
    pool.join()

    yeah = {}
    for filtered in transfers:
        if filtered.values() != [[]]:
            yeah.update(filtered)

out = open('aggregated/mad-90_threshold.pkl', 'wb')
pkl.dump(yeah, out)
out.close()

reference_tree = ete3.Tree('rooted_partitions-with_named_branches.treefile', format=1)
def visualize_tree(group, transfers, tree, output_folder='index_transfer_trees/mad-90_threshold'):
    gene_tree = tree.copy()
    tmp_names = {}

    for leaf in gene_tree.get_leaves():
        genome_name, gene_name = leaf.name.split('_')
        leaf.add_feature('genome', genome_name)

    count = 0
    for mappings in transfers:
        donor_map     = gene_tree.search_nodes(name=mappings['donor_map']    )[0]
        recipient_map = gene_tree.search_nodes(name=mappings['recipient_map'])[0]

        donor_branch = reference_tree.search_nodes(name=mappings[ 'donor'])[0]
        if 'annotation' not in donor_branch.features:
            donor_branch.add_feature('annotation', {'donor':[], 'recipient':[]})
        recipient_branch = reference_tree.search_nodes(name=mappings['recipient'])[0]
        if 'annotation' not in recipient_branch.features:
            recipient_branch.add_feature('annotation', {'donor':[], 'recipient':[]})

        donor_branch.annotation['donor'].append([group, count])
        recipient_branch.annotation['recipient'].append([group, count])

        #
        # set donor name
        if not donor_map.name in tmp_names:
            tmp_names[donor_map.name] = {'roles':[], 'hgt_count':[]}
        tmp_names[donor_map.name]['roles'].append('donor%i' %count)
        tmp_names[donor_map.name]['hgt_count'].append(str(count))

        #
        # set recipient name
        if not recipient_map.name in tmp_names:
            tmp_names[recipient_map.name] = {'roles':[], 'hgt_count':[]}
        tmp_names[recipient_map.name]['roles'].append('recipient%i' %count)
        tmp_names[recipient_map.name]['hgt_count'].append(str(count))

        count += 1
        del(donor_map, recipient_map)

    features = {}
    for node_name, feats in tmp_names.items():
        node   = gene_tree.search_nodes(name=node_name)[0]
        tmp_id = str(random.random())
        features[tmp_id] = '[&support=%.2f,hgt_role="%s",hgt_count="%s"]' %(node.support, ', '.join(feats['roles']), ', '.join(feats['hgt_count']))
        node.name = tmp_id


    out  = open('%s/%s.Figtree.tree' %(output_folder, group), 'wb')
    out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(gene_tree))
    for node in gene_tree.traverse():
        if node.is_leaf():
            taxid         = genbank_summary.loc[node.genome, 'taxid']
            lineage       = {j: i for i, j in ncbi.get_rank(ncbi.get_lineage(taxid)).items()}
            lineage_names = ncbi.get_taxid_translator(lineage.values())

            out.write('\t%s ' %(node.name))
            comment = ['source_name="%s"' %ncbi.get_taxid_translator([taxid]).itervalues().next()]
            for rank in ['phylum', 'class', 'order', 'family', 'genus']:
                if rank in lineage:
                    comment.append('tax_%s="%s"' %(rank, lineage_names[lineage[rank]]))
            out.write('[&%s]\n' %' '.join(comment))

        else:
            if node.support and not node.name:
                node.name = '[&support=%.2f]' %node.support

    newick_text = gene_tree.write(format=1)
    newick_text = re.sub('_&support_(\d+\.\d\d)_', '[&support=\\1]', newick_text)
    for key, value in features.items():
        newick_text = newick_text.replace(key, value)
    out.write(';\nend;\n')
    out.write('begin trees;\n\ttree tree_1 = [&R] %s\nend;' %newick_text)
    out.close()

#secure_copy = reference_tree.copy('deepcopy')

out  = open('index_transfer_trees/mad-90_threshold/species_tree.Figtree.tree', 'wb')
out.write("#NEXUS\nbegin taxa;\n\tdimensions ntax=%i;\n\ttaxlabels\n" %len(reference_tree))
tmp_names = {}
for node in reference_tree.traverse():
    if node.is_leaf():
        taxid = genbank_summary.loc[node.name, 'taxid']
        lineage = {j: i for i, j in ncbi.get_rank(ncbi.get_lineage(taxid)).items()}
        lineage_names = ncbi.get_taxid_translator(lineage.values())

        out.write('\t%s ' % (node.name))
        comment = ['source_name="%s"' % ncbi.get_taxid_translator([taxid]).itervalues().next()]
        for rank in ['class', 'phylum', 'order']:
            if rank in lineage:
                comment.append('tax_%s="%s"' % (rank, lineage_names[lineage[rank]]))
        out.write('[&%s]\n' %' '.join(comment))
    else:
        if 'annotation' in node.features:
            tmp_id            = str(random.random())
            tmp_names[tmp_id] = {'ranger_name':node.name}
            node.name         = tmp_id
            for group, count in node.annotation['donor']:
                if group not in tmp_names[tmp_id]:
                    tmp_names[tmp_id][group] = {'donor':[], 'recipient':[]}
                tmp_names[tmp_id][group]['donor'].append(str(count))
            for group, count in node.annotation['recipient']:
                if group not in tmp_names[tmp_id]:
                    tmp_names[tmp_id][group] = {'donor':[], 'recipient':[]}
                tmp_names[tmp_id][group]['recipient'].append(str(count))

newick_text = reference_tree.write(format=1)
for key, value in tmp_names.items():
    comment = []
    for group, roles in value.items():
        if group == 'ranger_name':
            continue
        donor_fragment = ''
        if roles['donor']:
            donor_fragment = 'donor%s' %','.join(roles['donor'])

        recipient_fragment = ''
        if roles['recipient']:
            recipient_fragment = 'recipient%s' %','.join(roles['recipient'])

        comment.append('%s_role="%s %s"' %(group, donor_fragment, recipient_fragment) )

    newick_text = newick_text.replace(key, '[&ranger_name="%s",%s]' %(value['ranger_name'], ','.join(comment)))
out.write(';\nend;\n')
out.write('begin trees;\n\ttree tree_1 = [&R] %s\nend;' %newick_text)
out.close()
