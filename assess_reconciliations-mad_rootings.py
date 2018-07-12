#!/usr/bin/env python2.7
#coding: utf-8

import ete3
import os
import re
import multiprocessing
import pickle as pkl
import linecache
import pandas as pd
import numpy as np
import plotly
import plotly.plotly as ptl
from plotly import graph_objs as go
ptl.sign_in('lthiberiol', 'm15ikp59lt')

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

########################################################################################################################
#                                                                                                                      #
#                                                                                                                      #
#                                                                                                                      #
#                                                                                                                      #
########################################################################################################################

def assess_dtl_dist((group, (transfer_data, gene_tree))):
    dtl_distances = {'donor':[], 'recipient':[]}
    donor_trees     = []
    recipient_trees = []
    for transfer in transfer_data:
        recipient_branch = gene_tree.search_nodes(name=transfer['recipient_map'])[0]
        donor_branch     = recipient_branch.get_sisters()[0]


        donor_trees.append(    donor_branch.write(    format=9))
        recipient_trees.append(recipient_branch.write(format=9))

    #
    # donor compatibility assessment
    os.system( 'cp ranger_input_trees/species_tree.template tmp_ranger-%s.input' %(multiprocessing.current_process().name))
    out = open('tmp_ranger-%s.input' %multiprocessing.current_process().name, 'a')
    out.write('\n'.join(donor_trees))
    out.close()
    os.system('/home/thiberio/ranger/CorePrograms/Ranger-DTL.linux -q -i tmp_ranger-%s.input -o tmp_ranger-%s.output' %(multiprocessing.current_process().name, multiprocessing.current_process().name))
    dtl_distances['donor'].extend([int(reconciliation_cost) for reconciliation_cost in re.findall('^The minimum reconciliation cost is: (\d+)', open('tmp_ranger-%s.output' %multiprocessing.current_process().name).read(), re.M)])

    #
    # recipient compatibility assessment
    os.system( 'cp ranger_input_trees/species_tree.template tmp_ranger-%s.input' %(multiprocessing.current_process().name))
    out = open('tmp_ranger-%s.input' %multiprocessing.current_process().name, 'a')
    out.write('\n'.join(recipient_trees))
    out.close()
    os.system('/home/thiberio/ranger/CorePrograms/Ranger-DTL.linux -q -i tmp_ranger-%s.input -o tmp_ranger-%s.output' %(multiprocessing.current_process().name, multiprocessing.current_process().name))
    dtl_distances['recipient'].extend([int(reconciliation_cost) for reconciliation_cost in re.findall('^The minimum reconciliation cost is: (\d+)', open('tmp_ranger-%s.output' %multiprocessing.current_process().name).read(), re.M)])

    return {group:dtl_distances}

dtl_distances = pkl.load(open('aggregated/mad-90_threshold-donor-recipient_branches_dtl.pkl'))

transfer_distances       = {group:[] for group in transfers.keys()}
for group, (transfer_data, gene_tree) in transfers.items():
    for transfer in transfer_data:

        donor_branch     = reference_tree.search_nodes(name=transfer['donor']    )[0]
        recipient_branch = reference_tree.search_nodes(name=transfer['recipient'])[0]

        transfer_distances[group].append(donor_branch.get_distance(recipient_branch, topology_only=False))

tracer = {'color':[], 'x':[], 'y':[], 'text':[], 'linecolor':'transparent', 'marker':'circle', 'color_means':'Branch length between Donor/Recipient', 'marker_size':10}
for group in dtl_distances.keys():
    for position in range(len(dtl_distances[group]['recipient'])):
        tracer['x'    ].append(dtl_distances[group]['recipient'][position])
        tracer['y'    ].append(dtl_distances[group]['donor'][position])
        tracer['text' ].append('%s-#%i' %(group, position))
        tracer['color'].append(transfer_distances[group][position])

color_range = np.linspace(0, np.max(tracer['color']), 100)
color_bins  = np.digitize(tracer['color'], color_range)
bins        = []
for n in range(100):
    bins.append(go.Scatter(x=[], y=[], mode='markers', text=[], name=str(round(color_range[n], 2)), hoverinfo='text', showlegend=False,
                        marker=dict(size=tracer['marker_size'], color=[], colorscale='RdBu', cmax=np.max(tracer['color']), cmin=np.min(tracer['color']), symbol=tracer['marker'], opacity=1.,
               )))

for position in range(len(tracer['color'])):
    current_bin = color_bins[position] - 1
    for feature in ['x', 'y', 'text']:
        bins[current_bin][feature].append(tracer[feature][position])
    bins[current_bin]['marker']['color'].append(tracer['color'][position])

steps = [dict(label='All',
                method='restyle',
                args=[
                    'visible', [True] * (len(bins) + 1)
                ])
]
for i in range(len(bins)):
    step = dict(label=bins[i]['name'],
                method='restyle',
                args=[
                    'visible', [False] * i + [True] * (len(bins) - i)
                ])
    step['args'][1].append(True)
    steps.append(step)
slider = dict(steps=steps, currentvalue={'prefix':'%s: ' %tracer['color_means']}, pad={'t':50})
bins.append(go.Scatter(x=[np.min(tracer['x']), np.max(tracer['x'])], y=[np.min(tracer['y']), np.max(tracer['y'])], showlegend=False, mode='markers',
                       marker=dict(size=10, color=[0.5], colorscale='RdBu', cmax=np.max(tracer['color']), cmin=np.min(tracer['color']), symbol=tracer['marker'], opacity=0,
                                    colorbar=dict(title=tracer['color_means']))
                       ))

plot_data = go.Data(bins)
layout    = go.Layout(title='Donor/Recipient subtree reconciliation costs', hovermode='closest', width=1200, height=1000,
                      xaxis=dict(title='Recipient subtree reconciliation costs'),
                      yaxis=dict(title='Donor subtree reconciliation costs'),
#                      legend=dict(orientation='h'),
                      sliders=[slider])
fig       = go.Figure(data=plot_data, layout=layout)
plot      = plotly.offline.plot(fig, filename='/Library/WebServer/Documents/indexTransfers/mad-donor_VS_recipient_DTL-transfer_distance.html', auto_open=False)








































