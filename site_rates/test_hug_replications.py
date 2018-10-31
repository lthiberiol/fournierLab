import os
import ete3
from datasketch import MinHash, MinHashLSH

os.chdir('/work/site_rate/hug_et_al')

hug_tree = ete3.Tree('ribosomal_concat_ml_tree-no_comments.tre')
iqtree   = ete3.Tree('replication_test.treefile')


missing_hug    = {}
missing_iqtree = {}

for element in set(hug_tree.get_leaf_names()).difference(iqtree.get_leaf_names()):
    missing_hug[element] = MinHash(num_perm=1024)
    for word in element.split('_'):
        missing_hug[element].update(word.encode('utf8'))

for element in set(iqtree.get_leaf_names()).difference(hug_tree.get_leaf_names()):
    missing_iqtree[element] = MinHash(num_perm=1024)
    for word in element.split('_'):
        missing_iqtree[element].update(word.encode('utf8'))

lsh_hug = MinHashLSH(threshold=0.75, num_perm=1024)
tmp     = [lsh_hug.insert(key, value) for key, value in missing_hug.items()]

lsh_iqtree = MinHashLSH(threshold=0.75, num_perm=1024)
tmp        = [lsh_iqtree.insert(key, value) for key, value in missing_iqtree.items()]

hug_matches = set()
for hug_name, value in missing_hug.items():
    iqtree_values = lsh_iqtree.query(value)
    if not iqtree_values:
        continue
    hug_matches.add((hug_name, iqtree_values[0]))

iqtree_matches = set()
for iqtree_name, value in missing_iqtree.items():
    hug_values = lsh_hug.query(value)
    if not hug_values:
        continue
    iqtree_matches.add((hug_values[0], iqtree_name))

reciprocal_matches = hug_matches.intersection(iqtree_matches)
for hug_match, iqtree_match in reciprocal_matches:
    if hug_match in missing_hug:
        missing_hug.pop(hug_match)
    if iqtree_match in missing_iqtree:
        missing_iqtree.pop(iqtree_match)

lsh_hug = MinHashLSH(threshold=0.5, num_perm=1024)
tmp     = [lsh_hug.insert(key, value) for key, value in missing_hug.items()]

lsh_iqtree = MinHashLSH(threshold=0.5, num_perm=1024)
tmp        = [lsh_iqtree.insert(key, value) for key, value in missing_iqtree.items()]

hug_matches = set()
for hug_name, value in missing_hug.items():
    iqtree_values = lsh_iqtree.query(value)
    if not iqtree_values:
        continue
    hug_matches.add((hug_name, iqtree_values[0]))

iqtree_matches = set()
for iqtree_name, value in missing_iqtree.items():
    hug_values = lsh_hug.query(value)
    if not hug_values:
        continue
    iqtree_matches.add((hug_values[0], iqtree_name))

reciprocal_matches.update(hug_matches.intersection(iqtree_matches))

reciprocal_matches.add(('Bacteria_Omnitrophica_WOR-2_uncultured_SMTZ_29',
                        'Bacteria_WOR_2_uncultured_SMTZ_29'))
reciprocal_matches.add(('Bacteria_CPR_Peregrinibacteria_CG_PER_02',
                        'Bacteria_Peregrinibacteria_CG2_30_FULL_Peregrinibacteria_PER_44_17'))
reciprocal_matches.add(('Bacteria_Fibrobacteres_Acidobacteria_Fibrobacteres_Fibrobacteria_Fibrobacterales_CG_Fibrob_01',
                        'Bacteria_Fibrobacteres_Acidobacteria_Fibrobacteres_Fibrobacteria_Fibrobacterales_CG2_30_FULL_Fibrobacteres_45_31'))

for hug_match, iqtree_match in reciprocal_matches:
    if hug_match in missing_hug:
        missing_hug.pop(hug_match)
    if iqtree_match in missing_iqtree:
        missing_iqtree.pop(iqtree_match)

for hug_match, iqtree_match in reciprocal_matches:
    leaf  = hug_tree.get_leaves_by_name(hug_match)[0]
    leaf.name = iqtree_match

hug_tree.write(outfile='ribosomal_concat_ml_tree-no_comments.tre', format=5, dist_formatter='%.20f')

rf       = hug_tree.robinson_foulds(iqtree, unrooted_trees=True)
