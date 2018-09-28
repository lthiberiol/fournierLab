#!/usr/bin/env python2
#coding: utf-8

import ete3
import os
import pandas as pd
from Bio import SeqIO, SearchIO, AlignIO, Align, Alphabet

os.chdir('/work/site_rate')
ncbi = ete3.NCBITaxa()

#
# already done!
#
header = '''assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name 
            infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name 
            submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'''.split()
assembly_summary_outgroup                    = pd.read_table('outgroups.tab', comment='#', header=None, names=header,
                                                             dtype={'taxid':str, 'infraspecific_name':str})
assembly_summary_outgroup['refseq_category'] = assembly_summary_outgroup['refseq_category'].str.lower()
assembly_summary_outgroup['assembly_level']  = assembly_summary_outgroup['assembly_level'].str.lower()
assembly_summary_outgroup['genome_rep']      = assembly_summary_outgroup['genome_rep'].str.lower()
assembly_summary_outgroup.set_index('assembly_accession', inplace=True)

sampled_genomes = pd.read_table('sampled_genomes.tab', index_col=0)

tree = ete3.Tree('engaging.treefile', format=1)
if tree.check_monophyly(assembly_summary_outgroup.index, 'name'):
    outgroup = tree.get_common_ancestor(assembly_summary_outgroup.index.tolist())
tree.set_outgroup(outgroup)

sampled_genomes = sampled_genomes.reindex(tree.get_leaf_names())

lineages = {}
for index, row in sampled_genomes.iterrows():
    tmp_lineage = {j: int(i) for i, j in ncbi.get_rank(ncbi.get_lineage(row.taxid)).items()}
    lineages[index] = tmp_lineage

lineages = pd.DataFrame.from_dict(lineages, dtype=int).T
lineages.drop('no rank', axis=1, inplace=True)
lineages.dropna(how='any', subset=['genus'], inplace=True)

archaea_groups = {28889: u'Crenarchaeota', 28890: u'Euryarchaeota', 651137: u'Thaumarchaeota'}
for taxid, taxon in archaea_groups.items():
    leaves = lineages[lineages.phylum==taxid].index.tolist()
    isItMonophyletic, fuckingUp, cladeType = tree.check_monophyly(leaves, 'name')
    if isItMonophyletic:
        print taxon

########################################################################################################################
# parse site rates
#
########################################################################################################################
alignment = AlignIO.read('ribosomal_concat.aln', 'fasta')
#
# GAMMA distribution
#
ratesG = pd.read_table('ratesG', comment='#')
simulated_partitions = {}
for partition_number in ratesG.Part.unique():
    partition                              = ratesG[ratesG.Part == partition_number]
    simulated_partitions[partition_number] = {category:{} for category in partition.Cat.unique()}

    for category in partition.Cat.unique():
        sites                          = partition[partition.Cat == category]
        simulated_partitions[partition_number][category] = {block.name:[] for block in alignment}
        for sequence in alignment:
            simulated_partitions[partition_number][category][sequence.name].append(''.join([sequence[position] for position in sites.index]))

full_sequences = {}
for category in range(1,9):
    full_sequences[category] = {}
    for sequence in alignment:
        full_sequences[category][sequence.name] = ''

for partition_number, categories in simulated_partitions.items():
    print partition_number
    partition = ratesG[ratesG.Part == partition_number]

    for category, sequences in categories.items():
        for header,sequence in sequences.items():
            if header in 'GCF_001315945.1 GCF_001316065.1'.split():
                continue
            full_sequence = ''
            while len(full_sequence) <= partition.shape[0]:
                full_sequence += sequence[0]
            full_sequences[category][header] += full_sequence[:partition.shape[0]]

for category in range(1,9):
    out = open('simulated_alignmentsG/%i.aln' %category, 'wb')
    for sequence in alignment:
        if sequence.name in 'GCF_001315945.1 GCF_001316065.1'.split():
            continue
        out.write('>%s\n%s\n' %(sequence.name, full_sequences[category][sequence.name]))
    out.close()
