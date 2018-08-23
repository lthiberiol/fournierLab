#!/usr/bin/env python2
#coding: utf-8

import ete3
import os
import itertools
import pandas as pd
import subprocess
from Bio import SeqIO, SearchIO, AlignIO, Align, Alphabet
import re
import multiprocessing
from copy import deepcopy
from popen2 import popen2
import pickle as pkl
import numpy as np
import random
from commands import getoutput

########################################################################################################################
aln_alphabet = Alphabet.Gapped(Alphabet.IUPAC.ambiguous_dna)

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

def qsub_hmmsearch(xargs_input_file, num_nodes=1, num_threads=1, memmory_size='1GB', time_limit='170:00:00',
                   working_dir='~', genome_dir='.', profile_dir='.'):
    # Open a pipe to the qsub command.
    qsub_out, qsub_in = popen2('sbatch')

    # Customize your options here
    job_string = """#!/bin/bash
#SBATCH -p sched_mit_g4nier               #greg's partition on the cluster
#SBATCH --nodes={num_nodes}               #number of nodes
#SBATCH --ntasks={num_threads}            #number of cores
#SBATCH --mem={memmory_size}              #max amount of memory
#SBATCH --time={time_limit}               #wall time
#SBATCH -J hmmsearch                      #job name
#SBATCH --output=hmmsearch.out            #output name
#SBATCH --error=hmmsearch.err             #error file name
#SBATCH --mail-user=lthiberiol@gmail.com  #if you want emails upon start/finish

module load engaging/hmmer/3.1b2
cd {working_dir}

cat {input_file} | xargs -n 2 -P {simultaneous_processes} sh -c 'hmmsearch --acc --noali --cpu 1 -o $1_-_$2.hmm_out {profile_dir}/$1 {genome_dir}/$2' sh
""".format(input_file=xargs_input_file, num_nodes=num_nodes, num_threads=num_threads, memmory_size=memmory_size,
           time_limit=time_limit, working_dir=working_dir, simultaneous_processes=num_threads*num_nodes,
           genome_dir=genome_dir, profile_dir=profile_dir)

    # Send job_string to qsub
    qsub_in.write(job_string)
    qsub_in.close()

    # Print your job and the system response to the screen as it's submitted
    return qsub_out.read()

def run_hmmsearch(xargs_input_file, num_threads=1, working_dir='~', genome_dir='.', profile_dir='.'):
    # Customize your options here
    return os.system("cat {input_file} | xargs -n 2 -P {num_threads} sh -c 'hmmsearch --acc --noali --cpu 1 -o $1_-_$2.hmm_out {profile_dir}/$1 {genome_dir}/$2' sh".format(
        input_file=xargs_input_file, num_threads=num_threads, working_dir=working_dir,
        genome_dir=genome_dir, profile_dir=profile_dir))

def run_mafft(input_file):
    output = open('%s.aln' %input_file, 'w')
    signal = subprocess.call(['mafft', '--auto', '--reorder', input_file], stdout=output)
    output.close()
    return signal
########################################################################################################################

os.chdir('/work/site_rate')

ncbi = ete3.NCBITaxa()

#
# already done!
#
header = 'assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'.split()
assembly_summary                    = pd.read_table('assembly_summary.txt', comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str})
assembly_summary['refseq_category'] = assembly_summary['refseq_category'].str.lower()
assembly_summary['assembly_level']  = assembly_summary['assembly_level'].str.lower()
assembly_summary['genome_rep']      = assembly_summary['genome_rep'].str.lower()
assembly_summary.set_index('assembly_accession', inplace=True)

assembly_summary_outgroup                    = pd.read_table('outgroups.tab', comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str})
assembly_summary_outgroup['refseq_category'] = assembly_summary_outgroup['refseq_category'].str.lower()
assembly_summary_outgroup['assembly_level']  = assembly_summary_outgroup['assembly_level'].str.lower()
assembly_summary_outgroup['genome_rep']      = assembly_summary_outgroup['genome_rep'].str.lower()
assembly_summary_outgroup.set_index('assembly_accession', inplace=True)

lineages = {}
for index, row in assembly_summary.iterrows():
    tmp_lineage = {j: int(i) for i, j in ncbi.get_rank(ncbi.get_lineage(row.taxid)).items()}
    lineages[index] = tmp_lineage
for index, row in assembly_summary_outgroup.iterrows():
    tmp_lineage = {j: int(i) for i, j in ncbi.get_rank(ncbi.get_lineage(row.taxid)).items()}
    lineages[index] = tmp_lineage

lineages = pd.DataFrame.from_dict(lineages, dtype=int).T
lineages.drop('no rank', axis=1, inplace=True)
lineages.dropna(how='any', subset=['genus'], inplace=True)

genus_groups    = lineages.groupby(by='genus')
genus_names     = ncbi.get_taxid_translator(lineages.genus.unique())
sampled_genomes = {}
for genus, indexes in genus_groups.groups.items():
    if indexes.intersection(assembly_summary_outgroup.index).shape[0]:
        for index in indexes:
            sampled_genomes[genus_names[genus]] = index
        continue

    tmp_df = assembly_summary.loc[indexes, 'refseq_category assembly_level genome_rep'.split()].copy()
    for filter in itertools.product(['reference genome', 'representative genome', 'na'],
                                     ['complete genome', 'chromosome', 'scaffold', 'contig'],
                                     ['full', 'partial']):
        filter = pd.Series(index='refseq_category assembly_level genome_rep'.split(), data=filter)
        tmp_df = tmp_df[tmp_df==filter].dropna(how='any')
        if not tmp_df.shape[0]:
            tmp_df = assembly_summary.loc[indexes, 'refseq_category assembly_level genome_rep'.split()].copy()
            continue
        else:
            sampled_genomes[genus_names[genus]] = tmp_df.iloc[0].name
            break

assembly_summary = assembly_summary.append(assembly_summary_outgroup)
assembly_summary = assembly_summary.reindex(index=sampled_genomes.values()).join(lineages)
assembly_summary.to_csv('sampled_genomes.tab', sep='\t')
assembly_summary = pd.read_table('sampled_genomes.tab', index_col=0, comment='#', dtype={'taxid':str, 'infraspecific_name':str})

for index, row in assembly_summary.iterrows():
    subprocess.call(['wget', '-P', 'genomes/', '%s/*_protein.faa.gz' %row.ftp_path])

successfully_downloaded = []
for filename in os.listdir('genomes/'):
    if os.path.getsize('genomes/%s' %filename):
        successfully_downloaded.append(re.match('(GCF_\d{9}\.\d+)_', filename).group(1))

fail_download = assembly_summary.index.difference(successfully_downloaded)
duplicated = [genome for genome in set(successfully_downloaded) if successfully_downloaded.count(genome) > 1]

subprocess.call(['find', 'genomes/', '-name', '*.gz', '-exec', 'gunzip', '-d', '{}', ';'])

os.mkdir('faas')
for filename in os.listdir('genomes/'):
    assembly_acc = re.match('(GC[AF]_\d{9}\.\d+)_', filename).group(1)

    faa                = SeqIO.parse('genomes/%s' %filename, 'fasta')
    sequences_to_write = []
    for protein in faa:
        protein.id = '%s|%s' %(assembly_acc, protein.id)
        sequences_to_write.append(protein)

    SeqIO.write(sequences_to_write, 'faas/%s.faa' %assembly_acc, 'fasta')

with cd('hmm_search/'):
    out = open('tmp_xargs_list', 'wb')
    for line in itertools.product(os.listdir('../hmm_models_same_as_smc'), os.listdir('../faas')):
        out.write('%s\n' %' '.join(line))
    out.close()

    print run_hmmsearch('tmp_xargs_list', num_threads=20, working_dir='hmm_search',
                     genome_dir='/work/site_rate/faas', profile_dir='/work/site_rate/hmm_models_same_as_smc')

failed_searches = [filename for filename in os.listdir('hmm_search/') if filename.endswith('.faa.hmm_out') and not os.path.getsize('hmm_search/%s' %filename)]

profiles = [hmm.replace('.hmm', '') for hmm in os.listdir('hmm_models_same_as_smc')]
genomes  = [faa.replace('.faa', '') for faa in os.listdir('faas')]

#
# fill the "marker_sequences" DataFrame with the ID of the sequences matching the HMM profiles in each genome
#
marker_sequences = pd.DataFrame(columns=profiles, index=genomes)
for search in os.listdir('hmm_search/'):
    if search.startswith('.') or not search.endswith('.faa.hmm_out'):
        continue
    profile, genome = search.split('_-_')
    profile         = profile.replace('.hmm', '')
    genome          = genome.replace('.faa.hmm_out', '')
    result          = SearchIO.read( 'hmm_search/%s' %search, 'hmmer3-text' )
    for hit in result.iterhits():
        if hit.evalue <= 1e-10:
            marker_sequences.loc[genome, profile] = hit.id.split('|')[1]
            break
marker_sequences.drop(['GCF_001315945.1', 'GCF_001316065.1'], axis=0, inplace=True)
marker_sequences.dropna(how='any', thresh=117, axis=1, inplace=True)

#
# create file handles for each marker gene, where we are gonna add the sequences from each genome
gene_family_handles = {profile_name:open('gene_families/%s.fasta' %profile_name, 'wb')
                        for profile_name in marker_sequences.columns}
for index, row in marker_sequences.iterrows():
    print index
    fasta = SeqIO.parse( 'faas/%s.faa' %index, 'fasta' )
    for protein in fasta:
        if any(row == protein.id.split('|')[1]):
            gene_family = row[row == protein.id.split('|')[1]].iteritems().next()[0]
            SeqIO.write(protein, gene_family_handles[gene_family], 'fasta')

for gene_family, handle in gene_family_handles.items():
    handle.close()
print '\t** FASTAs of gene families created!'

pool = multiprocessing.Pool(processes=4)
results = pool.map_async(run_mafft, [handle.name for handle in gene_family_handles.values()])
pool.close()
pool.join()

missing_genes = {} # just to keep track of the number of missing marker genes in each genome
concatenation = {}
for genome in genomes:
    missing_genes[genome]             = 0
    concatenation[genome]             = Align.SeqRecord( Align.Seq('', aln_alphabet) )
    concatenation[genome].name        = genome
    concatenation[genome].id          = genome
    concatenation[genome].description = genome

#
# fill the handles with the marker sequences from each genome
total_genes     = 0.0 # keep track of the number of genes added to the concatenation
partition_file  = open('partition_file.txt', 'wb')
partition_start = 1
for group, handle in gene_family_handles.items():

    alignment = '%s.aln' %handle.name

    if not os.path.isfile(alignment):
        continue

    total_genes += 1
    tmp_aln    = AlignIO.read(alignment, 'fasta')
    aln_length = tmp_aln.get_alignment_length() # get the expected size of the alignment so you can compare if all have the same size

    observed_genome = []
    fucked_up        = False
    for block in list(tmp_aln):
        # if this alignment has a different size from the rest, something is reaaaaaly wrong!
        if len(block) != aln_length:
            fucked_up = True
            break

        genome = block.id.split('|')[0]

        if genome in observed_genome:
            fucked_up = True
            break
        else:
            observed_genome.append(genome)
            concatenation[genome] += deepcopy( block.seq )

    partition_file.write('LG, %s = %i-%i\n' % (group, partition_start, partition_start + aln_length - 1))
    partition_start += aln_length

    #
    # add gaps for those genomes missing this gene (same size as the expected alignment)
    for genome in set(genomes).difference( observed_genome ):
        concatenation[genome] += Align.Seq( '-' * aln_length, aln_alphabet )
        missing_genes[genome] += 1

    if fucked_up:
        sys.exit( '\t**Problem with MSA concatenation: %s' %alignment )
partition_file.close()

#
# write the final concatenation
AlignIO.write(Align.MultipleSeqAlignment(concatenation.values()),'ribosomal_concat.aln', 'fasta')
