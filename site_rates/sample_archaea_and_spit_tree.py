#!/usr/bin/env python2
#coding: utf-8

import ete3
import os
import re
import pickle as pkl
import itertools
import pandas as pd
import subprocess
import re
from Bio import SeqIO
from popen2 import popen2
import numpy as np
import random
from commands import getoutput
import multiprocessing

########################################################################################################################
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
########################################################################################################################

os.chdir('/work/site_rate')

ncbi = ete3.NCBITaxa()

header = 'assembly_accession bioproject biosample wgs_master refseq_category taxid species_taxid organism_name infraspecific_name isolate version_status assembly_level release_type genome_rep seq_rel_date asm_name submitter gbrs_paired_asm paired_asm_comp ftp_path excluded_from_refseq relation_to_type_material'.split()
assembly_summary                     = pd.read_table('assembly_summary.txt', comment='#', header=None, names=header, dtype={'taxid':str, 'infraspecific_name':str})
assembly_summary['refseq_category']  = assembly_summary['refseq_category'].str.lower()
assembly_summary['assembly_level']   = assembly_summary['assembly_level'].str.lower()
assembly_summary['genome_rep']       = assembly_summary['genome_rep'].str.lower()
assembly_summary.set_index('assembly_accession', inplace=True)

lineages = {}
for index, row in assembly_summary.iterrows():
    tmp_lineage = {j: int(i) for i, j in ncbi.get_rank(ncbi.get_lineage(row.taxid)).items()}
    lineages[index] = tmp_lineage

lineages = pd.DataFrame.from_dict(lineages, dtype=int).T
lineages.drop('no rank', axis=1, inplace=True)
lineages.dropna(how='any', subset=['genus'], inplace=True)

genus_groups    = lineages.groupby(by='genus')
genus_names     = ncbi.get_taxid_translator(lineages.genus.unique())
sampled_genomes = {}
for genus, indexes in genus_groups.groups.items():
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

assembly_summary = assembly_summary.reindex(index=sampled_genomes.values()).join(lineages)
assembly_summary.to_csv('sampled_genomes.tab', sep='\t')

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













































