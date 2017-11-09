#!/usr/bin/env python
#encoding: utf-8

############################################################
#                                                          #
# Script to assess taxonomy available in Silva and define  #
#     genomes to sampled...                                #
#                                                          #
#                                       L. Thibério Rangel #
#                                     lthiberiol@gmail.com #
#                                                          #
############################################################

from sys import exit

if __name__ != '__main__':
    exit()

import ete3
import re
import time
import xml.etree.ElementTree as et
import pandas as pd
import os
import pickle as pkl
from copy import deepcopy
import numpy as np
import multiprocessing
import random
import itertools
import ftplib as ftp
from Bio import SeqIO, SearchIO, AlignIO, Align, Alphabet
from popen2 import popen2
from commands import getoutput

aln_alphabet = Alphabet.Gapped( Alphabet.IUPAC.ambiguous_dna )

def qsub_raxml( fasta_file ):
    # Open a pipe to the qsub command.
    qsub_out, qsub_in = popen2('qsub')
     
    # Customize your options here
    job_name = "raxml_%s" %fasta_file
    walltime = "170:00:00"
    processors = "nodes=1:ppn=40,vmem=20gb"
    command = "/usr/lib/raxml/raxmlHPC-PTHREADS-SSE3 -m PROTGAMMAAUTO -p 12345 -s %s -n %s -N 5 -T 40" %( fasta_file, fasta_file.replace('.aln', '') ) 
 
    job_string = """#!/bin/bash
    #PBS -N %s
    #PBS -m n
    #PBS -l walltime=%s
    #PBS -l %s
    #PBS -W group_list=hgt
    #PBS -o ./output/%s.out
    #PBS -e ./error/%s.err
    cd $PBS_O_WORKDIR
    %s""" % (job_name, walltime, processors, job_name.split('/')[-1], job_name.split('/')[-1], command)
     
    # Send job_string to qsub
    qsub_in.write(job_string)
    qsub_in.close()
     
    # Print your job and the system response to the screen as it's submitted
    return qsub_out.read()

def qsub_mafft( fasta_file ):
    # Open a pipe to the qsub command.
    qsub_out, qsub_in = popen2('qsub')
     
    # Customize your options here
    job_name = "mafft_%s" %fasta_file
    walltime = "2:00:00"
    processors = "nodes=1:ppn=1"
    command = "mafft --reorder --auto %s > %s" %( fasta_file, fasta_file.replace('.fasta', '.aln') )
 
    job_string = """#!/bin/bash
    #PBS -N %s
    #PBS -m n
    #PBS -l walltime=%s
    #PBS -l %s
    #PBS -W group_list=hgt
    #PBS -o ./output/%s.out
    #PBS -e ./error/%s.err
    cd $PBS_O_WORKDIR
    %s
    echo 'yeah!'""" % (job_name, walltime, processors, job_name.split('/')[-1], job_name.split('/')[-1], command)
     
    # Send job_string to qsub
    qsub_in.write(job_string)
    qsub_in.close()
     
    # Print your job and the system response to the screen as it's submitted
    return qsub_out.read()

hmm_folder        = '/mnt/work2/hgt/greg/ribo_db'
faa_folder        = '/mnt/work2/hgt/greg/fournierLab/formated_faas'
hmm_search_folder = '/mnt/work2/hgt/greg/fournierLab/ribo_searches'

profiles = [hmm.replace('.hmm', '') for hmm in os.listdir(hmm_folder)]
genomes  = [faa.replace('.faa', '') for faa in os.listdir(faa_folder)]

marker_sequences = pd.DataFrame( columns=profiles, index=genomes )
for search in os.listdir( hmm_search_folder ):
    if search.startswith('.') or search.endswith('.swp'):
        continue
    genome, profile = search.split('-search_for-')
    result          = SearchIO.read( '%s/%s' %(hmm_search_folder, search), 'hmmer3-text' )
    for hit in result.iterhits():
        if hit.evalue <= 1e-20:
            marker_sequences.loc[genome, profile] = hit.id.split('|')[1]
            break
marker_sequences.drop('Ribosomal_L16-archaea', axis=1, inplace=True)
marker_sequences.dropna(how='all', thresh=len(genomes)*0.9, axis=1, inplace=True)
print '\t** Matrix filled!'

gene_family_handles = { profile_name:open('gene_families/%s.fasta' %profile_name, 'wb') 
                            for profile_name in marker_sequences.columns
                      }
for index, row in marker_sequences.iterrows():
    print index
    fasta = SeqIO.parse( '%s/%s.faa' %(faa_folder, index), 'fasta' )
    for protein in fasta:
        if any(row == protein.id.split('|')[1]):
            gene_family = row[row == protein.id.split('|')[1]].iteritems().next()[0]
            SeqIO.write( protein, gene_family_handles[gene_family], 'fasta' )

for gene_family, handle in gene_family_handles.items():
    handle.close()
print '\t** FASTAs of gene families created!'

gene_family_handles = {}
for fasta in os.listdir('gene_families/'):
    if fasta.endswith('.fasta'):
        gene_family_handles[ fasta.replace('.fasta', '.aln') ] = open('gene_families/%s' %fasta)

job_ids = []
for handle in gene_family_handles.values():
    job_ids.append( qsub_mafft( handle.name ).strip() )

print '\t\t**MAFFT job ids: %s' %', '.join(job_ids) 

error_jobs = []
while job_ids:
    for job_id in job_ids:
        job_state = getoutput( 'qstat -f %s | grep job_state' %job_id ).strip()
        if job_state.endswith('C'):
            print '\t\t**%s done!' %job_id
            job_ids.remove( job_id )
        elif not job_state.endswith('R'):
            print '\t\t**%s is weird... you should check!' %job_id
            error_jobs.append( job_id )
            job_ids.remove( job_id )
    time.sleep( 5 )
print '\t**MAFFT alignments done!'

missing_genes = {}
concatenation = {}
for genome in genomes:
    missing_genes[genome]             = 0
    concatenation[genome]             = Align.SeqRecord( Align.Seq('', aln_alphabet) )
    concatenation[genome].name        = genome
    concatenation[genome].id          = genome
    concatenation[genome].description = genome

total_genes = 0.0
for alignment in os.listdir('gene_families'):
    if not alignment.endswith('.aln'):
        continue

    total_genes += 1
    tmp_aln    = AlignIO.read( 'gene_families/%s' %alignment, 'fasta' )
    aln_length = tmp_aln.get_alignment_length()

    observed_genome = []
    fucked_up        = False
    for block in tmp_aln.get_all_seqs():
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

    for genome in set(genomes).difference( observed_genome ):
        concatenation[genome] += Align.Seq( '-' * aln_length, aln_alphabet )
        missing_genes[genome] += 1

    if fucked_up:
        print 'problem with mla concatenation: %s' %alignment
        break

for gene_family, num_missing_genes in missing_genes.items():
    if num_missing_genes/total_genes > 0.2:
        print '**\t\t%s: excluded from analysis!' %gene_family
        concatenation.pop( gene_family )

AlignIO.write( Align.MultipleSeqAlignment( concatenation.values() ), 'ribosomal_concat.aln', 'fasta' )
qsub_raxml( 'ribosomal_concat.aln' )
