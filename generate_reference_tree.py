#!/usr/bin/env python
#encoding: utf-8

############################################################
#                                                          #
# Script to generate phylogenetic tree based in ribosomal  #
#     genes hmm database                                   #
#                                                          #
#                                       L. Thibério Rangel #
#                                     lthiberiol@gmail.com #
#                                                          #
############################################################

import sys

if __name__ != '__main__':
    sys.exit()

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

#
# parse configuration file parameters
#
with open(sys.argv[1]) as configuration_file:

    #
    # try to find the block containing script's parameters
    parameter_block   = re.search( '^generate_reference_tree.py\s?\{([\s\S]*)\}', configuration_file.read(), re.M ) 
    if not parameter_block:
        sys.exit('Could not find parameters respective to "generate_reference_tree.py", sorry...')

    #
    # go through the block lines and atribute parameters to their respective variables
    for line in parameter_block.group(1).strip().split('\n'):
        if line.startswith('#') or not line.strip():
            continue

        key, value = line.split('=')
        key        = key.strip()
        value      = value.split('#')[0].strip() # ignore everything after a "#"

        if key == 'hmm_folder':
            hmm_folder = value
        elif key == 'faa_folder':
            faa_folder = value
        elif key == 'hmm_search_folder':
            hmm_search_folder = value
        elif key == 'output_folder':
            output_folder = value
        elif key == 'num_threads':
            num_threads = int(value)
        elif key == 'cluster_mem':
            cluster_memory = int(value)

    # check if all required variables were defined
    if not set( 'hmm_folder faa_folder  hmm_search_folder output_folder num_threads cluster_memory'.split() ).issubset( locals() ):
        sys.exit( "Could not define values for the current parameters: %s" ', '.join( set( 'hmm_folder faa_folder  hmm_search_folder output_folder num_threads cluster_memory'.split() ).difference( locals() ) ) )

aln_alphabet = Alphabet.Gapped( Alphabet.IUPAC.ambiguous_dna )

def qsub_raxml( fasta_file ):
    # Open a pipe to the qsub command.
    qsub_out, qsub_in = popen2('qsub')
     
    # Customize your options here
    job_name = "raxml_%s" %fasta_file
    walltime = "170:00:00"
    processors = "nodes=1:ppn=%i,vmem=%igb" %(num_threads, cluster_memory)
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

def qsub_mafft( job_array_file_list, num_threads, num_alignments ):
    # Open a pipe to the qsub command.
    qsub_out, qsub_in = popen2('qsub')
     
    # Customize your options here
    job_name = "mafft_job_array"
    walltime = "2:00:00"
    processors = "nodes=1:ppn=1"
    command = 'mafft --reorder --auto $(awk "NR==$PBS_ARRAYID" %s) > $(awk "NR==$PBS_ARRAYID" %s).aln' %( job_array_file_list, job_array_file_list )
 
    job_string = """#!/bin/bash
    #PBS -N %s
    #PBS -m n
    #PBS -l walltime=%s
    #PBS -l %s
    #PBS -W group_list=hgt
    #PBS -o ./output/%s.out
    #PBS -e ./error/%s.err
    #PBS -t 1-%i%%%i
    cd $PBS_O_WORKDIR
    %s""" % (job_name, walltime, processors, job_name.split('/')[-1], job_name.split('/')[-1], num_alignments, num_threads, command)
     
    # Send job_string to qsub
    qsub_in.write(job_string)
    qsub_in.close()
     
    # Print your job and the system response to the screen as it's submitted
    return qsub_out.read()

#
# get the list of files in each folder
profiles = [hmm.replace('.hmm', '') for hmm in os.listdir(hmm_folder)]
genomes  = [faa.replace('.faa', '') for faa in os.listdir(faa_folder)]

#
# fill the "marker_sequences" DataFrame with the ID of the sequences matching the HMM profiles in each genome
#
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
#
# remove "Ribosomal_L16-archaea", we are dealing with BACTERIA (comment the next line if you wanna search for aerchaeas)
marker_sequences.drop('Ribosomal_L16-archaea', axis=1, inplace=True)
#
# remove any marker gene abscent in more than 10% of the sampled genomes
marker_sequences.dropna(how='all', thresh=len(genomes)*0.9, axis=1, inplace=True)
print '\t** Matrix filled!'

#
# create file handles for each marker gene, where we are gonna add the sequences from each genome
gene_family_handles = { profile_name:open('%s/%s.fasta' %(output_folder, profile_name), 'wb') 
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

#
# create a file listing the name of each created file handle
#   that will be given to the MAFFT job aray
job_array_file_list = open('job_array_file_list.tmp', 'wb')
job_array_file_list.write( '\n'.join( [handle.name for handle in gene_family_handles.values()] ) )
job_array_file_list.close() 

#
# feed the list of file handles to MAFFT job array
job_id     = qsub_mafft( 'job_array_file_list.tmp', num_threads, len(gene_family_handles) )
print '\t\t**MAFFT job ids: %s' %', '.join(job_ids) 

#
# wait all alignments to be done!
job_status = re.findall('job_state = (\S+)', getoutput( 'qstat -f -t %s | grep job_state' %job_id ).strip() )
waiting_status = set( 'R Q H W'.split() )
while waiting_status.isdisjoint( job_status ):
    time.sleep(5) # wait 5 seconds before checking again...
    job_status = re.findall('job_state = (\S+)', getoutput( 'qstat -f -t %s | grep job_state' %job_id ).strip() )
print '\t**MAFFT alignments done!'

#
# start to concatenate each MSA created in a single one
#
# create empty handles for each genome
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
total_genes = 0.0 # keep track of the number of genes added to the concatenation
for alignment in os.listdir(output_folder):
    if not alignment.endswith('.aln'):
        continue

    total_genes += 1
    tmp_aln    = AlignIO.read( '%s/%s' %(output_folder, alignment), 'fasta' )
    aln_length = tmp_aln.get_alignment_length() # get the expected size of the alignment so you can compare if all have the same size

    observed_genome = []
    fucked_up        = False
    for block in tmp_aln.get_all_seqs():
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

    #
    # add gaps for those genomes missing this gene (same size as the expected alignment)
    for genome in set(genomes).difference( observed_genome ):
        concatenation[genome] += Align.Seq( '-' * aln_length, aln_alphabet )
        missing_genes[genome] += 1

    if fucked_up:
        sys.exit( '\t**Problem with MSA concatenation: %s' %alignment )

#
# remove genomes missing more than 20% of the marker genes
for genome, num_missing_genes in missing_genes.items():
    if num_missing_genes/total_genes > 0.2:
        print '**\t\t%s: excluded from analysis!' %genome
        concatenation.pop( genome )

#
# write the final concatenation
AlignIO.write( Align.MultipleSeqAlignment( concatenation.values() ), 'ribosomal_concat.aln', 'fasta' )

#
# submit it to RAxML
qsub_raxml( 'ribosomal_concat.aln' )
