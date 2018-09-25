#!/usr/bin/env python

import subprocess
import os
import time
from commands import getoutput
import re

os.chdir('/nobackup1b/users/thiberio/phylobayes/test_submission')

#
# required parameters
#
# name of the alignment file in phylip format
alignment     = '001027.nex'
#
# Just just monitor chain, or start from beginning?
monitor_only = True
#
# burnin for for convergence assessment
burnin        = 1000
#
# interval between tested trees for convergence assessment
sampling_each = 100
#
# total number of threads to be used
total_threads = 21

#
# phylobayes parameters
#
phylobayes_threads = (total_threads-1) / 2
def run_chain(alignment, chain_name, num_threads=phylobayes_threads, memmory_size=5, time_limit=170):
    sbatch = 'sbatch -p sched_mit_g4nier -n {num_threads} -N 1 \
--mem={memmory_size}GB --time={num_hours}:00:00 \
-J {job_name} -o {job_output} -e {error_output} \
--wrap="; {command}"'.format(
        num_threads=num_threads, memmory_size=memmory_size, num_hours=time_limit,
        job_name=chain_name, job_output='%s.log' %chain_name, error_output='%s.err' %chain_name,
        command='mpirun -np {num_threads} pb_mpi -d {aln} -lg -cat -dgam 1 {chain_name}'.format(
            num_threads=num_threads, aln=alignment, chain_name=chain_name
        )
    )

    job = subprocess.Popen([sbatch], stdout=subprocess.PIPE, shell=True)
    return job.communicate()[0]

chains = []
#
# if more than two chains should be run, increase the array bellow...
for chain_count in [1,2]:
    # save chain names!
    chain_name = '%s.%i' %(os.path.splitext(alignment)[0], chain_count)
    chains.append(chain_name)

    if monitor_only:
        continue

    #
    # run those bastards!
    print run_chain(alignment, chain_name)

yet_to_converge = True
last_check      = 0
while yet_to_converge:
    #
    # wait time between cycles
    time.sleep(3600)

    #
    # check number of calculated trees
    num_trees = int(re.match('(\d+)\s',
                             getoutput('wc -l %s.treelist' %chains[0])
                            ).group(1))
    #
    # more trees than minimum (burn-in)?
    if num_trees < burnin:
        #
        # if not, sleep again...
        print 'Still within burn-in, %i to go...' %(burnin-num_trees)
        continue

    #
    # more trees than minimum since last check (interval)?
    if num_trees - last_check < sampling_each:
        #
        # if not, sleep again...
        print "Haven't slept snough, %i to go..." %(sampling_each-(num_trees-last_check))
        continue

    last_check = num_trees
    #
    # run bpcomp
    bpcomp     = getoutput('bpcomp -x {burnin} {interval} {chain1} {chain2}'.format(
        burnin=burnin, interval=sampling_each, chain1=chains[0], chain2=chains[1])
    )

    #
    # was any result reported?
    if 'empty tree collection' in bpcomp:
        #
        # if not, sleep again...
        print "Chain still too short..."
        continue

    #
    # capture reported maxdiff...
    maxdiff = float(re.search('^maxdiff\s+:\s+(\d+)$', bpcomp, re.M).group(1))
    # is it enough?!?!
    if maxdiff <= 0.3:
        #
        # if yes, end the loop!!! \o\ \o| |o| |o/ /o/
        yet_to_converge = False
    else:
        #
        # if not, sleep again...
        print "Yet to converge, maxdiff is %.4f..." %maxdiff

#
# my work is done, just end the chains!!!!
for chain_name in chains:
    handle = open('%s.run' %chain_name, 'w')
    handle.write('0')
    handle.close()

#
# see ya!
#
