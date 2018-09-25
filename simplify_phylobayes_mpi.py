#!/usr/bin/env python

import subprocess
import os
import time
from commands import getoutput
import re

os.chdir('/nobackup1b/users/thiberio/phylobayes/test_submission')

alignment = '001027.nex'

def run_chain(alignment, chain_name, num_threads=10, memmory_size=5, time_limit=170):
    sbatch = 'sbatch -p sched_mit_g4nier -n {num_threads} -N 1 \
--mem={memmory_size}GB --time={num_hours}:00:00 \
-J {job_name} -o {job_output} -e {error_output} \
--wrap="; {command}"'.format(
        num_threads=num_threads, memmory_size=memmory_size, num_hours=time_limit,
        job_name=chain_name, job_output='%s.log' %chain_name, error_output='%s.err' %chain_name,
        command='mpirun -np {num_threads} pb_mpi -d {aln} -lg -cat -dgam 1 {chain_name}'.format(
            num_threads=num_threads - 1, aln=alignment, chain_name=chain_name
        )
    )

    job = subprocess.Popen([sbatch], stdout=subprocess.PIPE, shell=True)
    return job.communicate()[0]

chains = []
for chain_count in [1,2]:
    chain_name = '%s.%i' %(os.path.splitext(alignment)[0], chain_count)
    chains.append(chain_name)

    print run_chain(alignment, chain_name)

time.sleep(60)

yet_to_converge = True
burnin          = 1000
sampling_each   = 100
last_check      = 0
while yet_to_converge:
    num_trees = int(re.match('(\d+)\s',
                             getoutput('wc -l %s.treelist' %chains[0])
                            ).group(1))
    if num_trees < burnin + sampling_each:
        print 'Still within burn-in, %i to go...' %(burnin-num_trees)
        time.sleep(3600)
        continue

    if num_trees - last_check < sampling_each:
        print "Haven't slept snough, %i to go..." %(sampling_each-(num_trees-last_check))
        time.sleep(3600)
        continue

    last_check = num_trees
    bpcomp     = getoutput('bpcomp -x {burnin} {interval} {chain1} {chain2}'.format(
        burnin=burnin, interval=sampling_each, chain1=chains[0], chain2=chains[1])
    )

    maxdiff = float(re.search('^maxdiff\s+:\s+(\d+)$', bpcomp, re.M).group(1))
    if maxdiff <= 0.3:
        yet_to_converge = False
    else:
        print "Yet to converge, maxdiff is %.4f..." %maxdiff
        time.sleep(3600)


