import os
import pandas as pd
import re
from scipy.spatial.distance import squareform
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression, HuberRegressor
import numpy as np
from sklearn.metrics import mean_squared_error

from scipy.stats import pearsonr, linregress
import seaborn as sns
from Bio import SeqIO, SearchIO, AlignIO, Align, Alphabet
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import train_test_split
import random

import matplotlib.colors as colors
import matplotlib.cm as cmx
import subprocess
import itertools
from Bio import SeqIO, SearchIO, AlignIO, Align, Alphabet
import multiprocessing
from copy import deepcopy

os.chdir('/work/site_rate/ml_dists')

ncbi = ete3.NCBITaxa()
outgroup = pd.read_table('../outgroups.tab', index_col=0, comment='#', header=None)

def read_mldist(filename):
    raw_file = open(filename).read()
    genomes  = re.findall('^(GC[FA]_\d+\.\d+)', raw_file, re.M)
    raw_file = raw_file.split('\n')
    raw_file.pop(0)

    out = open('%s.tab' %filename, 'wb')
    out.write('\t%s\n' % '\t'.join(genomes))
    for line in raw_file:
        out.write('%s\n' % '\t'.join(line.split()))
    out.close()

    tmp = pd.read_table('%s.tab' %filename, index_col=0)

    tmp.drop(labels=outgroup.index, axis=0, inplace=True)
    tmp.drop(labels=outgroup.index, axis=1, inplace=True)

    return tmp.copy()

original           = read_mldist('original_sequence.mldist')
condensed_original = squareform(original.values)
x_plot             = np.linspace(condensed_original.min(),
                                 condensed_original.max())

#pw_dist_bins = [np.percentile(condensed_original, decile) for decile in range(20, 81, 20)]
pw_dist_bins          = np.linspace(condensed_original.min(), condensed_original.max(), 6)
binning               = np.digitize(condensed_original, pw_dist_bins)
binning[binning == 6] = 5
colors                = '#6B242E #30862D #335A99 #782D86 #FF7733'.split()

simulations = {}
for category in range(1,9):
    simulations[category] = read_mldist('%i.mldist' %category)

for category, df in simulations.items():
    print 'Site-rate category %i' %category

    condensed_simulation = squareform(df.values)

    fig, ax = plt.subplots()
    ax.scatter(condensed_original, condensed_simulation, color='black', edgecolor='none', alpha=0.3)

    for bin in set(binning):

        x = condensed_original[  binning == bin]
        y = condensed_simulation[binning == bin]

        x_plot = np.linspace(x.min(),
                             x.max())

        lm_yx = LinearRegression(fit_intercept=False)
        lm_yx.fit(y.reshape(-1,1),
                  x)

        lm_xy = LinearRegression(fit_intercept=False)
        lm_xy.fit(x.reshape(-1,1),
                  y)

        prediction_y = lm_xy.predict(x_plot.reshape(-1,1))
        mse          = mean_squared_error(lm_yx.predict(condensed_simulation.reshape(-1,1)), condensed_original)

        print '\tbin %i mean of squared residuals: %.4f' %(bin, mse)

        ax.plot(x_plot, prediction_y, linewidth=3, color=colors[bin-1], label='mse=%.3f' %mse)

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.035), ncol=5, frameon=False)
    fig.set_size_inches(10, 6)
    fig.tight_layout()
    fig.savefig('%i.png' %category, dpi=300)
    plt.close()
    break