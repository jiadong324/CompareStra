#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/5/16

'''

import os
import numpy as np
import pandas as pd
import math
import matplotlib.pylab as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


from helpers.Constant import *

sns.set_theme(style="ticks", font="Arial", font_scale=1.0)

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.minor.width"] = 2
plt.rcParams["xtick.labelsize"] = 11

plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.minor.width"] = 2
plt.rcParams["ytick.labelsize"] = 11

plt.rcParams["axes.linewidth"] = 2


def plot_regioned_unique_loci(workdir, datasets, asm_methods, figdir):

    colnames = ['caller', 'asm_method', 'dataset', 'svtype', 'region', 'rptype', 'assm_unique', 'intersects', 'align_unique']

    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    bar_labels = {0: 'Assm uniques', 1: 'Concordant', 2: 'Align uniques'}

    bar_width = 0.3

    fig, axes = plt.subplots(len(CALLERS), len(datasets), sharex='col', sharey='row', figsize=(14, 8))
    r1 = np.arange(len(region_labels))
    r2 = [r + bar_width for r in r1]
    # r3 = [r + bar_width for r in r2]
    xticks = [r + (bar_width + 0.04) / 2 for r in r1]

    rs = [r1, r2]


    results = pd.read_csv(f'{workdir}/cross_loci_byregions.tsv', sep='\t', header=0, names=colnames)

    comasm_region_num = [[{caller: [[0 for j in range(len(region_labels))] for i in range(3)] for caller in CALLERS} for i in range(len(datasets))] for i in range(len(asm_methods))]

    for idx, row in results.iterrows():
        caller, asm_caller, dataset, region, assm_unique, intersect, align_unique = row['caller'], row['asm_method'], row['dataset'], row['region'], \
                                                                        int(row['assm_unique']), int(row['intersects']), int(row['align_unique'])
        label_idx = region_labels.index(region)
        dataset_idx = datasets.index(dataset)
        if asm_caller not in asm_methods:
            continue
        asm_caller_idx = asm_methods.index(asm_caller)

        comasm_region_num[asm_caller_idx][dataset_idx][caller][0][label_idx] = assm_unique
        comasm_region_num[asm_caller_idx][dataset_idx][caller][1][label_idx] = intersect
        comasm_region_num[asm_caller_idx][dataset_idx][caller][2][label_idx] = align_unique

    for row_idx, caller in enumerate(CALLERS):
        for col_idx, dataset in enumerate(datasets):
            for asm_caller_idx in range(len(asm_methods)):
                plot_data = comasm_region_num[asm_caller_idx][col_idx][caller]

                this_ax = axes[row_idx][col_idx]

                for i in range(3):
                    color = COMCOLOR['Assm']
                    if i == 1:
                        color = COMCOLOR['Intersect']
                    elif i == 2:
                        color = COMCOLOR['Align']

                    this_num = plot_data[i]
                    if i == 0:
                        this_ax.bar(rs[asm_caller_idx], this_num, width=bar_width, color=color, label=bar_labels[i])
                    else:
                        bottom = []
                        for m in range(0, i):
                            bottom.append([val for val in plot_data[m]])
                        bottom_sum = [sum(x) for x in zip(*bottom)]
                        this_ax.bar(rs[asm_caller_idx], this_num, width=bar_width, bottom=bottom_sum, color=color, label=bar_labels[i])


                if row_idx == 0:
                    this_ax.set_title(PLATMAP[dataset], fontsize=13)

                if col_idx == 0:
                    this_ax.set_ylabel(f'{TOOLMAP[caller]} SV', fontsize=13)

                this_ax.set_xticks(xticks)
                this_ax.set_xticklabels(region_labels, rotation=90, fontsize=12)

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{figdir}/regioned_unique_loci.pdf')