#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/3/17

'''

import pysam
import vcf
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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


def plot_sv_num(workdir, datasets):

    sv_counts = {'svimasm': [], 'pav': []}


    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['total'])
        if dataset in datasets:
        # minimap2_sv_count[caller][dataset_idx] = svcount
            sv_counts[caller].append(svcount)

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['total'])
        if dataset in datasets:
            sv_counts[caller].append(svcount)
        # minimap2_sv_count[caller][dataset_idx] = svcount

    for caller, count in sv_counts.items():
        std = np.std(count)
        print(f'{caller}: {std}')

def plot_sv_by_regions(workdir, datasets, figdir):
    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    pav_sv_counts = pd.read_csv(f'{workdir}/pav_sv_counts_region.tsv', sep='\t', header=0)
    svimasm_sv_counts = pd.read_csv(f'{workdir}/svimasm_sv_counts_region.tsv', sep='\t', header=0)

    sv_region_dict = {'pav': [[0 for i in range(len(datasets))] for j in range(len(region_labels))],
                      'svimasm':[[0 for i in range(len(datasets))] for j in range(len(region_labels))]}

    for idx, row in pav_sv_counts.iterrows():
        dataset, region, count = row['dataset'], row['region'], int(row['count'])
        print(dataset)
        dataset_idx = datasets.index(dataset)
        region_idx = region_labels.index(region)
        sv_region_dict['pav'][region_idx][dataset_idx] += count

    for idx, row in svimasm_sv_counts.iterrows():
        dataset, region, count = row['dataset'], row['region'], int(row['count'])

        dataset_idx = datasets.index(dataset)
        region_idx = region_labels.index(region)
        sv_region_dict['svimasm'][region_idx][dataset_idx] += count

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 4))
    xticks = np.arange(len(datasets))
    bar_width = 0.5

    for i, caller in enumerate(['pav', 'svimasm']):
        ax = axes[i]
        ax.set_title(TOOLMAP[caller], fontsize=13)
        this_caller_data = sv_region_dict[caller]
        for j, region_label in enumerate(region_labels):
            this_region_count = this_caller_data[j]
            if j == 0:
                ax.bar(xticks, this_region_count, label=region_label, width=bar_width, edgecolor='w')
            else:
                bottoms = []
                for m in range(0, j):
                    bottoms.append([val for val in this_caller_data[m]])
                bottom_sum = [sum(x) for x in zip(*bottoms)]
                ax.bar(xticks, this_region_count, bottom=bottom_sum, label=region_label, width=bar_width, edgecolor='w')
        ax.legend()
        ax.set_xticks(xticks)
        ax.set_xticklabels([PLATMAP[plat] for plat in datasets], fontsize=12, rotation=90)


    plt.tight_layout()
    plt.show()

def plot_insdel(workdir, callers, datasets):

    sv_counts = {'svimasm': [[0 for i in range(len(datasets))], [0 for i in range(len(datasets))]],
                 'pav': [[0 for i in range(len(datasets))], [0 for i in range(len(datasets))]]}

    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        caller, dataset, aligner, svcount, ins_num, del_num = row['caller'], row['dataset'], row['aligner'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        if dataset in datasets:
            dataset_idx = datasets.index(dataset)
            # minimap2_sv_count[caller][dataset_idx] = svcount
            sv_counts[caller][0][dataset_idx] += ins_num * 100 / svcount
            sv_counts[caller][1][dataset_idx] += del_num * 100 / svcount

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        caller, dataset, aligner, svcount, ins_num, del_num = row['caller'], row['dataset'], row['aligner'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        if dataset in datasets:
            dataset_idx = datasets.index(dataset)
            # minimap2_sv_count[caller][dataset_idx] = svcount
            sv_counts[caller][0][dataset_idx] += ins_num * 100 / svcount
            sv_counts[caller][1][dataset_idx] += del_num * 100 / svcount

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 4))
    xticks = np.arange(len(datasets))
    for i, caller in enumerate(callers):
        ax = axes[i]
        ax.set_title(TOOLMAP[caller], fontsize=13)
        ins_pcrt = sv_counts[caller][0]
        del_pcrt = sv_counts[caller][1]
        ax.plot(xticks, del_pcrt, label='DEL', color=SVTYPECOLORS['DEL'], marker='s', lw=2, markersize=9)
        ax.plot(xticks, ins_pcrt, label='INS', color=SVTYPECOLORS['INS'], marker='X', lw=2, markersize=9)

        ax.set_xticks(xticks)
        ax.set_xticklabels([PLATMAP[data] for data in datasets], rotation=90, fontsize=12)

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)])

        ax.legend()

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{iMACASMFIGURE}/insdel_pcrt.pdf')