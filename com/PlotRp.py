#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/2/24

'''


import os
import numpy as np
import pandas as pd
import math
import matplotlib.pylab as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from upsetplot import from_memberships, plot, UpSet

from helpers.Reader import *
from helpers.Functions import *


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


def recall_precision(workdir, lowcov_datasets, figdir):

    svtypes = ['ins', 'del']

    fig, axes = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10, 8))

    for col_idx, svtype in enumerate(svtypes):
        for caller in CALLERS:
            recalls = []
            precisions = []
            labels = []
            for dataset in lowcov_datasets:

                merged_mat = f'{workdir}/minimap2_{dataset}/comasm_insdel/{caller}-pav.minimap2-ngmlr.{svtype}.mat.txt'
                values = read_survivor_gencomp(merged_mat)
                recall = values[0][1] * 100 / values[0][0]
                precision = values[0][1] * 100 / values[1][1]
                recalls.append(recall)
                precisions.append(precision)

                labels.append(PLATMAP[dataset])


            xticks = np.arange(len(labels))
            recall_ax = axes[0][col_idx]
            precision_ax = axes[1][col_idx]

            recall_ax.plot(xticks, recalls, marker=TOOLMARKERS[caller], markersize=9, lw=2, label=TOOLMAP[caller])
            precision_ax.plot(xticks, precisions, marker=TOOLMARKERS[caller], markersize=9, lw=2, label=TOOLMAP[caller])

            recall_ax.set_title(SVTYPEMAP[svtype], fontsize=13)

            if col_idx == 0:
                recall_ax.set_ylabel('Recall')
                precision_ax.set_ylabel('Precision')

            for ax in [recall_ax, precision_ax]:
                ax.set_ylim(0, 100)
                ax.set_yticks(np.linspace(0, 100, 5))
                ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)

                ax.set_xticks(xticks)
                ax.set_xticklabels(labels, rotation=90, fontsize=12)
                if col_idx == 1:
                    precision_ax.legend(ncol=2)

    plt.tight_layout()
    plt.show()
    # fig.savefig(f'{figdir}/rp_caller_insdel.pdf')


def recall_precision_of_merged(workdir, lowcov_datasets, highcov_datasets, figdir):

    svtypes = ['ins', 'del']

    supps = [1, 2, 3, 4]

    fig, axes = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10, 8))
    for supp in supps:
        for col_idx, svtype in enumerate(svtypes):
            recalls = []
            precisions = []
            labels = []
            for dataset in lowcov_datasets:
                merged_mat = f'{workdir}/merged_callers/{dataset}/supp_{supp}/com/{dataset}.pav.{svtype}.{supp}.mat.txt'
                values = read_survivor_gencomp(merged_mat)
                recall = values[0][1] * 100 / values[0][0]
                precision = values[0][1] * 100 / values[1][1]
                recalls.append(recall)
                precisions.append(precision)
                if '27X' in dataset:
                    labels.append(f'Low-{PLATMAP[dataset]}')
                else:
                    labels.append(PLATMAP[dataset])


            for dataset in highcov_datasets:
                merged_mat = f'{workdir}/merged_callers/{dataset}/supp_{supp}/com/{dataset}.pav.{svtype}.{supp}.mat.txt'
                values = read_survivor_gencomp(merged_mat)
                recall = values[0][1] * 100 / values[0][0]
                precision = values[0][1] * 100 / values[1][1]
                recalls.append(recall)
                precisions.append(precision)
                labels.append(f'High-{PLATMAP[dataset]}')

            recall_ax, precision_ax = axes[0][col_idx], axes[1][col_idx]
            if col_idx == 0:
                recall_ax.set_ylabel('Recall', fontsize=13)
                precision_ax.set_ylabel('Precision', fontsize=13)

            recall_ax.set_title(SVTYPEMAP[svtype], fontsize=13)
            xticks = np.arange(len(labels))

            recall_ax.plot(xticks, recalls, lw=2, marker='.', markersize=9, label=f'SUPP={supp}')
            precision_ax.plot(xticks, precisions, lw=2, marker='.', markersize=9, label=f'SUPP={supp}')

            recall_ax.legend()
            recall_ax.set_ylim(10, 100)
            recall_ax.set_yticks(np.linspace(10, 100, 4))
            recall_ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(10, 100, 4)])

            precision_ax.legend()
            precision_ax.set_ylim(40, 100)
            precision_ax.set_yticks(np.linspace(40, 100, 4))
            precision_ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 4)])

            for ax in [recall_ax, precision_ax]:
                ax.set_xticks(xticks)
                ax.set_xticklabels(labels, fontsize=12, rotation=90)

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{figdir}/rp_merged_insdel.pdf')


def recall_precision_of_acc_merged(workdir, lowcov_datasets, highcov_datasets, figdir):


    svtypes = ['ins', 'del']
    for col_idx, svtype in enumerate(svtypes):

        for dataset in lowcov_datasets:
            recalls = []
            precisions = []
            labels = []

            for i, caller in enumerate(CALLERS):
                labels.append([TOOLMAP[caller]])

                compare_outdir = f'{workdir}/minimap2_{dataset}/compav'
                mat_out = f'{compare_outdir}/{caller}.{dataset}.{svtype}.mat.txt'

                values = read_survivor_gencomp(mat_out)
                recall = values[0][1] * 100 / values[0][0]
                precision = values[0][1] * 100 / values[1][1]

                recalls.append(recall)
                precisions.append(precision)

            for supp in [2, 3]:
                compare_dir =  f'{workdir}/merged_callers/{dataset}/supp_{supp}/com/caller_combine'
                for file in os.listdir(compare_dir):
                    if 'mat' in file and svtype in file:
                        callers = [TOOLMAP[ele] for ele in file.split('.')[3].split('-')]
                        labels.append(callers)
                        this_mat_file = f'{compare_dir}/{file}'
                        values = read_survivor_gencomp(this_mat_file)
                        recall = values[0][1] * 100 / values[0][0]
                        precision = values[0][1] * 100 / values[1][1]
                        recalls.append(recall)
                        precisions.append(precision)

            labels.append([TOOLMAP[ele] for ele in CALLERS])
            all_callers_merged_mat = f'{workdir}/merged_callers/{dataset}/supp_4/com/{dataset}.pav.{svtype}.4.mat.txt'
            values = read_survivor_gencomp(all_callers_merged_mat)
            recall = values[0][1] * 100 / values[0][0]
            precision = values[0][1] * 100 / values[1][1]

            recalls.append(recall)
            precisions.append(precision)

            fig = plt.figure(figsize=(5, 3))

            upset_data = from_memberships(labels, recalls)
            ax = UpSet(upset_data).plot(fig)
            ax['intersections'].set_ylabel(f'Recall of {SVTYPEMAP[svtype]}', fontsize=13)
            plt.suptitle(f'{PLATMAP[dataset]}', fontsize=13)
            # plt.show()
            fig.savefig(f'{figdir}/{dataset}_{svtype}_acc_merge_caller_recall.pdf')

'''
def recall_precision(workdir, lowcov_datasets, highcov_datasets, figdir):
    svtypes = ['ins', 'del']

    all_datasets = []
    all_datasets.extend(lowcov_datasets)
    all_datasets.extend(highcov_datasets)

    fig, axes = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10, 8))

    for col_idx, svtype in enumerate(svtypes):
        for caller in CALLERS:
            recalls = []
            precisions = []
            labels = []
            for dataset in all_datasets:
                merged_mat = f'{workdir}/minimap2_{dataset}/compav/{caller}.{dataset}.{svtype}.mat.txt'
                values = read_survivor_gencomp(merged_mat)
                recall = values[0][1] * 100 / values[0][0]
                precision = values[0][1] * 100 / values[1][1]
                recalls.append(recall)
                precisions.append(precision)
                if dataset in lowcov_datasets:
                    if '27X' in dataset:
                        labels.append(f'Low-{PLATMAP[dataset]}')
                    else:
                        labels.append(PLATMAP[dataset])
                else:
                    labels.append(f'High-{PLATMAP[dataset]}')

            xticks = np.arange(len(labels))
            recall_ax = axes[0][col_idx]
            precision_ax = axes[1][col_idx]

            recall_ax.plot(xticks, recalls, marker=TOOLMARKERS[caller], markersize=9, lw=2, label=TOOLMAP[caller])
            precision_ax.plot(xticks, precisions, marker=TOOLMARKERS[caller], markersize=9, lw=2, label=TOOLMAP[caller])

            recall_ax.set_title(SVTYPEMAP[svtype], fontsize=13)

            if col_idx == 0:
                recall_ax.set_ylabel('Recall')
                precision_ax.set_ylabel('Precision')

            for ax in [recall_ax, precision_ax]:
                ax.set_ylim(0, 100)
                ax.set_yticks(np.linspace(0, 100, 5))
                ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)

                ax.set_xticks(xticks)
                ax.set_xticklabels(labels, rotation=90, fontsize=12)
                if col_idx == 1:
                    precision_ax.legend(ncol=2)

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{figdir}/rp_caller_insdel.pdf')

'''