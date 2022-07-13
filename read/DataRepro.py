#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/6/23

'''

import matplotlib.pyplot as plt
import pandas as pd
import pysam
import math
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import vcf

from helpers.Annot import *
from helpers.Functions import *
from com.CheckFDR import *
from helpers.Reader import *


sns.set_theme(style="ticks", font="Arial", font_scale=1.0)

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.minor.width"] = 1
plt.rcParams["xtick.labelsize"] = 11

plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.minor.width"] = 1
plt.rcParams["ytick.labelsize"] = 11

plt.rcParams["axes.linewidth"] = 2



def overview_dataset_repro(workdir, aligners):

    labels = [1, 2, 3, 4, 5, 6]


    xticks = np.arange(len(labels))
    fig, axes = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(9, 4))

    hue_color = [TOOLCOLORS[ele] for ele in CALLERS]

    for col_idx, aligner in enumerate(aligners):
        ax = axes[col_idx]
        supp_pcrt = []
        for caller in CALLERS:
            merged_vcf = f'{workdir}/read_dataset_repro/{caller}.{aligner}.jasmine.merged.vcf'
            supp_dict, merged_total = get_survivor_supp(merged_vcf)
            read_pcrt = []
            for supp in labels:
                count = 0
                if supp in supp_dict:
                    count = supp_dict[supp]
                supp_pcrt.append((count / merged_total, f'SUPP={supp}', caller))
                read_pcrt.append(count / merged_total)

            # ax.plot([r + 0.2 for r in xticks], read_pcrt, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller],lw=2, label=TOOLMAP[caller])

        df_supp_pcrt = pd.DataFrame(supp_pcrt, columns=['pcrt', 'supp', 'caller'])

        box = sns.boxplot(data=df_supp_pcrt, x='supp', y='pcrt', ax=ax)
        for i in range(len(labels)):
            mybox = box.artists[i]
            mybox.set_facecolor('white')
            mybox.set_edgecolor('black')

        sns.stripplot(data=df_supp_pcrt, x='supp', y='pcrt', hue='caller', hue_order=CALLERS, palette=hue_color, ax=ax, size=7)

        ax.set_xticks(xticks)
        ax.set_xticklabels([f'SUPP={supp}' for supp in labels], fontsize=13, rotation=90)

        ax.set_ylim(0, .8)
        ax.set_yticks(np.linspace(0, .8, 5))
        ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, .8, 5)], fontsize=12)

        ax.legend()
        ax.set_xlabel('')
        if col_idx == 0:
            ax.set_ylabel('% of concordant calls', fontsize=13)
        else:
            ax.set_ylabel('', fontsize=13)
        ax.set_title(aligner, fontsize=13)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)


    fig.tight_layout()

    plt.show()
    fig.savefig(f'{READFIGURE}/datasets_repro_of_aligners.pdf')


def plot_inacc_svs(workdir, aligners, overlap_pcrt):

    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    caller_colors = [TOOLCOLORS[caller] for caller in CALLERS]
    caller_legend = [Line2D([0], [0], lw=2, label=TOOLMAP[caller], color=TOOLCOLORS[caller]) for caller in CALLERS]

    fig, axes = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(10, 3))
    fig1, axes1 = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(10, 3))

    for col_idx, aligner in enumerate(aligners):

        inacc_svlen = {'Simple Repeats': [], 'Repeat Masked': [], 'Segment Dup': [], 'Unique': []}

        inacc_stras = {'Simple Repeats': [{}, {}], 'Repeat Masked': [{}, {}], 'Segment Dup': [{}, {}],
                       'Unique': [{}, {}]}

        inacc_svtypes = {'Simple Repeats': {}, 'Repeat Masked': {}, 'Segment Dup': {}, 'Unique': {}}
        # inacc_regions = [[0 for i in range(len(callers))] for j in range(len(region_labels))]

        inacc_counts = [[0 for i in range(len(region_labels))] for caller in CALLERS]

        read_concordants_num = 0

        for caller_idx, caller in enumerate(CALLERS):

            matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])

            for idx, row in matched_info.iterrows():
                start_std, end_std, svtype, svlen, sv_region, rppcrt = abs(float(row['START_STD'])), abs(float(row['END_STD'])), \
                                                               row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE'], float(row['PCRT'])
                region_idx = region_labels.index(sv_region)
                if int(row['SUPP']) == 6:
                    if svtype == 'BND':
                        continue
                    read_concordants_num += 1
                    if start_std > 50 and end_std > 50:
                        if rppcrt >= overlap_pcrt:
                            inacc_svlen[sv_region].append((svlen + 1, svtype, caller))
                            # inacc_regions[region_idx][caller_idx] += 1
                            inacc_counts[caller_idx][region_idx] += 1
                        else:
                            inacc_svlen['Unique'].append((svlen + 1, svtype, caller))
                            # inacc_regions[3][caller_idx] += 1
                            inacc_counts[caller_idx][region_idx] += 1

                        if svtype in inacc_stras[sv_region][1]:
                            inacc_stras[sv_region][1][svtype] += 1
                        else:
                            inacc_stras[sv_region][1][svtype] = 1

                        if svtype in inacc_svtypes[sv_region]:
                            inacc_svtypes[sv_region][svtype] += 1
                        else:
                            inacc_svtypes[sv_region][svtype] = 1

        size_ax = axes[col_idx]

        df_unique_regions = pd.DataFrame(inacc_svlen['Simple Repeats'], columns=['svlen', 'svtype', 'caller'])
        sns.kdeplot(data=df_unique_regions, x='svlen', hue='caller', hue_order=CALLERS, palette=caller_colors,
                    log_scale=True, ax=size_ax)

        size_ax.set_title(aligner, fontsize=13)
        size_ax.set_xlabel('')
        size_ax.set_ylabel('Density', fontsize=13)
        size_ax.set_xlim(10, 100000)
        size_ax.set_ylim(0, 0.4)
        size_ax.set_yticks(np.linspace(0, 0.4, 5))
        size_ax.set_yticklabels([0, 0.1, 0.2, 0.3, 0.4], fontsize=12)
        size_ax.legend(handles=caller_legend)
        size_ax.spines['top'].set_visible(False)
        size_ax.spines['right'].set_visible(False)


        count_ax = axes1[col_idx]
        yticks = np.arange(len(region_labels))
        for j, caller_counts in enumerate(inacc_counts):
            count_ax.plot(caller_counts, yticks, color=TOOLCOLORS[CALLERS[j]],
                          marker=TOOLMARKERS[CALLERS[j]], lw=2, label=CALLERS[j])

        count_ax.set_yticks(yticks)
        count_ax.set_yticklabels(region_labels, fontsize=12)
        count_ax.set_title(aligner, fontsize=13)

        count_ax.set_xlim(0, 1200)
        count_ax.set_xticks(np.linspace(0, 1200, 4))
        count_ax.set_xticklabels([int(val) for val in np.linspace(0, 1200, 4)], fontsize=12)
        count_ax.legend()
        count_ax.spines['top'].set_visible(False)
        count_ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig1.tight_layout()
    plt.show()
    # fig.savefig(f'{READFIGURE}/dataset_repro_inacc_segdup_svsize.pdf')
    fig1.savefig(f'{READFIGURE}/dataset_repro_inacc_count.pdf')



def plot_bpstd(workdir, aligners):

    bpstd_list = []
    for col_idx, aligner in enumerate(aligners):
        for caller in CALLERS:
            matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
            for idx, row in matched_info.iterrows():
                start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])),row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                if svtype in ['INS', 'DEL'] and int(row['SUPP']) == 6:
                    bpstd_list.append((start_std, TOOLMAP[caller], aligner))

    df_bpstd = pd.DataFrame(bpstd_list, columns=['bpstd', 'caller', 'aligner'])

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    aligner_legends = [Patch(label=aligner, color=ALIGNERCOLOR[aligner]) for aligner in aligners]

    sns.histplot(data=df_bpstd, x='bpstd', hue='aligner', log_scale=[False, True], bins=40,
                 hue_order=aligners, palette=[ALIGNERCOLOR[aligner] for aligner in aligners], ax=ax)

    ax.set_xlabel('Breakpoint std.(bp)', fontsize=13)
    ax.legend(handles=aligner_legends)

    ax.set_ylabel('# of SVs', fontsize=13)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    plt.show()
    fig.savefig(f'{READFIGURE}/all_dataset_repro_bpstd_distr_by_aligners.pdf')
