#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/5/11

'''

import matplotlib.pyplot as plt
import pandas as pd
import pysam
import math
from scipy import stats
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
plt.rcParams["xtick.minor.width"] = 2
plt.rcParams["xtick.labelsize"] = 11

plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.minor.width"] = 2
plt.rcParams["ytick.labelsize"] = 11

plt.rcParams["axes.linewidth"] = 2

def print_stra_compare_stats(workdir, aligners, datasets):
    aligner_stra_pcrt = [[] for i in range(len(aligners))]

    concordants = {'HiFi': [], 'ONT': []}
    assembly_unique = {'HiFi': [], 'ONT': []}
    read_unique = {'HiFi': [], 'ONT': []}

    assembly_unique_nums = []
    read_unique_nums = []

    union_counts = {dataset: [] for dataset in datasets}
    assm_unique_counts = {dataset: [] for dataset in datasets}
    read_unique_counts = {dataset: [] for dataset in datasets}

    for col_idx, aligner in enumerate(aligners):
        print(f'\nminimap2-{aligner} ===========')
        stra_matched = [[0 for i in range(len(datasets))] for j in range(3)]

        for dataset in datasets:
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            for caller in CALLERS:
                for asm_caller in ASMCALLERS:
                    merged_vcf = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{asm_caller}.{aligner}-minimap2.{dataset}.jasmine.merged.vcf'
                    suppvec_dict, merged_total = get_survivor_suppvec(merged_vcf)
                    dataset_idx = datasets.index(dataset)
                    union_counts[dataset].append(merged_total)

                    assm_unique_counts[dataset].append(suppvec_dict['10'])
                    read_unique_counts[dataset].append(suppvec_dict['01'])

                    for suppvec, count in suppvec_dict.items():
                        if suppvec == '10':
                            stra_matched[0][dataset_idx] += count
                        elif suppvec == '01':
                            stra_matched[2][dataset_idx] += count
                        else:
                            stra_matched[1][dataset_idx] += count

        print('Total ======')
        sums = [sum(x) for x in zip(*[stra_matched[0], stra_matched[1], stra_matched[2]])]
        assm_pcrt = np.mean([i / j for i, j in zip(stra_matched[0], sums)])
        read_pcrt = np.mean([i / j for i, j in zip(stra_matched[2], sums)])
        con_pcrt = np.mean([i / j for i, j in zip(stra_matched[1], sums)])

        assm_unique_count = sum(stra_matched[0])
        read_unique_count = sum(stra_matched[2])

        assembly_unique_nums.append(assm_unique_count)
        read_unique_nums.append(read_unique_count)

        print('\t% assembly unqiue: ', assm_pcrt)
        print('\t% concordants: ', con_pcrt)
        print('\t% read unqiue: ', read_pcrt)


        print('HiFi ======')
        hifi_sums = [sum(x) for x in zip(*[stra_matched[0][0:3], stra_matched[1][0:3], stra_matched[2][0:3]])]
        hifi_assm_pcrt = np.mean([i / j for i,j in zip(stra_matched[0][0:3], hifi_sums)])
        hifi_read_pcrt = np.mean([i / j for i,j in zip(stra_matched[2][0:3], hifi_sums)])
        hifi_con_pcrt = np.mean([i / j for i, j in zip(stra_matched[1][0:3], hifi_sums)])

        print('\t% assembly unqiue: ', hifi_assm_pcrt)
        print('\t% concordants: ', hifi_con_pcrt)
        print('\t% read unqiue: ', hifi_read_pcrt)

        concordants['HiFi'].append(hifi_con_pcrt)
        assembly_unique['HiFi'].append(hifi_assm_pcrt)
        read_unique['HiFi'].append(hifi_read_pcrt)

        print('ONT ======')
        ont_sums = [sum(x) for x in zip(*[stra_matched[0][3:6], stra_matched[1][3:6], stra_matched[2][3:6]])]
        ont_assm_pcrt = np.mean([i / j for i, j in zip(stra_matched[0][3:6], ont_sums)])
        ont_read_pcrt = np.mean([i / j for i, j in zip(stra_matched[2][3:6], ont_sums)])
        ont_con_pcrt = np.mean([i / j for i, j in zip(stra_matched[1][3:6], ont_sums)])

        print('\t% assembly unqiue: ', ont_assm_pcrt)
        print('\t% concordants: ', ont_con_pcrt)
        print('\t% read unqiue: ', ont_read_pcrt)

        concordants['ONT'].append(ont_con_pcrt)
        assembly_unique['ONT'].append(ont_assm_pcrt)
        read_unique['ONT'].append(ont_read_pcrt)


    print('\nOverall uniques and concordants: ===')
    for plat in ['HiFi', 'ONT']:
        print(f'{plat} ======= ')
        print('\tassemly: ', np.mean(assembly_unique[plat]))
        print('\tread: ', np.mean(read_unique[plat]))
        print('\tconcordants: ', np.mean(concordants[plat]))

    medians = []
    for key, vals in union_counts.items():
        medians.append(np.median(vals))
        print(f'{key} median: {np.median(vals)}')
    print(f'All datasets median: {np.median(medians)}')

    # t_value, p_value = stats.ttest_ind()
    all_assm_uniques = []
    for key, vals in assm_unique_counts.items():
        all_assm_uniques.extend(vals)
        print(f'{key} median: {np.median(vals)}')
    print(f'Assembly-based unique median: {np.median(all_assm_uniques)}')

    all_read_uniques = []
    for key, vals in read_unique_counts.items():
        all_read_uniques.extend(vals)
        print(f'{key} median: {np.median(vals)}')
    print(f'Read-based unique median: {np.median(all_read_uniques)}')

def stra_concordant_overview(workdir, aligners, datasets):

    stra_counts = []
    stra_pcrts = []
    # stra_unique_pcrt = []
    merged_counts = []
    aligner_stra_pcrt = [[] for i in range(len(aligners))]

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for col_idx, aligner in enumerate(aligners):
        stra_matched = [[0 for i in range(len(datasets))] for j in range(3)]
        for dataset in datasets:
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            for caller in CALLERS:
                for asm_caller in ASMCALLERS:
                    for assembler in plat_assemblers[plat]:
                        merged_vcf = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.jasmine.merged.vcf'
                        suppvec_dict, merged_total = get_survivor_suppvec(merged_vcf)
                        dataset_idx = datasets.index(dataset)
                        merged_counts.append((plat, merged_total, PLATMAP[dataset], aligner, assembler,
                                              f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}', TOOLMAP[asm_caller]))
                        for suppvec, count in suppvec_dict.items():
                            if suppvec == '10':
                                stra_matched[0][dataset_idx] += count
                                stra_counts.append((plat, count, 'Assembly-Unique', PLATMAP[dataset], aligner, assembler,
                                                    f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))
                            elif suppvec == '01':
                                stra_matched[2][dataset_idx] += count
                                stra_counts.append((plat, count, 'Read-Unique', PLATMAP[dataset], aligner, assembler,
                                                    f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))
                            else:
                                stra_matched[1][dataset_idx] += count
                                stra_counts.append((plat, count, 'Concordant', PLATMAP[dataset], aligner, assembler,
                                                    f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}', TOOLMAP[asm_caller]))

                        stra_pcrts.append((plat, suppvec_dict['11'] / merged_total * 100, PLATMAP[dataset],
                                           aligner, assembler, f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}', TOOLMAP[asm_caller]))

                        # stra_unique_pcrt.append((plat, (suppvec_dict['10'] + suppvec_dict['01']) / merged_total * 100, PLATMAP[dataset], assembler,
                        #                    aligner, f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}', TOOLMAP[asm_caller]))

        for j in range(len(datasets)):
            stra_matched_pcrt = stra_matched[1][j] / (stra_matched[0][j] + stra_matched[1][j] + stra_matched[2][j])
            aligner_stra_pcrt[col_idx].append(stra_matched_pcrt)

    fig, ax = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(3, 4))
    df_merged_counts = pd.DataFrame(merged_counts, columns=['plat', 'count', 'dataset', 'aligner', 'assembler', 'callers', 'asmcaller'])
    # sns.stripplot(data=df_merged_counts, x='plat', y='count', hue='dataset', palette="Dark2", size=6, ax=axes[0])
    sns.violinplot(data=df_merged_counts, x='plat', y='count', palette=[PLATCOLORS['hifi'], PLATCOLORS['ont']],  ax=ax)


    ax.set_ylabel('Number of merged SVs (x$10^3$)', fontsize=13)
    ax.set_ylim(20000, 50000, 4)
    ax.set_yticks(np.linspace(20000, 50000, 4))
    ax.set_yticklabels([int(val) for val in np.linspace(20, 50, 4)], fontsize=12)
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=13)

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=1.5)
    # ax.legend(title='')

    # for i in [0, 1]:
    #     mybox = box.artists[i]
    #     mybox.set_facecolor('white')
    #     mybox.set_edgecolor('black')
    #
    # sns.stripplot(data=df_merged_counts, x='plat', y='count', hue='asmcaller', hue_order=['PAV', 'SVIM-asm'],
    #               palette=[TOOLCOLORS['pav'], TOOLCOLORS['svimasm']], size=6, ax=axes[1])
    #
    # box = sns.boxplot(data=df_merged_counts, x='plat', y='count', ax=axes[1])
    # for i in [0, 1]:
    #     mybox = box.artists[i]
    #     mybox.set_facecolor('white')
    #     mybox.set_edgecolor('black')
    #
    #
    # for i, ax in enumerate(axes):
    #     if i == 0:
    #         ax.set_ylabel('Number of nonredundant SVs (x$10^3$)', fontsize=13)
    #     else:
    #         ax.set_ylabel('')
    #
    #     ax.set_ylim(25000, 40000, 4)
    #     ax.set_yticks(np.linspace(25000, 40000, 4))
    #     ax.set_yticklabels([int(val) for val in np.linspace(25, 40, 4)], fontsize=12)
    #     ax.set_xlabel('')
    #     ax.set_xticklabels(ax.get_xticklabels(), fontsize=13)
    #
    #     ax.spines['left'].set_visible(False)
    #     ax.spines['right'].set_visible(False)
    #     ax.grid(axis='y', ls='--', color='grey', lw=1.5)
    #     ax.legend(title='')

    fig.tight_layout()

    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    df_stra_pcrts = pd.DataFrame(stra_pcrts, columns=['plat', 'pcrt', 'dataset', 'aligner', 'assembler','callers', 'asmcaller'])
    # sns.stripplot(data=df_stra_pcrts, x='plat', y='pcrt', hue='dataset', palette="Dark2", size=6, ax=axes)
    # assemblers = ['flye', 'shasta', 'hifiasm']
    # assembler_color = [ASS]
    sns.stripplot(data=df_stra_pcrts, x='plat', y='pcrt', hue='assembler', palette="Dark2", size=6, ax=axes[0])
    box = sns.boxplot(data=df_stra_pcrts, x='plat', y='pcrt', ax=axes[0])
    for i in [0, 1]:
        mybox = box.artists[i]
        mybox.set_facecolor('white')
        mybox.set_edgecolor('black')

    aligner_color = [ALIGNERCOLOR[ele] for ele in aligners]
    sns.stripplot(data=df_stra_pcrts, x='plat', y='pcrt', hue='aligner', hue_order=aligners, palette=aligner_color, size=6, ax=axes[1])
    box = sns.boxplot(data=df_stra_pcrts, x='plat', y='pcrt', ax=axes[1])
    for i in [0, 1]:
        mybox = box.artists[i]
        mybox.set_facecolor('white')
        mybox.set_edgecolor('black')


    for i, ax in enumerate(axes):
        ax.set_xlabel('')
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=13)
        if i == 0:
            ax.set_ylabel('Percent of concordant SVs', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_ylim(20, 80, 4)
        ax.set_yticks(np.linspace(20, 80, 4))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(20, 80, 4)], fontsize=12)

        ax.legend(title='')
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig1.tight_layout()

    # df_stra_unique_pcrts = pd.DataFrame(stra_unique_pcrt, columns=['plat', 'pcrt', 'dataset', 'aligner', 'callers'])
    # sns.boxplot(data=df_stra_unique_pcrts, y='pcrt', x='plat', palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=axes[2])
    # sns.stripplot(data=df_stra_unique_pcrts, y='pcrt', x='plat', hue='dataset', palette="Set2", ax=axes[2])
    # axes[2].set_xlabel('')
    # axes[2].set_xticklabels(axes[1].get_xticklabels(), fontsize=13)
    #
    # axes[2].set_ylabel('Percent of concordant SVs', fontsize=13)
    # # axes[2].set_ylim(40, 100, 4)
    # # axes[2].set_yticks(np.linspace(40, 100, 4))
    # # axes[2].set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 4)], fontsize=12)
    # axes[2].legend(title='')



    fig2, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(4, 4))
    hifi_datasets = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb']
    ont_datasets = ['ont_9kb', 'ont_19kb', 'ont_30kb']

    hifi_xticks = np.arange(len(hifi_datasets))
    ont_xticks = np.arange(len(ont_datasets))

    for i in range(len(aligners)):
        axes[0].plot(hifi_xticks, aligner_stra_pcrt[i][0:3], label=aligners[i], lw=2, marker='o')
        axes[1].plot(ont_xticks, aligner_stra_pcrt[i][3:6], label=aligners[i], lw=2, marker='o')

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_xticks(hifi_xticks)
            ax.set_ylim(0.4, 0.7)
            ax.set_yticks(np.linspace(0.4, 0.7, 4))
            ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0.4, 0.7, 4)], fontsize=12)
            ax.set_ylabel('Percent of concordant SVs', fontsize=13)
            ax.set_xticklabels([PLATMAP[val] for val in hifi_datasets], rotation=90, fontsize=13)
        else:
            ax.legend()
            ax.set_ylabel('', fontsize=13)
            ax.set_xticks(ont_xticks)
            ax.set_xticklabels([PLATMAP[val] for val in ont_datasets], rotation=90, fontsize=13)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig2.tight_layout()

    plt.show()

    fig.savefig(f'{Figure4}/strategy_merged_count.pdf')
    fig1.savefig(f'{Figure4}/strategy_matched_pcrt.pdf')
    fig2.savefig(f'{Figure4}/strategy_matched_call_pcrt.pdf')



def plot_stra_compare_stats(workdir, aligners, datasets):

    fig, axes = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(11, 4))
    xticks = np.arange(len(datasets))

    aligner_stra_pcrt = [[] for i in range(len(aligners))]
    stra_counts = []
    stra_pcrts = []
    merged_counts = []


    for col_idx, aligner in enumerate(aligners):
        stra_matched = [[0 for i in range(len(datasets))] for j in range(3)]
        for dataset in datasets:
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            for caller in CALLERS:
                for asm_caller in ASMCALLERS:
                    merged_vcf = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{asm_caller}.{aligner}-minimap2.{dataset}.jasmine.merged.vcf'
                    suppvec_dict, merged_total = get_survivor_suppvec(merged_vcf)
                    dataset_idx = datasets.index(dataset)
                    merged_counts.append((plat, merged_total, PLATMAP[dataset], f'minimap2-{aligner}', f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))
                    for suppvec, count in suppvec_dict.items():
                        if suppvec == '10':
                            stra_matched[0][dataset_idx] += count
                            stra_counts.append((plat, count, 'Assembly-Unique', PLATMAP[dataset], f'minimap2-{aligner}', f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))
                        elif suppvec == '01':
                            stra_matched[2][dataset_idx] += count
                            stra_counts.append((plat, count, 'Read-Unique', PLATMAP[dataset], f'minimap2-{aligner}',
                                                f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))
                        else:
                            stra_matched[1][dataset_idx] += count
                            stra_counts.append((plat, count, 'Concordant', PLATMAP[dataset], f'minimap2-{aligner}',
                                                f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))

                    stra_pcrts.append((plat, suppvec_dict['11'] / merged_total * 100, PLATMAP[dataset], f'minimap2-{aligner}', f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))



        labels = ['Assembly-Unique', 'Concordant', 'Read-Unique']
        this_ax = axes[col_idx]
        for j in range(len(datasets)):
            stra_matched_pcrt = stra_matched[1][j] / (stra_matched[0][j] + stra_matched[1][j] + stra_matched[2][j])
            aligner_stra_pcrt[col_idx].append(stra_matched_pcrt)

        for i in range(3):
            this_data = stra_matched[i]
            if i == 0:
                this_ax.bar(xticks, this_data, label=labels[i])
            else:
                bottoms = []
                for j in range(0, i):
                    bottoms.append(stra_matched[j])
                bottom_sum = [sum(x) for x in zip(*bottoms)]
                this_ax.bar(xticks, this_data, bottom=bottom_sum, label=labels[i])

        this_ax.set_title(f'minimap2-{aligner}')
        this_ax.set_ylim(0, 350000, 6)
        this_ax.set_yticks(np.linspace(0, 350000, 6))
        this_ax.set_yticklabels([int(val) for val in np.linspace(0, 35, 6)])
        if col_idx == 0:
            this_ax.set_ylabel('# of SVs ($10^4$)')

        this_ax.set_xticks(xticks)
        this_ax.set_xticklabels([PLATMAP[val] for val in datasets], rotation=90)
        this_ax.legend()


    fig1, ax1 = plt.subplots(1, 1, figsize=(5, 4))
    for i in range(len(aligners)):
        ax1.plot(xticks, aligner_stra_pcrt[i], label=f'minimap2-{aligners[i]}',lw=2, marker='o')

    ax1.legend()
    ax1.set_ylim(0.4, 0.7)
    ax1.set_yticks(np.linspace(0.4, 0.7, 4))
    ax1.set_yticklabels([f'{int(val * 100)}' for val in np.linspace(0.4, 0.7, 4)])
    ax1.set_ylabel('% of SVs')
    ax1.set_xticks(xticks)
    ax1.set_xticklabels([PLATMAP[val] for val in datasets], rotation=90)

    fig.tight_layout()
    fig1.tight_layout()


    fig2, axes = plt.subplots(1, 2, figsize=(8, 4))
    df_merged_counts = pd.DataFrame(merged_counts, columns=['plat', 'count', 'dataset', 'aligner', 'callers'])
    sns.swarmplot(data=df_merged_counts, x='plat', y='count', hue='dataset', palette="Set2", ax=axes[0])
    sns.boxplot(data=df_merged_counts, x='plat', y='count', palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=axes[0])
    # ax.set_yticklabels(ax.get_yticklabels(), rotation=30)
    axes[0].set_ylabel('Number of SVs (x$10^3$)', fontsize=13)
    axes[0].set_ylim(25000, 40000, 4)
    axes[0].set_yticks(np.linspace(25000, 40000, 4))
    axes[0].set_yticklabels([int(val) for val in np.linspace(25, 40, 4)], fontsize=12)
    axes[0].set_xlabel('')
    axes[0].set_xticklabels(axes[0].get_xticklabels(), fontsize=13)

    df_stra_pcrts = pd.DataFrame(stra_pcrts, columns=['plat', 'pcrt', 'dataset', 'aligner', 'callers'])
    sns.swarmplot(data=df_stra_pcrts, x='plat', y='pcrt', hue='dataset', palette="Set2", ax=axes[1])
    sns.boxplot(data=df_stra_pcrts, x='plat', y='pcrt', palette=[PLATCOLORS['HiFi'], PLATCOLORS['ONT']], ax=axes[1])
    axes[1].set_xlabel('')
    axes[1].set_xticklabels(axes[1].get_xticklabels(), fontsize=13)

    axes[1].set_ylabel('Percent of SVs', fontsize=13)
    axes[1].set_ylim(20, 80, 4)
    axes[1].set_yticks(np.linspace(20, 80, 4))
    axes[1].set_yticklabels([f'{int(val)}%' for val in np.linspace(20, 80, 4)], fontsize=12)

    # for i in range(len(aligners)):
    #     axes[2].plot(aligner_stra_pcrt[i], xticks, label=f'{aligners[i]}',lw=2, marker='o')
    #
    # axes[2].legend(loc='lower right')
    # axes[2].set_xlim(0.4, 0.7)
    # axes[2].set_xticks(np.linspace(0.4, 0.7, 4))
    # axes[2].set_xticklabels([f'{int(val * 100)}%' for val in np.linspace(0.4, 0.7, 4)], fontsize=12)
    # axes[2].set_xlabel('% of SVs', fontsize=13)

    # axes[2].set_xticks(xticks)
    # axes[2].set_xticklabels([PLATMAP[val] for val in datasets], rotation=90)

    fig2.tight_layout()



    # fig3, ax = plt.subplots(1, 1, figsize=(6, 4))
    # df_stra_counts = pd.DataFrame(stra_counts, columns=['plat', 'count', 'stra', 'dataset', 'aligners', 'callers'])
    # sns.boxplot(data=df_stra_counts, y='dataset',)

    plt.show()

    # fig.savefig(f'{Figure4}/strategy_matched_call_counts.pdf')
    # fig1.savefig(f'{Figure4}/strategy_matched_call_pcrt.pdf')
    fig2.savefig(f'{Figure4}/strategy_matched_overview.pdf')


def stra_concordant_bp_overview(workdir, aligners, datasets):
    shift_labels = ['0', '0,10', '10,50', '>50']


    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}
    ins_shift_dict = {'Left': [[] for i in range(len(shift_labels))], 'Right': [[] for i in range(len(shift_labels))]}
    del_shift_dict = {'Left': [[] for i in range(len(shift_labels))], 'Right': [[] for i in range(len(shift_labels))]}

    for col_idx, aligner in enumerate(aligners):
        for dataset in datasets:
            ins_left_bpstd_dict = {shift: 0 for shift in shift_labels}
            ins_right_bpstd_dict = {shift: 0 for shift in shift_labels}
            del_left_bpstd_dict = {shift: 0 for shift in shift_labels}
            del_right_bpstd_dict = {shift: 0 for shift in shift_labels}

            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            ins_total = 0
            del_total = 0
            for caller in CALLERS:
                for asm_caller in ASMCALLERS:
                    for assembler in plat_assemblers[plat]:
                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if svtype == 'INS':
                                ins_total += 1
                                if start_std == 0:
                                    ins_left_bpstd_dict['0'] += 1
                                elif start_std <= 10:
                                    ins_left_bpstd_dict['0,10'] += 1
                                elif start_std <= 50:
                                    ins_left_bpstd_dict['10,50'] += 1
                                else:
                                    ins_left_bpstd_dict['>50'] += 1

                                if end_std == 0:
                                    ins_right_bpstd_dict['0'] += 1
                                elif end_std <= 10:
                                    ins_right_bpstd_dict['0,10'] += 1
                                elif start_std <= 50:
                                    ins_right_bpstd_dict['10,50'] += 1
                                else:
                                    ins_right_bpstd_dict['>50'] += 1

                            if svtype == 'DEL':
                                del_total += 1
                                if start_std == 0:
                                    del_left_bpstd_dict['0'] += 1
                                elif start_std <= 10:
                                    del_left_bpstd_dict['0,10'] += 1
                                elif start_std <= 50:
                                    del_left_bpstd_dict['10,50'] += 1
                                else:
                                    del_left_bpstd_dict['>50'] += 1

                                if end_std == 0:
                                    del_right_bpstd_dict['0'] += 1
                                elif end_std <= 10:
                                    del_right_bpstd_dict['0,10'] += 1
                                elif start_std <= 50:
                                    del_right_bpstd_dict['10,50'] += 1
                                else:
                                    del_right_bpstd_dict['>50'] += 1

                for shift_idx, shift in enumerate(shift_labels):

                    # leftbp_pcrt_list.append((ins_left_bpstd_dict[shift] / ins_total * 100, plat, 'INS', aligner, shift))
                    # leftbp_pcrt_list.append((del_left_bpstd_dict[shift] / del_total * 100, plat, 'DEL', aligner, shift))
                    #
                    # rightbp_pcrt_list.append((ins_right_bpstd_dict[shift] / ins_total * 100, plat, 'INS', aligner, shift))
                    # rightbp_pcrt_list.append((del_right_bpstd_dict[shift] / del_total * 100, plat, 'DEL', aligner, shift))

                    ins_shift_dict['Left'][shift_idx].append(ins_left_bpstd_dict[shift])
                    ins_shift_dict['Right'][shift_idx].append(ins_right_bpstd_dict[shift])

                    del_shift_dict['Left'][shift_idx].append(del_left_bpstd_dict[shift])
                    del_shift_dict['Right'][shift_idx].append(del_right_bpstd_dict[shift])

    # hue_order = ['INS', 'DEL']
    # hue_color = [SVTYPECOLORS[val] for val in hue_order]
    # svtype_legends = [Patch(facecolor=SVTYPECOLORS[val], label=val) for val in ['INS', 'DEL']]
    #
    #
    # df_leftbp = pd.DataFrame(leftbp_pcrt_list, columns=['pcrt', 'plat', 'svtype', 'aligner', 'shift'])
    # df_rightbp = pd.DataFrame(rightbp_pcrt_list, columns=['pcrt', 'plat', 'svtype', 'aligner', 'shift'])
    #
    # fig, axes = plt.subplots(1, 2,  figsize=(6, 4))
    #
    # sns.boxplot(data=df_leftbp, y='shift', x='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, ax=axes[0])
    # axes[0].set_ylabel('Left breakpoint std. (bp)', fontsize=13)
    # sns.boxplot(data=df_rightbp, y='shift', x='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, ax=axes[1])
    # axes[1].set_ylabel('Right breakpoint std. (bp)', fontsize=13)
    #
    # for i, ax in enumerate(axes):
    #
    #     ax.set_xlim(0, 80)
    #     ax.set_xticks(np.linspace(0, 80, 3))
    #     ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 80, 3)], fontsize=12)
    #     ax.set_xlabel('Percent of SVs', fontsize=13)
    #
    #     ax.legend(handles=svtype_legends)
    #
    #     ax.set_yticks(np.arange(3))
    #     ax.set_yticklabels(['0', '0~10', '10~50'], fontsize=12)


    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    shift_colors = [SHIFTCOLORDICT[shift] for shift in shift_labels]
    size = 0.3
    for key, values in ins_shift_dict.items():
        shift_avg = []
        for i, counts in enumerate(values):
            shift_avg.append(np.mean(counts))

        if key == 'Left':
            ax.pie(shift_avg, autopct='%1.1f%%', radius=1, colors=shift_colors,
                        startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))
        else:
            ax.pie(shift_avg, autopct='%1.1f%%', radius=1 - size, colors=shift_colors,
                   startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))

    fig.tight_layout()

    fig1, ax = plt.subplots(1, 1, figsize=(4, 4))
    shift_colors = [SHIFTCOLORDICT[shift] for shift in shift_labels]
    size = 0.3
    for key, values in del_shift_dict.items():
        shift_avg = []
        for i, counts in enumerate(values):
            shift_avg.append(np.mean(counts))

        if key == 'Left':
            ax.pie(shift_avg, autopct='%1.1f%%', radius=1, colors=shift_colors,
                   startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))
        else:
            ax.pie(shift_avg, autopct='%1.1f%%', radius=1 - size, colors=shift_colors,
                   startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))
    fig1.tight_layout()

    fig.savefig(f'{Figure4}/ins_bpstd_avg_pieplot.pdf')
    fig1.savefig(f'{Figure4}/del_bpstd_avg_pieplot.pdf')

    plt.show()

def bpstd_avg_radar_plot(df, ax):

    categories = list(df)[1:]
    N = len(categories)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    # Draw one axe per variable + add labels
    plt.xticks(angles[:-1], categories, fontsize=13)

    # Draw ylabels
    ax.set_rlabel_position(45)
    plt.yticks([10, 30, 50], ["10%", "30%", "50%"], color="grey", size=12)
    plt.ylim(0, 50)

    values = df.loc[0].drop('group').values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, linewidth=2, linestyle='solid', label='INS', color=SVTYPECOLORS['INS'])
    ax.fill(angles, values, SVTYPECOLORS['INS'], alpha=0.1)

    values = df.loc[1].drop('group').values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, linewidth=2, linestyle='solid', label='DEL', color=SVTYPECOLORS['DEL'])
    ax.fill(angles, values, SVTYPECOLORS['DEL'], alpha=0.1)

def stra_concordant_bp_scatter(workdir, aligners, datasets):

    datasets_dict = {'HiFi': ['hifi_10kb', 'hifi_15kb', 'hifi_18kb'],
                     'ONT': ['ont_9kb', 'ont_19kb', 'ont_30kb']}

    shift_labels = ['0,10', '>10']

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    for dataset in datasets:
        acc_ins_bpstd_dict = {shift: 0 for shift in shift_labels}
        acc_del_bpstd_dict = {shift: 0 for shift in shift_labels}

        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        this_dataset_idx = datasets_dict[plat].index(dataset)

        ins_total = 0
        del_total = 0

        this_aligner_acc_ins = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}
        this_aligner_acc_del = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}

        for col_idx, aligner in enumerate(aligners):
            for assembler in plat_assemblers[plat]:

                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:
                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if svtype == 'INS':
                                ins_total += 1
                                if start_std <= 10 and end_std <= 10:
                                    acc_ins_bpstd_dict['0,10'] += 1

                            if svtype == 'DEL':
                                del_total += 1
                                if start_std <= 10 and end_std <= 10:
                                    acc_del_bpstd_dict['0,10'] += 1

                this_aligner_acc_ins[plat][this_dataset_idx] += acc_ins_bpstd_dict['0,10'] / ins_total * 100
                this_aligner_acc_del[plat][this_dataset_idx] += acc_del_bpstd_dict['0,10'] / del_total * 100

            ax.scatter(this_aligner_acc_ins['HiFi'], this_aligner_acc_ins['ONT'], color=SVTYPECOLORS['INS'], edgecolor='black', marker=ALIGNERMK[aligner], s=100)
            ax.scatter(this_aligner_acc_del['HiFi'], this_aligner_acc_del['ONT'], color=SVTYPECOLORS['DEL'], edgecolor='black', marker=ALIGNERMK[aligner], s=100)

    all_legends = [Patch(facecolor=SVTYPECOLORS['INS'], label='INS'), Patch(facecolor=SVTYPECOLORS['DEL'], label='DEL')]
    aligner_legends = [Line2D([0], [0], marker=ALIGNERMK[aligner], color='w', mfc='white', mec='black', markersize=10, label=aligner) for aligner in aligners]

    ax.set_xlim(40, 100)
    ax.set_xticks(np.linspace(40, 100, 4))
    ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 4)], fontsize=12)
    ax.set_xlabel('HiFi datasets', fontsize=13)

    ax.set_ylabel('ONT datasets', fontsize=13)
    ax.set_ylim(40, 100)
    ax.set_yticks(np.linspace(40, 100, 4))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 4)], fontsize=12)

    ax.plot([40, 100], [40, 100], color='#c96150', ls='--')

    all_legends.extend(aligner_legends)

    ax.legend(handles=all_legends)
    fig.tight_layout()
    fig.savefig(f'{Figure4}/insdel_bpstd_scatter.pdf')

    plt.show()


def stra_concordant_insdel_bpstd(workdir, aligners, datasets):

    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']
    bpstd_list = []
    inacc_bpstd_list = []

    shift_labels = ['0,10', '10,50', '50,100', '>100']
    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for dataset in datasets:
        ins_bpstd_dict = {shift: 0 for shift in shift_labels}
        del_bpstd_dict = {shift: 0 for shift in shift_labels}

        inacc_ins = {region: 0 for region in regions}
        inacc_del = {region: 0 for region in regions}

        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        ins_total = 0
        del_total = 0

        for col_idx, aligner in enumerate(aligners):
            for assembler in plat_assemblers[plat]:
                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:
                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']


                            if start_std <= 10 and end_std <= 10:
                                if svtype == 'INS':
                                    ins_total += 1
                                    ins_bpstd_dict['0,10'] += 1
                                if svtype == 'DEL':
                                    del_total += 1
                                    del_bpstd_dict['0,10'] += 1
                            else:
                                if svtype == 'INS':
                                    ins_total += 1
                                    inacc_ins[sv_region] += 1
                                if svtype == 'DEL':
                                    del_total += 1
                                    inacc_del[sv_region] += 1

                bpstd_list.append((ins_bpstd_dict['0,10'] / ins_total * 100, plat, aligner, assembler, 'INS'))
                bpstd_list.append((del_bpstd_dict['0,10'] / del_total * 100, plat, aligner, assembler, 'DEL'))
                for region in regions:
                    inacc_bpstd_list.append((inacc_ins[region] / (ins_total - ins_bpstd_dict['0,10']) * 100, plat, aligner, assembler, region, 'INS'))
                    inacc_bpstd_list.append((inacc_del[region] / (del_total - del_bpstd_dict['0,10']) * 100, plat, aligner, assembler, region, 'DEL'))

    hue_order = ['INS', 'DEL']
    hue_color = [SVTYPECOLORS[val] for val in hue_order]
    legends = [Patch(facecolor=SVTYPECOLORS[val], label=val) for val in hue_order]

    df_bpstd = pd.DataFrame(bpstd_list, columns=['pcrt', 'plat', 'aligner', 'assembler', 'svtype'])
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sns.lineplot(data=df_bpstd, x='aligner', y='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color,
                 err_style='band', style='svtype', lw=2, markers=True, ms=10, ax=ax)

    ax.set_ylabel('Percent of SVs', fontsize=13)
    ax.set_ylim(40, 100)
    ax.set_yticks(np.linspace(40, 100, 4))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 4)], fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticks(np.arange(len(aligners)))
    ax.set_xticklabels(aligners, fontsize=13)
    ax.legend(handles=legends)
    ax.set_xlabel('')
    fig.tight_layout()

    fig1, ax = plt.subplots(1, 1, figsize=(3, 4))
    sns.barplot(data=df_bpstd, x='plat', y='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, capsize=.2, ax=ax)
    ax.set_ylabel('Percent of SVs', fontsize=13)
    ax.set_ylim(40, 100)
    ax.set_yticks(np.linspace(40, 100, 4))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 4)], fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=13)
    ax.legend(handles=legends)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('')
    fig1.tight_layout()

    df_inacc_bpstd = pd.DataFrame(inacc_bpstd_list, columns=['pcrt','plat','aligner', 'assembler', 'region', 'svtype'])
    fig2, ax = plt.subplots(1, 1, figsize=(4, 4))
    sns.barplot(data=df_inacc_bpstd, y='region', x='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, capsize=.2, ax=ax)
    ax.set_ylabel('')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)
    ax.set_xlim(0, 100)
    ax.set_xticks(np.linspace(0, 100, 5))
    ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
    ax.set_xlabel('Percent of SVs')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(handles=legends)
    fig2.tight_layout()

    plt.show()
    fig.savefig(f'{Figure4}/strategy_matched_insdel_bpstd_byaligner.pdf')
    fig1.savefig(f'{Figure4}/strategy_matched_insdel_bpstd_byplat.pdf')
    fig2.savefig(f'{Figure4}/strategy_matched_insdel_inacc_bpstd_byregions.pdf')


def stra_concordant_features(workdir, aligners, datasets):

    # regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    matched_types_list = []

    matched_svsize_list = []
    matched_regions_list = []

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for col_idx, aligner in enumerate(aligners):
            for assembler in plat_assemblers[plat]:
                match_types_dict = {}
                matched_regions = {'Simple Repeats': 0, 'Repeat Masked': 0, 'Segment Dup': 0, 'Unique': 0}
                matched_total = 0
                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:

                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])
                        matched_total += len(df_matched)

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                            if svtype in match_types_dict:
                                match_types_dict[svtype] += 1
                            else:
                                match_types_dict[svtype] = 1
                            matched_regions[sv_region] += 1
                            if svlen <= 10000:
                                matched_svsize_list.append((svlen + 1, plat, PLATMAP[dataset], aligner, assembler, f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))

                for type, count in match_types_dict.items():
                    matched_types_list.append((type, count / matched_total * 100, plat))

                for region, count in matched_regions.items():
                    matched_regions_list.append((region, count / matched_total * 100, plat, aligner, assembler))

    fig, axes = plt.subplots(1, 3, sharex='col', sharey='row', figsize=(6, 4))
    df_regions = pd.DataFrame(matched_regions_list, columns=['region', 'pcrt', 'plat', 'aligner', 'assembler'])

    plat_order = ['HiFi', 'ONT']
    plat_color = [PLATCOLORS[plat] for plat in plat_order]
    plat_legends = [Patch(facecolor=PLATCOLORS[plat], label=plat) for plat in plat_order]
    sns.barplot(data=df_regions, y='region', x='pcrt', hue='plat', hue_order=plat_order, palette=plat_color, ax=axes[0])
    axes[0].legend(handles=plat_legends)
    axes[0].set_yticklabels(axes[0].get_yticklabels(), fontsize=13)


    aligner_color = [ALIGNERCOLOR[aligner] for aligner in aligners]
    aligner_legends = [Patch(facecolor=ALIGNERCOLOR[val], label=val) for val in aligners]

    sns.barplot(data=df_regions, y='region', x='pcrt', hue='aligner', palette=aligner_color, ax=axes[1])
    axes[1].legend(handles=aligner_legends)

    sns.barplot(data=df_regions, y='region', x='pcrt', hue='assembler', ax=axes[2])
    axes[2].legend(title='')

    for i, ax in enumerate(axes):
        ax.set_xlabel('Percent of SVs', fontsize=13)
        ax.set_xlim(0, 100)
        ax.set_xticks(np.linspace(0, 100, 3))
        ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 3)], fontsize=12)
        ax.set_ylabel('')
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.tight_layout()

    fig1, axes = plt.subplots(2, 1, figsize=(6, 4))
    df_svsize = pd.DataFrame(matched_svsize_list, columns=['size', 'plat', 'dataset', 'aligner', 'assembler', 'callers'])
    sns.histplot(data=df_svsize[df_svsize['size'] < 1000], x='size', hue='plat',
                 palette=plat_color, ax=axes[0], bins=50, kde=True)
    axes[0].legend(handles=plat_legends)

    sns.histplot(data=df_svsize[df_svsize['size'] >= 1000], x='size', hue='plat',
                 palette=plat_color, ax=axes[1], bins=50, kde=True)
    axes[1].legend(handles=plat_legends)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('# of SVs ($10^4$)', fontsize=13)
            ax.set_ylim(0, 600000)
            ax.set_yticks(np.linspace(0, 600000, 4))
            ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 60, 4)])
        else:
            ax.set_ylim(0, 60000)
            ax.set_yticks(np.linspace(0, 60000, 4))
            ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 60, 4)])
            ax.set_ylabel('# of SVs ($10^3$)', fontsize=13)
        ax.set_xlabel('')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # ax.set_xlim(10, 1000)

    fig1.tight_layout()

    plt.show()

    fig.savefig(f'{Figure4}/strategy_matched_byregions.pdf')
    fig1.savefig(f'{Figure4}/strategy_matched_svsize.pdf')
    # fig4.savefig(f'{Figure4}/strategy_matched_svtypes.pdf')


def wgs_highconf_acc_inacc_overview(workdir, aligners, datasets):

    bpstd_count = {'INS': {'Acc': [], 'Inacc': []}, 'DEL': {'Acc': [], 'Inacc': []}}
    highconf_bpstd_count = {'INS': {'Acc': [], 'Inacc': []}, 'DEL': {'Acc': [], 'Inacc': []}}

    bpstd_pcrt = {'INS': {'Acc': [], 'Inacc': []}, 'DEL': {'Acc': [], 'Inacc': []}}
    highconf_bpstd_pcrt = {'INS': {'Acc': [], 'Inacc': []}, 'DEL': {'Acc': [], 'Inacc': []}}
    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for dataset in datasets:
        ins_bpstd_dict = {'Acc': 0, 'Inacc': 0}
        del_bpstd_dict = {'Acc': 0, 'Inacc': 0}

        highconf_ins_bpstd_dict = {'Acc': 0, 'Inacc': 0}
        highconf_del_bpstd_dict = {'Acc': 0, 'Inacc': 0}

        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        ins_total = 0
        del_total = 0

        highconf_ins_total = 0
        highconf_del_total = 0

        for col_idx, aligner in enumerate(aligners):
            for assembler in plat_assemblers[plat]:

                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:
                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(
                                float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if start_std <= 10 and end_std <= 10:
                                if svtype == 'INS':
                                    ins_total += 1
                                    ins_bpstd_dict['Acc'] += 1
                                if svtype == 'DEL':
                                    del_total += 1
                                    del_bpstd_dict['Acc'] += 1
                            else:
                                if svtype == 'INS':
                                    ins_total += 1
                                    ins_bpstd_dict['Inacc'] += 1
                                if svtype == 'DEL':
                                    del_total += 1
                                    del_bpstd_dict['Inacc'] += 1

                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered_in_highconf/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(
                                float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if start_std <= 10 and end_std <= 10:
                                if svtype == 'INS':
                                    highconf_ins_total += 1
                                    highconf_ins_bpstd_dict['Acc'] += 1
                                if svtype == 'DEL':
                                    highconf_del_total += 1
                                    highconf_del_bpstd_dict['Acc'] += 1
                            else:
                                if svtype == 'INS':
                                    highconf_ins_total += 1
                                    highconf_ins_bpstd_dict['Inacc'] += 1
                                if svtype == 'DEL':
                                    highconf_del_total += 1
                                    highconf_del_bpstd_dict['Inacc'] += 1

        bpstd_count['INS']['Acc'].append(ins_bpstd_dict['Acc'])
        bpstd_count['INS']['Inacc'].append(ins_bpstd_dict['Inacc'])

        bpstd_count['DEL']['Acc'].append(del_bpstd_dict['Acc'])
        bpstd_count['DEL']['Inacc'].append(del_bpstd_dict['Inacc'])

        highconf_bpstd_count['INS']['Acc'].append(highconf_ins_bpstd_dict['Acc'])
        highconf_bpstd_count['INS']['Inacc'].append(highconf_ins_bpstd_dict['Inacc'])

        highconf_bpstd_count['DEL']['Acc'].append(highconf_del_bpstd_dict['Acc'])
        highconf_bpstd_count['DEL']['Inacc'].append(highconf_del_bpstd_dict['Inacc'])

        bpstd_pcrt['INS']['Acc'].append(ins_bpstd_dict['Acc'] / ins_total * 100)
        bpstd_pcrt['INS']['Inacc'].append(ins_bpstd_dict['Inacc'] / ins_total * 100)

        bpstd_pcrt['DEL']['Acc'].append(del_bpstd_dict['Acc'] / del_total * 100)
        bpstd_pcrt['DEL']['Inacc'].append(del_bpstd_dict['Inacc'] / del_total * 100)

        highconf_bpstd_pcrt['INS']['Acc'].append(highconf_ins_bpstd_dict['Acc'] / highconf_ins_total * 100)
        highconf_bpstd_pcrt['INS']['Inacc'].append(highconf_ins_bpstd_dict['Inacc'] / highconf_ins_total * 100)

        highconf_bpstd_pcrt['DEL']['Acc'].append(highconf_del_bpstd_dict['Acc'] / highconf_del_total * 100)
        highconf_bpstd_pcrt['DEL']['Inacc'].append(highconf_del_bpstd_dict['Inacc'] / highconf_del_total * 100)


    bpstd_avg_count = []
    bpstd_avg_pcrt = {'INS': [], 'DEL': []}

    for svtype in ['INS', 'DEL']:
        for flag in ['Acc', 'Inacc']:
            # bpstd_avg_count.append((np.mean(bpstd_count[svtype][flag]), flag,  f'{svtype}-WGS'))
            # bpstd_avg_count.append((np.mean(highconf_bpstd_count[svtype][flag]), flag, f'{svtype}-HighConf'))

            print(f'{flag}-{svtype}-WGS: {np.mean(bpstd_count[svtype][flag])}')
            print(f'{flag}-{svtype}-HighConf: {np.mean(highconf_bpstd_count[svtype][flag])}')

            bpstd_avg_pcrt[svtype].append((np.mean(bpstd_pcrt[svtype][flag]), flag, f'{svtype}-WGS'))
            bpstd_avg_pcrt[svtype].append((np.mean(highconf_bpstd_pcrt[svtype][flag]), flag, f'{svtype}-HighConf'))

    # df_bpstd_avg = pd.DataFrame(bpstd_avg_count, columns=['count', 'flag', 'region'])


    hue_order = ['Acc', 'Inacc']
    legends = [Patch(facecolor='#5ca36d', label='Accurate'), Patch(facecolor='#dc802f', label='Inaccurate')]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 3))
    # bars = sns.barplot(data=df_bpstd_avg, x='region', y='count', hue='flag', hue_order=hue_order, palette=['#5ca36d', '#dc802f'], ax=ax)
    #
    # for i in [0, 2, 4, 6]:
    #     mybar = bars.patches[i]
    #     mybar.set_ls('-')
    #     mybar.set_edgecolor('black')
    #
    # for j in [1, 3, 5, 7]:
    #     mybar = bars.patches[j]
    #     mybar.set_ls('--')
    #     mybar.set_edgecolor('black')
    #
    # ax.set_ylabel('Number of SVs (x$10^3$)', fontsize=13)
    # ax.set_xlabel('')
    # ax.set_xticks([0.5, 2.5])
    # ax.set_xticklabels(['HiFi', 'ONT'], fontsize=13)
    #
    #
    # ax.set_ylim(0, 60000)
    # ax.set_yticks(np.linspace(0, 60000, 5))
    # ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 60, 5)], fontsize=12)

    # ax.legend(handles=legends)
    #
    # ax1 = ax.twinx()
    fig_idx = 0
    for svtype, avg_pcrt in bpstd_avg_pcrt.items():
        df_bpstd_pcrt = pd.DataFrame(avg_pcrt, columns=['pcrt', 'flag', 'region'])
        sns.lineplot(data=df_bpstd_pcrt, x='region', y='pcrt', hue='flag', style='flag', dashes=False, markers=True, ms=8, lw=2,
                     hue_order=hue_order, palette=['#5ca36d', '#dc802f'], ax=axes[fig_idx])
        axes[fig_idx].set_title(svtype, fontsize=13)
        # axes[0].plot([1.5, 2.5], avg_pcrt, lw=2, marker='o', markersize=8, label='')

        fig_idx += 1


    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('Percent of SVs', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
        ax.set_xlabel('')
        ax.set_xticks(np.arange(2))
        ax.set_xticklabels(['WGS', 'HighConf'], fontsize=13)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

        ax.legend(title='')


    fig.tight_layout()

    fig.savefig(f'{Figure4}/acc_concordant_pcrt_count.pdf')
    plt.show()


def wgs_highconf_acc_insdel(workdir, aligners, datasets):

    bpstd_pcrt = []
    highconf_bpstd_pcrt = []
    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for col_idx, aligner in enumerate(aligners):
            for assembler in plat_assemblers[plat]:
                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:
                        ins_bpstd_dict = {'Acc': 0, 'Inacc': 0}
                        del_bpstd_dict = {'Acc': 0, 'Inacc': 0}

                        highconf_ins_bpstd_dict = {'Acc': 0, 'Inacc': 0}
                        highconf_del_bpstd_dict = {'Acc': 0, 'Inacc': 0}

                        ins_total = 0
                        del_total = 0

                        highconf_ins_total = 0
                        highconf_del_total = 0

                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(
                                float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if start_std <= 10 and end_std <= 10:
                                if svtype == 'INS':
                                    ins_total += 1
                                    ins_bpstd_dict['Acc'] += 1
                                if svtype == 'DEL':
                                    del_total += 1
                                    del_bpstd_dict['Acc'] += 1
                            else:
                                if svtype == 'INS':
                                    ins_total += 1
                                if svtype == 'DEL':
                                    del_total += 1

                        bpstd_pcrt.append(('INS', ins_bpstd_dict['Acc'] / ins_total * 100, assembler, aligner, caller, asm_caller))
                        bpstd_pcrt.append(('DEL', del_bpstd_dict['Acc'] / del_total * 100, assembler, aligner, caller, asm_caller))

                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered_in_highconf/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(
                                float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if start_std <= 10 and end_std <= 10:
                                if svtype == 'INS':
                                    highconf_ins_total += 1
                                    highconf_ins_bpstd_dict['Acc'] += 1
                                if svtype == 'DEL':
                                    highconf_del_total += 1
                                    highconf_del_bpstd_dict['Acc'] += 1
                            else:
                                if svtype == 'INS':
                                    highconf_ins_total += 1
                                if svtype == 'DEL':
                                    highconf_del_total += 1

                        highconf_bpstd_pcrt.append(('INS', highconf_ins_bpstd_dict['Acc'] / highconf_ins_total * 100, assembler, aligner, caller, asm_caller))
                        highconf_bpstd_pcrt.append(('DEL', highconf_del_bpstd_dict['Acc'] / highconf_del_total * 100, assembler, aligner, caller, asm_caller))

    df_bpstd = pd.DataFrame(bpstd_pcrt, columns=['svtype', 'pcrt', 'assembler', 'aligner', 'caller', 'asmcaller'])
    df_highconf_bpstd = pd.DataFrame(highconf_bpstd_pcrt, columns=['svtype', 'pcrt', 'assembler', 'aligner', 'caller', 'asmcaller'])
    aligner_color = [ALIGNERCOLOR[ele] for ele in aligners]

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    sns.boxplot(data=df_bpstd[df_bpstd['svtype'] == 'INS'], x='assembler', y='pcrt', hue='aligner', hue_order=aligners, palette=aligner_color, ax=axes[0])
    axes[0].set_title('INS', fontsize=13)
    sns.boxplot(data=df_bpstd[df_bpstd['svtype'] == 'DEL'], x='assembler', y='pcrt', hue='aligner', hue_order=aligners, palette=aligner_color, ax=axes[1])
    axes[1].set_title('DEL', fontsize=13)
    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of SVs', fontsize=13)
        else:
            ax.set_ylabel('')
        ax.legend('')
        ax.set_xlabel('')
        ax.set_xticks(np.arange(3))
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=12)

        ax.set_ylim(25, 100)
        ax.set_yticks(np.linspace(25, 100, 4))
        ax.set_yticklabels([int(val) for val in np.linspace(25, 100, 4)], fontsize=12)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()

    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    sns.boxplot(data=df_highconf_bpstd[df_highconf_bpstd['svtype'] == 'INS'], x='assembler', y='pcrt', hue='aligner', hue_order=aligners,
                palette=aligner_color, ax=axes[0])
    axes[0].set_title('INS', fontsize=13)
    sns.boxplot(data=df_highconf_bpstd[df_highconf_bpstd['svtype'] == 'DEL'], x='assembler', y='pcrt', hue='aligner', hue_order=aligners,
                palette=aligner_color, ax=axes[1])

    axes[1].set_title('DEL', fontsize=13)
    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of SVs', fontsize=13)
        else:
            ax.set_ylabel('')
        ax.legend('')
        ax.set_xlabel('')
        ax.set_xticks(np.arange(3))
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=12)

        ax.set_ylim(25, 100)
        ax.set_yticks(np.linspace(25, 100, 4))
        ax.set_yticklabels([int(val) for val in np.linspace(25, 100, 4)], fontsize=12)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=1.5)

    fig.tight_layout()
    fig1.tight_layout()

    fig.savefig(f'{Figure4}/wgs_acc_insdel_pcrt_aligner_assembler.pdf')
    fig1.savefig(f'{Figure4}/highconf_acc_insdel_pcrt_aligner_assembler.pdf')
    plt.show()


def inacc_concordant_svsize(workdir, aligners, datasets):

    inacc_bpstd_list = {'HiFi': [], 'ONT': []}
    # highconf_inacc_bpstd_list = {'HiFi': [], 'ONT': []}

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for col_idx, aligner in enumerate(aligners):
            for assembler in plat_assemblers[plat]:
                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:
                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(
                                float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if start_std <= 10 and end_std <= 10:
                                continue

                            if svtype == 'INS' and sv_region == 'Simple Repeats':
                                inacc_bpstd_list[plat].append((svlen, 'INS', sv_region, aligner, assembler))
                            if svtype == 'DEL' and sv_region == 'Simple Repeats':
                                inacc_bpstd_list[plat].append((svlen, 'DEL', sv_region, aligner, assembler))

                        # matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered_in_highconf/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        # df_matched = pd.read_csv(matched_file, sep='\t', header=[0])
                        # for idx, row in df_matched.iterrows():
                        #     start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(
                        #         float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                        #
                        #     if start_std <= 10 and end_std <= 10:
                        #         continue
                        #
                        #     if svtype == 'INS' and sv_region == 'Simple Repeats':
                        #         highconf_inacc_bpstd_list[plat].append((svlen, 'INS', sv_region, aligner, assembler))
                        #
                        #     if svtype == 'DEL' and sv_region == 'Simple Repeats':
                        #         highconf_inacc_bpstd_list[plat].append((svlen, 'DEL', sv_region, aligner, assembler))

    df_inacc_hifi = pd.DataFrame(inacc_bpstd_list['HiFi'], columns=['svlen', 'svtype', 'region', 'aligner', 'assembler'])
    df_inacc_ont = pd.DataFrame(inacc_bpstd_list['ONT'], columns=['svlen', 'svtype', 'region', 'aligner', 'assembler'])


    svtype_order = ['INS', 'DEL']
    # svtype_color = [SVTYPECOLORS[svtype] for svtype in svtype_order]
    # svtype_legends = [Patch(facecolor=SVTYPECOLORS['INS'], label='INS'), Patch(facecolor=SVTYPECOLORS['DEL'], label='DEL')]

    fig, axes = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(7, 4))
    sns.histplot(data=df_inacc_hifi[df_inacc_hifi['svlen'] < 1000], x='svlen', hue='aligner', bins=100, kde=True, ax=axes[0])

    axes[0].set_ylabel('HiFi (x$10^3$)', fontsize=13)

    sns.histplot(data=df_inacc_ont[df_inacc_ont['svlen'] < 1000], x='svlen', hue='aligner', bins=100, kde=True, ax=axes[1])
    axes[1].set_ylabel('ONT (x$10^3$)', fontsize=13)

    for ax in axes:
        ax.set_xlabel('')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend()

        ax.set_ylim(0, 40000)
        ax.set_yticks(np.linspace(0, 40000, 3))
        ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 40, 3)], fontsize=12)

        ax.set_xticks(np.linspace(0, 1000, 6))
        ax.set_xticklabels([int(val) for val in np.linspace(0, 1000, 6)], fontsize=12)

    fig.tight_layout()


    fig.savefig(f'{Figure4}/wgs_inacc_concordant_simrep_svsize.pdf')
    plt.show()

def wgs_highconf_acc_concordant_factors(workdir, aligners, datasets):

    bpstd_pcrt_list = []

    shift_labels = ['0,10', '10,50', '50,100', '>100']
    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for dataset in datasets:
        ins_bpstd_dict = {shift: 0 for shift in shift_labels}
        del_bpstd_dict = {shift: 0 for shift in shift_labels}

        highconf_ins_bpstd_dict = {shift: 0 for shift in shift_labels}
        highconf_del_bpstd_dict = {shift: 0 for shift in shift_labels}

        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        ins_total = 0
        del_total = 0

        highconf_ins_total = 0
        highconf_del_total = 0

        for col_idx, aligner in enumerate(aligners):
            for assembler in plat_assemblers[plat]:
                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:
                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(
                                float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if start_std <= 10 and end_std <= 10:
                                if svtype == 'INS':
                                    ins_total += 1
                                    ins_bpstd_dict['0,10'] += 1
                                if svtype == 'DEL':
                                    del_total += 1
                                    del_bpstd_dict['0,10'] += 1
                            else:
                                if svtype == 'INS':
                                    ins_total += 1
                                if svtype == 'DEL':
                                    del_total += 1

                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered_in_highconf/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(
                                float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if start_std <= 10 and end_std <= 10:
                                if svtype == 'INS':
                                    highconf_ins_total += 1
                                    highconf_ins_bpstd_dict['0,10'] += 1
                                if svtype == 'DEL':
                                    highconf_del_total += 1
                                    highconf_del_bpstd_dict['0,10'] += 1
                            else:
                                if svtype == 'INS':
                                    highconf_ins_total += 1
                                if svtype == 'DEL':
                                    highconf_del_total += 1

                bpstd_pcrt_list.append((ins_bpstd_dict['0,10'] / ins_total * 100, plat, aligner, assembler, 'INS', 'WGS'))
                bpstd_pcrt_list.append((del_bpstd_dict['0,10'] / del_total * 100, plat, aligner, assembler, 'DEL', 'WGS'))

                bpstd_pcrt_list.append((highconf_ins_bpstd_dict['0,10'] / highconf_ins_total * 100, plat, aligner, assembler, 'INS', 'HighConf'))
                bpstd_pcrt_list.append((highconf_del_bpstd_dict['0,10'] / highconf_del_total * 100, plat, aligner, assembler, 'DEL', 'HighConf'))

    hue_order = ['INS', 'DEL']
    hue_color = [SVTYPECOLORS[val] for val in hue_order]
    legends = [Patch(facecolor=SVTYPECOLORS[val], label=val) for val in hue_order]
    assemblers = ['flye', 'shasta', 'hifiasm']
    assem_colors = [ASSMBLERCOLOR[ele] for ele in assemblers]

    df_bpstd = pd.DataFrame(bpstd_pcrt_list, columns=['pcrt', 'plat', 'aligner', 'assembler', 'svtype', 'region'])

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    sns.boxplot(data=df_bpstd[df_bpstd['region'] == 'WGS'], x='plat', y='pcrt', hue='svtype', hue_order=hue_order,
                palette=hue_color, ax=axes[0])
    sns.stripplot(data=df_bpstd[df_bpstd['region'] == 'WGS'], x='plat', y='pcrt', hue_order=assemblers, palette=assem_colors,
                  hue='assembler', size=6, ax=axes[0])

    sns.boxplot(data=df_bpstd[df_bpstd['region'] == 'HighConf'], x='plat', y='pcrt', hue='svtype', hue_order=hue_order,
                palette=hue_color, ax=axes[1])
    sns.stripplot(data=df_bpstd[df_bpstd['region'] == 'HighConf'], x='plat', y='pcrt', hue_order=assemblers, palette=assem_colors,
                  hue='assembler', size=6, ax=axes[1])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('Percent of SVs', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_ylim(40, 100)
        ax.set_yticks(np.linspace(40, 100, 4))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 4)], fontsize=12)

        ax.legend(handles=legends)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel('')

    fig.tight_layout()

    aligner_color = [ALIGNERCOLOR[aligner] for aligner in aligners]

    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 4))
    sns.boxplot(data=df_bpstd[df_bpstd['region']=='WGS'], x='plat', y='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, ax=axes[0])
    sns.stripplot(data=df_bpstd[df_bpstd['region'] == 'WGS'], x='plat', y='pcrt', hue='aligner', hue_order=aligners, size=6,
                palette=aligner_color, ax=axes[0])

    sns.boxplot(data=df_bpstd[df_bpstd['region'] == 'HighConf'], x='plat', y='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, ax=axes[1])
    sns.stripplot(data=df_bpstd[df_bpstd['region'] == 'HighConf'], x='plat', y='pcrt', hue='aligner', hue_order=aligners, size=6,
                palette=aligner_color, ax=axes[1])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('Percent of SVs', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_ylim(40, 100)
        ax.set_yticks(np.linspace(40, 100, 4))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 4)], fontsize=12)

        ax.legend(handles=legends)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel('')

    fig1.tight_layout()


    fig.savefig(f'{Figure4}/wgs_highconf_acc_concordant_byassembler.pdf')
    fig1.savefig(f'{Figure4}/wgs_highconf_acc_concordant_byaligners.pdf')

    plt.show()

def highconf_stra_concordant_bp_overview(workdir, aligners, datasets):
    shift_labels = ['0', '0,10', '10,50', '>50']
    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    ins_shift_dict = {'Left': [[] for i in range(len(shift_labels))], 'Right': [[] for i in range(len(shift_labels))]}
    del_shift_dict = {'Left': [[] for i in range(len(shift_labels))], 'Right': [[] for i in range(len(shift_labels))]}

    for col_idx, aligner in enumerate(aligners):
        for dataset in datasets:
            ins_left_bpstd_dict = {shift: 0 for shift in shift_labels}
            ins_right_bpstd_dict = {shift: 0 for shift in shift_labels}
            del_left_bpstd_dict = {shift: 0 for shift in shift_labels}
            del_right_bpstd_dict = {shift: 0 for shift in shift_labels}

            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            ins_total = 0
            del_total = 0

            for caller in CALLERS:
                for asm_caller in ASMCALLERS:
                    for assembler in plat_assemblers[plat]:
                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered_in_highconf/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if svtype == 'INS':
                                ins_total += 1
                                if start_std == 0:
                                    ins_left_bpstd_dict['0'] += 1
                                elif start_std <= 10:
                                    ins_left_bpstd_dict['0,10'] += 1
                                elif start_std <= 50:
                                    ins_left_bpstd_dict['10,50'] += 1
                                else:
                                    ins_left_bpstd_dict['>50'] += 1

                                if end_std == 0:
                                    ins_right_bpstd_dict['0'] += 1
                                elif end_std <= 10:
                                    ins_right_bpstd_dict['0,10'] += 1
                                elif start_std <= 50:
                                    ins_right_bpstd_dict['10,50'] += 1
                                else:
                                    ins_right_bpstd_dict['>50'] += 1

                            if svtype == 'DEL':
                                del_total += 1
                                if start_std == 0:
                                    del_left_bpstd_dict['0'] += 1
                                elif start_std <= 10:
                                    del_left_bpstd_dict['0,10'] += 1
                                elif start_std <= 50:
                                    del_left_bpstd_dict['10,50'] += 1
                                else:
                                    del_left_bpstd_dict['>50'] += 1

                                if end_std == 0:
                                    del_right_bpstd_dict['0'] += 1
                                elif end_std <= 10:
                                    del_right_bpstd_dict['0,10'] += 1
                                elif start_std <= 50:
                                    del_right_bpstd_dict['10,50'] += 1
                                else:
                                    del_right_bpstd_dict['>50'] += 1

                for shift_idx, shift in enumerate(shift_labels):
                    # leftbp_pcrt_list.append((ins_left_bpstd_dict[shift] / ins_total * 100, plat, 'INS', aligner, shift))
                    # leftbp_pcrt_list.append((del_left_bpstd_dict[shift] / del_total * 100, plat, 'DEL', aligner, shift))
                    #
                    # rightbp_pcrt_list.append((ins_right_bpstd_dict[shift] / ins_total * 100, plat, 'INS', aligner, shift))
                    # rightbp_pcrt_list.append((del_right_bpstd_dict[shift] / del_total * 100, plat, 'DEL', aligner, shift))

                    ins_shift_dict['Left'][shift_idx].append(ins_left_bpstd_dict[shift])
                    ins_shift_dict['Right'][shift_idx].append(ins_right_bpstd_dict[shift])

                    del_shift_dict['Left'][shift_idx].append(del_left_bpstd_dict[shift])
                    del_shift_dict['Right'][shift_idx].append(del_right_bpstd_dict[shift])

    # hue_order = ['INS', 'DEL']
    # hue_color = [SVTYPECOLORS[val] for val in hue_order]
    # svtype_legends = [Patch(facecolor=SVTYPECOLORS[val], label=val) for val in ['INS', 'DEL']]
    #
    #
    # df_leftbp = pd.DataFrame(leftbp_pcrt_list, columns=['pcrt', 'plat', 'svtype', 'aligner', 'shift'])
    # df_rightbp = pd.DataFrame(rightbp_pcrt_list, columns=['pcrt', 'plat', 'svtype', 'aligner', 'shift'])
    #
    # fig, axes = plt.subplots(1, 2,  figsize=(6, 4))
    #
    # sns.boxplot(data=df_leftbp, y='shift', x='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, ax=axes[0])
    # axes[0].set_ylabel('Left breakpoint std. (bp)', fontsize=13)
    # sns.boxplot(data=df_rightbp, y='shift', x='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, ax=axes[1])
    # axes[1].set_ylabel('Right breakpoint std. (bp)', fontsize=13)
    #
    # for i, ax in enumerate(axes):
    #
    #     ax.set_xlim(0, 80)
    #     ax.set_xticks(np.linspace(0, 80, 3))
    #     ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 80, 3)], fontsize=12)
    #     ax.set_xlabel('Percent of SVs', fontsize=13)
    #
    #     ax.legend(handles=svtype_legends)
    #
    #     ax.set_yticks(np.arange(3))
    #     ax.set_yticklabels(['0', '0~10', '10~50'], fontsize=12)

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    shift_colors = [SHIFTCOLORDICT[shift] for shift in shift_labels]
    size = 0.3
    for key, values in ins_shift_dict.items():
        shift_avg = []
        for i, counts in enumerate(values):
            shift_avg.append(np.mean(counts))

        if key == 'Left':
            ax.pie(shift_avg, autopct='%1.1f%%', radius=1, colors=shift_colors,
                   startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))
        else:
            ax.pie(shift_avg, autopct='%1.1f%%', radius=1 - size, colors=shift_colors,
                   startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))

    fig.tight_layout()

    fig1, ax = plt.subplots(1, 1, figsize=(4, 4))
    shift_colors = [SHIFTCOLORDICT[shift] for shift in shift_labels]
    size = 0.3
    for key, values in del_shift_dict.items():
        shift_avg = []
        for i, counts in enumerate(values):
            shift_avg.append(np.mean(counts))

        if key == 'Left':
            ax.pie(shift_avg, autopct='%1.1f%%', radius=1, colors=shift_colors,
                   startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))
        else:
            ax.pie(shift_avg, autopct='%1.1f%%', radius=1 - size, colors=shift_colors,
                   startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))
    fig1.tight_layout()

    fig.savefig(f'{Figure4}/highconf_ins_bpstd_avg_pieplot.pdf')
    fig1.savefig(f'{Figure4}/highconf_del_bpstd_avg_pieplot.pdf')

    plt.show()

def highconf_stra_concordant_insdel_bpstd(workdir, aligners, datasets):

    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']
    bpstd_list = []
    inacc_bpstd_list = []

    shift_labels = ['0,10', '10,50', '50,100', '>100']
    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    for dataset in datasets:
        ins_bpstd_dict = {shift: 0 for shift in shift_labels}
        del_bpstd_dict = {shift: 0 for shift in shift_labels}

        inacc_ins = {region: 0 for region in regions}
        inacc_del = {region: 0 for region in regions}

        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        ins_total = 0
        del_total = 0

        for col_idx, aligner in enumerate(aligners):
            for assembler in plat_assemblers[plat]:

                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:
                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered_in_highconf/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.concordants.info.tsv'
                        df_matched = pd.read_csv(matched_file, sep='\t', header=[0])

                        for idx, row in df_matched.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                            if start_std <= 10 and end_std <= 10:
                                if svtype == 'INS':
                                    ins_total += 1
                                    ins_bpstd_dict['0,10'] += 1
                                if svtype == 'DEL':
                                    del_total += 1
                                    del_bpstd_dict['0,10'] += 1
                            else:
                                if svtype == 'INS':
                                    ins_total += 1
                                    inacc_ins[sv_region] += 1
                                if svtype == 'DEL':
                                    del_total += 1
                                    inacc_del[sv_region] += 1

                bpstd_list.append((ins_bpstd_dict['0,10'] / ins_total * 100, plat, aligner, assembler, 'INS'))
                bpstd_list.append((del_bpstd_dict['0,10'] / del_total * 100, plat, aligner, assembler, 'DEL'))
                for region in regions:
                    inacc_bpstd_list.append((inacc_ins[region] / (ins_total - ins_bpstd_dict['0,10']) * 100, plat, aligner, assembler, region, 'INS'))
                    inacc_bpstd_list.append((inacc_del[region] / (del_total - del_bpstd_dict['0,10']) * 100, plat, aligner, assembler,region, 'DEL'))

    hue_order = ['INS', 'DEL']
    hue_color = [SVTYPECOLORS[val] for val in hue_order]
    legends = [Patch(facecolor=SVTYPECOLORS[val], label=val) for val in hue_order]

    df_bpstd = pd.DataFrame(bpstd_list, columns=['pcrt', 'plat', 'aligner', 'assembler', 'svtype'])
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sns.lineplot(data=df_bpstd, x='aligner', y='pcrt', hue='svtype', hue_order=hue_order,
                palette=hue_color, style='svtype', lw=2, markers=True, ms=10, ax=ax)
    ax.set_ylabel('Percent of SVs', fontsize=13)
    ax.set_ylim(40, 100)
    ax.set_yticks(np.linspace(40, 100, 4))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 4)], fontsize=12)
    ax.set_xticks(np.arange(len(aligners)))
    ax.set_xticklabels(aligners, fontsize=13)

    ax.legend(handles=legends)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('')
    fig.tight_layout()

    fig1, ax = plt.subplots(1, 1, figsize=(3, 4))
    sns.barplot(data=df_bpstd, x='plat', y='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, capsize=.2, ax=ax)
    ax.set_ylabel('Percent of SVs', fontsize=13)
    ax.set_ylim(40, 100)
    ax.set_yticks(np.linspace(40, 100, 4))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 4)])
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=13)
    ax.legend(handles=legends)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('')
    fig1.tight_layout()

    df_inacc_bpstd = pd.DataFrame(inacc_bpstd_list, columns=['pcrt','plat','aligner', 'assembler', 'region','svtype'])
    fig2, ax = plt.subplots(1, 1, figsize=(4, 4))
    sns.barplot(data=df_inacc_bpstd, y='region', x='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, capsize=.2, ax=ax)
    ax.set_ylabel('')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)
    ax.set_xlim(0, 100)
    ax.set_xticks(np.linspace(0, 100, 5))
    ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
    ax.set_xlabel('Percent of SVs')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(handles=legends)

    fig2.tight_layout()

    plt.show()
    fig.savefig(f'{Figure4}/highconf_strategy_matched_insdel_bpstd_byaligner.pdf')
    fig1.savefig(f'{Figure4}/highconf_strategy_matched_insdel_bpstd_byplat.pdf')
    fig2.savefig(f'{Figure4}/highconf_strategy_matched_insdel_inacc_bpstd_byregions.pdf')





def highconf_stra_compare_overview(workdir, aligners, datasets):

    highconf_stra_counts = []
    highconf_stra_pcrts = []
    highconf_merged_counts = []

    stra_counts = []
    stra_pcrts = []

    diff_region_stra_pcrts = []

    region = 'highconf'
    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    highconf_aligner_stra_pcrt = [[] for i in range(len(aligners))]
    aligner_stra_pcrt = [[] for i in range(len(aligners))]
    aligner_assembler_pcrt = {}

    for col_idx, aligner in enumerate(aligners):
        highconf_stra_matched = [[0 for i in range(len(datasets))] for j in range(3)]
        stra_matched = [[0 for i in range(len(datasets))] for j in range(3)]
        for dataset in datasets:
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            for caller in CALLERS:
                for asm_caller in ASMCALLERS:
                    for assembler in plat_assemblers[plat]:
                        merged_vcf = f'{workdir}/{plat}/{aligner}_{dataset}/filtered_in_{region}/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.jasmine.merged.vcf'
                        suppvec_dict, merged_total = get_survivor_suppvec(merged_vcf)
                        dataset_idx = datasets.index(dataset)
                        highconf_merged_counts.append((plat, merged_total, PLATMAP[dataset], aligner, assembler,
                                                       f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}', TOOLMAP[asm_caller]))
                        for suppvec, count in suppvec_dict.items():
                            if suppvec == '10':
                                highconf_stra_matched[0][dataset_idx] += count
                                highconf_stra_counts.append((plat, count, 'Assembly-Unique', PLATMAP[dataset],
                                                    f'minimap2-{aligner}', f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))
                            elif suppvec == '01':
                                highconf_stra_matched[2][dataset_idx] += count
                                highconf_stra_counts.append((plat, count, 'Read-Unique', PLATMAP[dataset], f'minimap2-{aligner}',
                                                    f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))
                            else:
                                highconf_stra_matched[1][dataset_idx] += count
                                highconf_stra_counts.append((plat, count, 'Concordant', PLATMAP[dataset], f'minimap2-{aligner}',
                                                    f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))

                        highconf_stra_pcrts.append((plat, suppvec_dict['11'] / merged_total * 100, PLATMAP[dataset],
                                           aligner, f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))


                        aln_asm_key = f'{aligner}-{assembler}'
                        if aln_asm_key in aligner_assembler_pcrt:
                            aligner_assembler_pcrt[aln_asm_key].append(suppvec_dict['11'] / merged_total * 100)
                        else:
                            aligner_assembler_pcrt[aln_asm_key] = [suppvec_dict['11'] / merged_total * 100]

                        diff_region_stra_pcrts.append((f'{plat}-HighConf', suppvec_dict['11'] / merged_total * 100, PLATMAP[dataset],
                                           aligner, f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))

                        merged_vcf = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.jasmine.merged.vcf'
                        suppvec_dict, merged_total = get_survivor_suppvec(merged_vcf)
                        dataset_idx = datasets.index(dataset)

                        for suppvec, count in suppvec_dict.items():
                            if suppvec == '10':
                                stra_matched[0][dataset_idx] += count
                                stra_counts.append((plat, count, 'Assembly-Unique', PLATMAP[dataset], aligner,
                                                    f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))
                            elif suppvec == '01':
                                stra_matched[2][dataset_idx] += count
                                stra_counts.append((plat, count, 'Read-Unique', PLATMAP[dataset], aligner,
                                                    f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))
                            else:
                                stra_matched[1][dataset_idx] += count
                                stra_counts.append((plat, count, 'Concordant', PLATMAP[dataset], aligner,
                                                    f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))

                        stra_pcrts.append((plat, suppvec_dict['11'] / merged_total * 100, PLATMAP[dataset],
                                           aligner, assembler, f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))

                        diff_region_stra_pcrts.append((f'{plat}-WGS', suppvec_dict['11'] / merged_total * 100, PLATMAP[dataset],
                                           aligner, f'{TOOLMAP[caller]}-{TOOLMAP[asm_caller]}'))

        for j in range(len(datasets)):
            highconf_matched_pcrt = highconf_stra_matched[1][j] / (highconf_stra_matched[0][j] + highconf_stra_matched[1][j] + highconf_stra_matched[2][j])
            highconf_aligner_stra_pcrt[col_idx].append(highconf_matched_pcrt * 100)

            matched_pcrt = stra_matched[1][j] / (stra_matched[0][j] + stra_matched[1][j] + stra_matched[2][j])
            aligner_stra_pcrt[col_idx].append(matched_pcrt * 100)

    fig, axes = plt.subplots(1, 2, sharey='row', sharex='col', figsize=(5, 4))
    df_highconf_merged_counts = pd.DataFrame(highconf_merged_counts, columns=['plat', 'count', 'dataset', 'aligner', 'assembler', 'callers', 'asmcaller'])

    sns.stripplot(data=df_highconf_merged_counts, y='count', x='plat', hue='dataset', palette="Dark2", ax=axes[0])
    box = sns.boxplot(data=df_highconf_merged_counts, y='count', x='plat', ax=axes[0])
    for i in [0, 1]:
        mybox = box.artists[i]
        mybox.set_facecolor('white')
        mybox.set_edgecolor('black')

    sns.stripplot(data=df_highconf_merged_counts, x='plat', y='count', hue='asmcaller', hue_order=['PAV', 'SVIM-asm'],
                  palette=[TOOLCOLORS['pav'], TOOLCOLORS['svimasm']], size=6, ax=axes[1])
    box = sns.boxplot(data=df_highconf_merged_counts, x='plat', y='count', ax=axes[1])
    for i in [0, 1]:
        mybox = box.artists[i]
        mybox.set_facecolor('white')
        mybox.set_edgecolor('black')

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('Number of nonredundant SVs (x$10^3$)', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_ylim(12000, 22000, 3)
        ax.set_yticks(np.linspace(12000, 22000, 3))
        ax.set_yticklabels([int(val) for val in np.linspace(12, 22, 3)], fontsize=12)
        ax.set_xlabel('')
        ax.set_xticks(np.arange(2))
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=13)
        ax.legend(title='')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig.tight_layout()

    region_legends = [Patch(facecolor='#cbd5e7', label='HighConf'), Patch(facecolor='#f6cae4', label='Whole Genome')]

    fig1, ax = plt.subplots(1, 1, figsize=(3, 4))
    df_region_stra_pcrts = pd.DataFrame(diff_region_stra_pcrts, columns=['plat', 'pcrt', 'dataset', 'aligner', 'callers'])
    box = sns.boxplot(data=df_region_stra_pcrts, y='pcrt', x='plat', ax=ax)
    # sns.stripplot(data=df_region_stra_pcrts, y='pcrt', x='plat', hue='dataset', palette="Dark2", ax=ax)

    for i in np.arange(4):
        mybox = box.artists[i]
        mybox.set_edgecolor('black')
        if i in [0, 2]:
            mybox.set_facecolor('#cbd5e7')
        else:
            mybox.set_facecolor('#f6cae4')

    ax.set_ylabel('Percent of concordant SVs', fontsize=13)
    ax.set_ylim(20, 100, 5)
    ax.set_yticks(np.linspace(20, 100, 5))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(20, 100, 5)], fontsize=12)
    ax.legend(handles=region_legends)

    ax.set_xticks([0.5, 2.5])
    ax.set_xticklabels(['HiFi', 'ONT'], fontsize=13)
    ax.set_xlabel('')


    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig1.tight_layout()

    fig2, ax = plt.subplots(1, 1, figsize=(5, 4))
    # barwidth = 0.3
    # r1 = np.arange(len(datasets))
    # r2 = [r + barwidth + 0.05 for r in r1]


    xticks = np.arange(len(datasets))
    for i in range(len(aligners)):
        # ax.bar(r1, highconf_aligner_stra_pcrt[i], label=f'{aligners[i]}', width=barwidth, color='#cbd5e7', edgecolor='black')
        # ax.bar(r2, aligner_stra_pcrt[i], label=f'{aligners[i]}', width=barwidth, color='#f6cae4', edgecolor='black')
        ax.scatter(highconf_aligner_stra_pcrt[i], aligner_stra_pcrt[i], color=ALIGNERCOLOR[aligners[i]],
                   marker=ALIGNERMK[aligners[i]], s=90, label=aligners[i])


    ax.legend()
    ax.set_ylim(30, 90)
    ax.set_yticks(np.linspace(30, 90, 4))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(30, 90, 4)], fontsize=12)
    ax.set_ylabel('Percent of SVs at WGS', fontsize=13)

    # ax.set_xticks(xticks)
    # ax.set_xticklabels([PLATMAP[val] for val in datasets], rotation=90, fontsize=13)

    ax.set_xlim(30, 90)
    ax.set_xticks(np.linspace(30, 90, 4))
    ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(30, 90, 4)], fontsize=12)
    ax.set_xlabel('Percent of SVs at HighConf', fontsize=13)

    ax.plot([30, 90], [30, 90], ls='--', color='#c96150')

    fig2.tight_layout()

    plt.show()
    fig.savefig(f'{Figure4}/highconf_strategy_matched_overview.pdf')
    fig1.savefig(f'{Figure4}/wgs_highconf_strategy_matched_pcrt_of_plats.pdf')
    fig2.savefig(f'{Figure4}/wgs_highconf_strategy_matched_pcrt_of_aligners.pdf')

def plot_stra_concordant_features_at_regions(workdir, aligners, datasets):
    ins_bpstd_list = []
    del_bpstd_list = []
    shift_labels = ['0,10', '10,50', '50,100', '>100']

    region = 'highconf'


    for col_idx, aligner in enumerate(aligners):

        for dataset in datasets:
            ins_bpstd_dict = {shift: 0 for shift in shift_labels}
            del_bpstd_dict = {shift: 0 for shift in shift_labels}
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            matched_total = 0
            for caller in CALLERS:
                for asm_caller in ASMCALLERS:
                    matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered_in_{region}/comstra/{caller}-{asm_caller}.{aligner}-minimap2.{dataset}.caller_concordants.info.tsv'
                    df_matched = pd.read_csv(matched_file, sep='\t', header=[0])
                    matched_total += len(df_matched)

                    for idx, row in df_matched.iterrows():
                        start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(
                            float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                        if start_std <= 10 and end_std <= 10:
                            if svtype == 'INS':
                                ins_bpstd_dict['0,10'] += 1
                            if svtype == 'DEL':
                                del_bpstd_dict['0,10'] += 1
                        elif start_std <= 50 and end_std <= 50:
                            if svtype == 'INS':
                                ins_bpstd_dict['10,50'] += 1
                            if svtype == 'DEL':
                                del_bpstd_dict['10,50'] += 1
                        elif start_std <= 100 and end_std <= 100:
                            if svtype == 'INS':
                                ins_bpstd_dict['50,100'] += 1
                            if svtype == 'DEL':
                                del_bpstd_dict['50,100'] += 1
                        else:
                            if svtype == 'INS':
                                ins_bpstd_dict['>100'] += 1
                            if svlen == 'DEL':
                                del_bpstd_dict['>100'] += 1

            for shift_label in shift_labels:
                ins_bpstd_list.append((shift_label, ins_bpstd_dict[shift_label] / matched_total, plat, aligner))
                del_bpstd_list.append((shift_label, del_bpstd_dict[shift_label] / matched_total, plat, aligner))

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 4))
    df_ins_bpstd = pd.DataFrame(ins_bpstd_list, columns=['label', 'pcrt', 'plat', 'aligner'])
    df_del_bpstd = pd.DataFrame(del_bpstd_list, columns=['label', 'pcrt', 'plat', 'aligner'])
    sns.boxplot(data=df_ins_bpstd, x='label', y='pcrt', hue='aligner', ax=axes[0])
    axes[0].set_title('INS', fontsize=13)
    sns.boxplot(data=df_del_bpstd, x='label', y='pcrt', hue='aligner', ax=axes[1])
    axes[1].set_title('DEL', fontsize=13)
    for ax in axes:
        ax.legend(title='')
        ax.set_ylabel('')
        ax.set_xlabel('')

    fig.tight_layout()

    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 4))
    sns.violinplot(data=df_ins_bpstd, x='label', y='pcrt', hue='plat', ax=axes[0])
    axes[0].set_title('INS', fontsize=13)
    sns.violinplot(data=df_del_bpstd, x='label', y='pcrt', hue='plat', ax=axes[1])
    axes[1].set_title('DEL', fontsize=13)

    for ax in axes:
        ax.legend(title='')
        ax.set_ylabel('')
        ax.set_xlabel('')

    fig1.tight_layout()

    plt.show()

    fig.savefig(f'{Figure4}/highconf_strategy_matched_insdel_bpstd_byaligner.pdf')
    fig1.savefig(f'{Figure4}/highconf_strategy_matched_insdel_bpstd_byplat.pdf')






'''
def overview_stra_compare_insdel(workdir, aligners, datasets):

    svtypes = ['ins', 'del']

    for svtype in svtypes:

        # fig, axes = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(11, 4))
        xticks = np.arange(len(datasets))

        aligner_stra_pcrt = [[] for i in range(len(aligners))]

        for col_idx, aligner in enumerate(aligners):
            stra_matched = [[0 for i in range(len(datasets))] for j in range(3)]
            for dataset in datasets:
                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:
                        merged_vcf = f'{workdir}/HiFi/{aligner}_{dataset}/filtered/comstra_insdel/{caller}-{asm_caller}.{aligner}-minimap2.{dataset}.{svtype}.jasmine.merged.vcf'
                        if 'ont' in dataset:
                            merged_vcf = f'{workdir}/ONT/{aligner}_{dataset}/filtered/comstra_insdel/{caller}-{asm_caller}.{aligner}-minimap2.{dataset}.{svtype}.jasmine.merged.vcf'

                        suppvec_dict, merged_total = get_survivor_suppvec(merged_vcf)
                        dataset_idx = datasets.index(dataset)

                        for suppvec, count in suppvec_dict.items():
                            if suppvec == '10':
                                stra_matched[0][dataset_idx] += count
                            elif suppvec == '01':
                                stra_matched[2][dataset_idx] += count
                            else:
                                stra_matched[1][dataset_idx] += count

            for j in range(len(datasets)):
                stra_matched_pcrt = stra_matched[1][j] / (stra_matched[0][j] + stra_matched[1][j] + stra_matched[2][j])
                aligner_stra_pcrt[col_idx].append(stra_matched_pcrt)

            # labels = ['Assembly-Unique', 'Concordant', 'Read-Unique']
            # this_ax = axes[col_idx]
            # for i in range(3):
            #     this_data = stra_matched[i]
            #     if i == 0:
            #         this_ax.bar(xticks, this_data, label=labels[i])
            #     else:
            #         bottoms = []
            #         for j in range(0, i):
            #             bottoms.append(stra_matched[j])
            #         bottom_sum = [sum(x) for x in zip(*bottoms)]
            #         this_ax.bar(xticks, this_data, bottom=bottom_sum, label=labels[i])
            #
            # this_ax.set_title(f'minimap2-{aligner}')
            # this_ax.set_ylim(0, 350000, 6)
            # this_ax.set_yticks(np.linspace(0, 350000, 6))
            # this_ax.set_yticklabels([int(val) for val in np.linspace(0, 35, 6)])
            # if col_idx == 0:
            #     this_ax.set_ylabel('# of SVs ($10^4$)')
            #
            # this_ax.set_xticks(xticks)
            # this_ax.set_xticklabels([PLATMAP[val] for val in datasets], rotation=90)
            # this_ax.legend()


        fig1, ax1 = plt.subplots(1, 1, figsize=(5, 4))
        for i in range(len(aligners)):
            ax1.plot(xticks, aligner_stra_pcrt[i], label=f'minimap2-{aligners[i]}',lw=2, marker='o')

        ax1.legend()
        # if svtype == 'ins':
        ax1.set_ylim(0.3, 0.8)
        ax1.set_yticks(np.linspace(0.3, 0.8, 6))
        ax1.set_yticklabels([f'{int(val * 100)}' for val in np.linspace(0.3, 0.8, 6)])
        # else:
        #     ax1.set_ylim(0.4, 0.8)
        #     ax1.set_yticks(np.linspace(0.4, 0.8, 4))
        #     ax1.set_yticklabels([f'{int(val * 100)}' for val in np.linspace(0.3, 0.6, 4)])

        ax1.set_ylabel(f'% of {SVTYPEMAP[svtype]}')
        ax1.set_xticks(xticks)
        ax1.set_xticklabels([PLATMAP[val] for val in datasets], rotation=90)

        # fig.tight_layout()
        fig1.tight_layout()
        plt.show()

        # fig.savefig(f'{Figure3}/strategy_matched_{svtype}_counts.pdf')
        fig1.savefig(f'{Figure3}/strategy_matched_{svtype}_pcrt.pdf')


def stra_insdel_bpstd(workdir, aligners, datasets):

    shift_labels = ['0,10', '10,50', '50,100', '>100']

    svtypes = ['del', 'ins']

    fig, axes = plt.subplots(1, 2, figsize=(7, 3))
    for col_idx, type_flag in enumerate(svtypes):
        this_ax = axes[col_idx]
        bpstd_list = []
        for dataset in datasets:
            total_matched = 0
            bpstd_dict = {shift: 0 for shift in shift_labels}
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'

            for aligner in aligners:
                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:
                        matched_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra_insdel/{caller}-{asm_caller}.{aligner}-minimap2.{dataset}.{type_flag}.caller_concordant.info.tsv'
                        df_matched_info = pd.read_csv(matched_file, sep='\t', header=[0])
                        total_matched += len(df_matched_info)
                        for idx, row in df_matched_info.iterrows():
                            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                            if start_std <= 10:
                                bpstd_dict['0,10'] += 1
                            elif start_std <= 50:
                                bpstd_dict['10,50'] += 1
                            elif start_std <= 100:
                                bpstd_dict['50,100'] += 1
                            else:
                                bpstd_dict['>100'] += 1

            for shift_label in shift_labels:
                bpstd_list.append((shift_label, bpstd_dict[shift_label] / total_matched, plat))

        df_bpstd = pd.DataFrame(bpstd_list, columns=['label', 'count', 'plat'])
        sns.lineplot(data=df_bpstd, x='label', y='count', hue='plat', style='plat', ax=this_ax, markers=True, lw=2, markersize=9)
        # sns.barplot(data=df_bpstd, x='label', y='count', hue='plat', ax=this_ax)

        this_ax.set_ylabel(f'% of {SVTYPEMAP[type_flag]}')
        this_ax.set_xlabel('')
        this_ax.legend(title='')

    fig.tight_layout()
    plt.show()

def compare_ins_types(workdir, aligners, datasets):
    disc_svtypes_list = []

    for caller in CALLERS:
        for aligner in aligners:
            for dataset in datasets:
                discordant_types = {'pav': {}, 'svimasm': {}}
                concordant_total = {'pav': 0, 'svimasm': 0}
                for asm_caller in ASMCALLERS:
                    concordant_file = f'{workdir}/HiFi/{aligner}_{dataset}/filtered/comstra_insdel/{caller}-{asm_caller}.{aligner}-minimap2.{dataset}.caller_concordants.info.tsv'
                    matched_ins_info = pd.read_csv(concordant_file, sep='\t', header=[0])
                    concordant_total[asm_caller] += len(matched_ins_info)

                    for idx, row in matched_ins_info.iterrows():
                        asm_type_str, algn_type_str = row['TYPE_MATCH'].split('-')[0], row['TYPE_MATCH'].split('-')[1]
                        asm_type_set = set()
                        for ele in asm_type_str.split(','):
                            asm_type_set.add(ele)
                        asm_type = list(asm_type_set)[0]

                        aln_type_set = set()
                        for ele in algn_type_str.split(','):
                            aln_type_set.add(ele)

                        algn_type = ','.join(list(aln_type_set))

                        if asm_type != algn_type:
                            if algn_type not in discordant_types[asm_caller]:
                                discordant_types[asm_caller][algn_type] = 1
                            else:
                                discordant_types[asm_caller][algn_type] += 1

                for asm_caller, count_dict in discordant_types.items():

                    for svtype, count in count_dict.items():
                        pcrt = count * 100 / concordant_total[asm_caller]
                        disc_svtypes_list.append((dataset, caller, aligner, asm_caller, svtype, round(pcrt, 2)))

    df_disc = pd.DataFrame(disc_svtypes_list, columns=['dataset', 'aln_caller', 'aligner', 'asm_caller', 'svtype', 'pcrt'])
    df_disc.to_csv(f'{workdir}/asm_align_disc_ins.tsv', sep='\t', header=True, index=False)
'''
