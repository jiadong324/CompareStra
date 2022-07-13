#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/1/31
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


def lineplot_caller_sv_count(workdir, datasets, aligners, figdir):

    sns.set_theme(style="ticks", rc=custom_params, font="Arial", font_scale=1.0)
    colnames = ['caller', 'dataset', 'aligner', 'all_sv_num', 'sv_exbnd_num', 'ins_num', 'del_num', 'inv_num', 'bnd_num']

    df_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', names=colnames, sep='\t')

    sv_count_dict = {aligner: {caller: [0 for i in range(len(datasets))] for caller in CALLERS} for aligner in aligners}


    caller_sv_num = {caller: [] for caller in CALLERS}

    for idx, row in df_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['all_sv_num'])

        if dataset in datasets:
            dataset_idx = datasets.index(dataset)
            sv_count_dict[aligner][caller][dataset_idx] = svcount
            caller_sv_num[caller].append(svcount)

    xticks = np.arange(len(datasets))
    fig, axes = plt.subplots(1, len(aligners), sharex='col', sharey='row', figsize=(12, 4))

    for aligner_idx, aligner in enumerate(aligners):
        ax = axes[aligner_idx]
        ax.set_title(aligner, fontsize=13)
        for caller, counts in sv_count_dict[aligner].items():
            ax.plot(xticks, counts, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], markersize=9, lw=2, label=TOOLMAP[caller])


        if aligner_idx == 0:
            ax.set_ylabel('# of SVs (x1000)', fontsize=13)
        if aligner_idx == 3:
            ax.legend(bbox_to_anchor=(1.05, 1), borderaxespad=1)

        ax.tick_params(width=2)

        for axis in ['bottom', 'left']:
            ax.spines[axis].set_linewidth(2)

        ax.set_ylim(12000, 27000)
        ax.set_yticks(np.linspace(12000, 27000, 4))
        ax.set_yticklabels([int(val) for val in np.linspace(12, 27, 4)], fontsize=12)

        ax.set_xlabel('')
        ax.set_xticks(xticks)
        ax.set_xticklabels([PLATMAP[ele] for ele in datasets], fontsize=12, rotation=90)

    plt.tight_layout()
    plt.show()

    fig.savefig(f'{figdir}/lineplot_sv_num.pdf')


    for caller, count_list in caller_sv_num.items():
        std = np.std(count_list)
        print(f'{caller}: {std}')


def sv_num_pcrt_by_regions(workdir, datasets, aligners):
    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    df_sv_counts = pd.read_csv(f'{workdir}/caller_sv_counts_region.tsv', sep='\t', header=0)

    sv_region_dict = {caller:[[0 for i in range(len(datasets))] for j in range(len(region_labels))] for caller in CALLERS}

    for idx, row in df_sv_counts.iterrows():
        caller, dataset, aligner, region, count = row['caller'], row['dataset'], row['aligner'], row['region'], int(row['count'])

        dataset_idx = datasets.index(dataset)
        region_idx = region_labels.index(region)
        sv_region_dict[caller][region_idx][dataset_idx] += count




    fig, axes = plt.subplots(1, 4, sharex='col', sharey='row',figsize=(8, 4))
    xticks = np.arange(len(datasets))
    bar_width = 0.5

    for col_idx, caller in enumerate(CALLERS):
        ax = axes[col_idx]
        this_caller_data = sv_region_dict[caller]
        for j, region_label in enumerate(region_labels):
            this_region_pcrt = [val / 4 for val in this_caller_data[j]]
            if j == 0:
                ax.bar(xticks, this_region_pcrt, label=region_label, width=bar_width, edgecolor='w')
            else:
                bottoms = []
                for m in range(0, j):
                    bottoms.append([val / 4 for val in this_caller_data[m]])
                bottom_sum = [sum(x) for x in zip(*bottoms)]
                ax.bar(xticks, this_region_pcrt, bottom=bottom_sum, label=region_label, width=bar_width, edgecolor='w')
        ax.set_title(TOOLMAP[caller], fontsize=13)
        if col_idx == 0:
            ax.legend()
        ax.set_xticks(xticks)
        ax.set_xticklabels([PLATMAP[plat] for plat in datasets], fontsize=12, rotation=90)

    plt.tight_layout()
    plt.show()


def stackedplot_num_svtypes(workdir, datasets, aligners, fig_path):

    sns.set_theme(style="ticks", rc=custom_params, font="Arial", font_scale=1.0)

    df_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', header=[0], sep='\t')

    caller_sv_counts = {caller: {aligner: {'INS': [0 for i in range(len(datasets))],
                                                     'DEL': [0 for i in range(len(datasets))],
                                                     'INV': [0 for i in range(len(datasets))],
                                                     'DUP': [0 for i in range(len(datasets))],
                                                     'Others': [0 for i in range(len(datasets))]} for aligner in aligners} for caller in CALLERS}

    for idx, row in df_sv_count.iterrows():
        all_sv_num, ins_num, del_num, inv_num, dup_num, aligner, dataset, caller = int(row['all_num']), \
                                                                            int(row['ins_num']), int(row['del_num']),int(row['inv_num']), int(row['dup_num']), row['aligner'], \
                                                                            row['dataset'], row['caller']

        dataset_idx = datasets.index(dataset)


        caller_sv_counts[caller][aligner]['INS'][dataset_idx] += ins_num
        caller_sv_counts[caller][aligner]['DEL'][dataset_idx] += del_num
        caller_sv_counts[caller][aligner]['INV'][dataset_idx] += inv_num
        caller_sv_counts[caller][aligner]['DUP'][dataset_idx] += dup_num
        caller_sv_counts[caller][aligner]['Others'][dataset_idx] += all_sv_num - ins_num - del_num - inv_num - dup_num


    svtype_legends = [Patch(facecolor=SVTYPECOLORS[svtype], edgecolor='black', label=svtype) for svtype in ['INS', 'DEL', 'INV', 'DUP', 'Others']]

    barwidth = 0.5
    xticks = np.arange(len(datasets))
    # r2 = [rr + barwidth + 0.06 for rr in r1]
    # rs = [r1, r2]

    fig, axes = plt.subplots(len(CALLERS), len(aligners), sharex='col', sharey='row', figsize=(12, 8))

    for row_idx, caller in enumerate(CALLERS):
        for aligner, count_dict in caller_sv_counts[caller].items():
            col_idx = aligners.index(aligner)
            ax = axes[row_idx][col_idx]

            ax.bar(xticks, count_dict['INS'], color=SVTYPECOLORS['INS'], edgecolor='black', width=barwidth)
            ax.bar(xticks, count_dict['DEL'], color=SVTYPECOLORS['DEL'], bottom=count_dict['INS'], edgecolor='black', width=barwidth)
            ax.bar(xticks, count_dict['DUP'], color=SVTYPECOLORS['DUP'], bottom=[i + j for i, j in zip(count_dict['INS'], count_dict['DEL'])], edgecolor='black', width=barwidth)
            ax.bar(xticks, count_dict['INV'], color=SVTYPECOLORS['INV'], bottom=[i + j + k for i, j, k in zip(count_dict['INS'], count_dict['DEL'], count_dict['DUP'])], edgecolor='black', width=barwidth)
            ax.bar(xticks, count_dict['Others'], color=SVTYPECOLORS['Others'], bottom=[i + j + k + m for i, j, k, m in zip(count_dict['INS'], count_dict['DEL'], count_dict['INV'], count_dict['DUP'])],
                   edgecolor='black', width=barwidth)


            ax.legend(handles=svtype_legends, ncol=2)
            ax.set_xticks(xticks)
            ax.set_xticklabels([PLATMAP[ele] for ele in datasets], rotation=90, fontsize=12)

            ax.set_ylim(0, 25000)
            ax.set_yticks(np.linspace(0, 25000, 6))
            ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 25, 6)], fontsize=12)
            if col_idx == 0:
                ax.set_ylabel(f'{TOOLMAP[caller]} SVs (x1000)', fontsize=12)

            if row_idx == 0:
                ax.set_title(aligner, fontsize=13)

    plt.tight_layout()
    plt.show()
    fig.savefig(fig_path)


def caller_other_svtypes(workdir, aligners, datasets):

    fig, axes = plt.subplots(1, len(CALLERS), figsize=(14, 5))

    for col_idx, caller in enumerate(CALLERS):
        this_ax = axes[col_idx]
        svtypes_list = []
        for dataset in datasets:
            svtypes_dict = {}
            for aligner_idx, aligner in enumerate(aligners):
                others_vcf = pysam.VariantFile(f'{workdir}/{aligner}_{dataset}/filtered/HG002.{caller}.others.vcf')

                for record in others_vcf:
                    svtype = record.info['SVTYPE']
                    if svtype in svtypes_dict:
                        svtypes_dict[svtype][aligner_idx] += 1
                    else:
                        svtypes_dict[svtype] = [0 for i in range(len(aligners))]
                        svtypes_dict[svtype][aligner_idx] += 1

            for svtype, counts in svtypes_dict.items():
                for i in range(len(aligners)):
                    svtypes_list.append((dataset, svtype, counts[i], aligners[i]))


        df_svtypes = pd.DataFrame(svtypes_list, columns=['dataset','svtype', 'count', 'aligner'])

        sns.boxplot(data=df_svtypes, x='svtype', y='count', hue='aligner', hue_order=aligners, palette=[ALIGNERCOLOR[ele] for ele in aligners] ,ax=this_ax)
        this_ax.legend(title='')
        this_ax.set_yscale('log')
        if col_idx == 0:
            this_ax.set_ylabel('# of SVs')
        else:
            this_ax.set_ylabel('')
        this_ax.set_xlabel('')
        this_ax.set_title(TOOLMAP[caller])
        this_ax.set_xticklabels(this_ax.get_xticklabels(), rotation=90)

        if col_idx == 0:
            this_ax.set_ylim(10, 10000)
        elif col_idx == 1:
            this_ax.set_ylim(1, 10000)
        elif col_idx == 2:
            this_ax.set_ylim(10, 1000)
        elif col_idx == 4:
            this_ax.set_ylim(100, 10000)
        else:
            this_ax.set_ylim(1, 1000)


    plt.tight_layout()
    plt.show()
    fig.savefig(f'{iMACALIGNFIGURE}/other_svtypes_count.pdf')


def plot_truvari_results(workdir, figdir):

    hifi_datasets = ['hifi_10kb', 'hifi_11kb', 'hifi_15kb', 'hifi_15kb_20kb_27X']
    ont_datasets = ['minion_27X', 'promethion_27X']
    mixed_datasets = ['hifi_10kb', 'hifi_11kb', 'hifi_15kb', 'hifi_15kb_20kb_27X', 'minion_27X', 'promethion_27X']

    result_file = f'{workdir}/recall_precision.txt'

    sns.set_theme(style="ticks", rc=custom_params, font="Arial", font_scale=1.0)
    df_results = pd.read_csv(result_file, sep='\t', names=['plat', 'caller', 'aligner', 'recall', 'precision'])

    boxplot_fscore_hifi_data = []
    boxplot_fscore_mixed_data = []

    boxplot_recall_hifi_data = []
    boxplot_recall_mixed_data = []

    boxplot_precision_hifi_data = []
    boxplot_precision_mixed_data = []


    for idx, row in df_results.iterrows():
        platform, caller, aligner = row['plat'], row['caller'], row['aligner']
        fscore = fmeasure(float(row['precision']), float(row['recall']))

        if platform in mixed_datasets:
            boxplot_fscore_mixed_data.append((platform, caller, aligner, fscore))
            boxplot_recall_mixed_data.append((platform, caller, aligner, float(row['recall']) * 100))
            boxplot_precision_mixed_data.append((platform, caller, aligner, float(row['precision']) * 100))

        if platform in hifi_datasets:
            boxplot_fscore_hifi_data.append((platform, caller, aligner, fscore))
            boxplot_recall_hifi_data.append((platform, caller, aligner, float(row['recall']) * 100))
            boxplot_precision_hifi_data.append((platform, caller, aligner, float(row['precision']) * 100))


    df_fscore_hifi_data = pd.DataFrame(boxplot_fscore_hifi_data, columns=['plat', 'caller', 'aligner', 'fscore'])
    df_fscore_mixed_data = pd.DataFrame(boxplot_fscore_mixed_data, columns=['plat', 'caller', 'aligner', 'fscore'])

    df_recall_hifi_data = pd.DataFrame(boxplot_recall_hifi_data, columns=['plat', 'caller', 'aligner', 'recall'])
    df_recall_mixed_data = pd.DataFrame(boxplot_recall_mixed_data, columns=['plat', 'caller', 'aligner', 'recall'])

    df_precision_hifi_data = pd.DataFrame(boxplot_precision_hifi_data, columns=['plat', 'caller', 'aligner', 'precision'])
    df_precision_mixed_data = pd.DataFrame(boxplot_precision_mixed_data, columns=['plat', 'caller', 'aligner', 'precision'])

    fig, axes = plt.subplots(3, 2, sharex='row', figsize=(10, 10))
    hue_order = ['minimap2', 'ngmlr']
    hue_color = [ALIGNERCOLOR[val] for val in hue_order]

    sns.violinplot(data=df_recall_hifi_data, x='caller', y='recall', hue='aligner', ax=axes[0][0],
                   hue_order=hue_order, palette=hue_color)
    sns.violinplot(data=df_recall_mixed_data, x='caller', y='recall', hue='aligner', ax=axes[0][1],
                   hue_order=hue_order, palette=hue_color)

    axes[0][0].set_title('HiFi', fontsize=13)
    axes[0][1].set_title('HiFi and ONT', fontsize=13)

    sns.violinplot(data=df_precision_hifi_data, x='caller', y='precision', hue='aligner', ax=axes[1][0],
                   hue_order=hue_order, palette=hue_color)
    sns.violinplot(data=df_precision_mixed_data, x='caller', y='precision', hue='aligner', ax=axes[1][1],
                   hue_order=hue_order, palette=hue_color)

    sns.violinplot(data=df_fscore_hifi_data, x='caller', y='fscore', hue='aligner', ax=axes[2][0],
                   hue_order=hue_order, palette=hue_color)
    sns.violinplot(data=df_fscore_mixed_data, x='caller', y='fscore', hue='aligner', ax=axes[2][1],
                   hue_order=hue_order, palette=hue_color)

    for i in range(3):
        for j in range(2):
            ax = axes[i][j]
            ax.set_xlabel('')
            if i < 2:
                ax.set_ylim(40, 100)
                ax.set_yticks(np.linspace(40, 100, 6))
                ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 100, 6)], fontsize=12)
            else:
                ax.set_ylim(0.4, 1)
                ax.set_yticks(np.linspace(0.4, 1, 6))
                ax.set_yticklabels(np.linspace(0.4, 1, 6), fontsize=12)


            ax.legend(loc='lower center')
            if i == 0 and j == 0:
                ax.set_ylabel('Recall', fontsize=13)
            elif i == 1 and j == 0:
                ax.set_ylabel('Precision', fontsize=13)
            elif i == 2 and j == 0:
                ax.set_ylabel('F-score', fontsize=13)
            else:
                ax.set_ylabel('')

            ax.set_xticklabels([TOOLMAP[caller] for caller in CALLERS], fontsize=12)

            ax.tick_params(width=2)
            for axis in ['bottom', 'left']:
                ax.spines[axis].set_linewidth(2)

    plt.tight_layout()
    plt.show()

    fig.savefig(f'{figdir}/evaluation.pdf')

def fmeasure(p, r):
    return 2*p*r / (p+r)



'''
def lineplot_sv_count_diffcov(workdir, figdir):

    datasets = ['hifi_15kb_20kb_27X', 'hifi_15kb_20kb', 'minion_27X', 'minion', 'promethion_27X', 'promethion']

    sns.set_theme(style="ticks", rc=custom_params, font="Arial", font_scale=1.0)
    colnames = ['caller', 'dataset', 'aligner', 'all_sv_num', 'sv_exbnd_num', 'ins_num', 'del_num', 'inv_num', 'bnd_num']

    aligners = ['minimap2', 'ngmlr']


    sv_count_dict = {'minimap2': {caller: [0 for i in range(len(datasets))] for caller in CALLERS},
                     'ngmlr': {caller: [0 for i in range(len(datasets))] for caller in CALLERS}}

    xticklabels = []
    for label in datasets:
        if '27X' in label:
            if 'hifi' in label:
                xticklabels.append('Low-HiFi')
            else:
                xticklabels.append(f'Low-{PLATMAP[label]}')
        else:
            if 'hifi' in label:
                xticklabels.append('High-HiFi')
            else:
                xticklabels.append(f'High-{PLATMAP[label]}')

    lowcov_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', names=colnames, sep='\t')
    for idx, row in lowcov_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['all_sv_num'])
        if dataset in datasets:
            dataset_idx = datasets.index(dataset)
            sv_count_dict[aligner][caller][dataset_idx] = svcount

    highcov_sv_count = pd.read_csv(f'{workdir}/caller_sv_count_highcov.tsv', names=colnames, sep='\t')
    for idx, row in highcov_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['all_sv_num'])
        if dataset in datasets:
            dataset_idx = datasets.index(dataset)
            sv_count_dict[aligner][caller][dataset_idx] = svcount

    xticks = np.arange(len(xticklabels))
    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(9, 4))

    for aligner_idx, aligner in enumerate(aligners):
        ax = axes[aligner_idx]
        ax.set_title(aligner, fontsize=13)
        for caller, counts in sv_count_dict[aligner].items():
            ax.plot(xticks, counts, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], markersize=9, lw=2, label=TOOLMAP[caller])

        if aligner_idx == 0:
            ax.set_ylabel('# of SVs (x1000)', fontsize=13)
        else:
            ax.legend(bbox_to_anchor=(1.05, 1), borderaxespad=0)

        ax.set_ylim(12000, 32000)
        ax.set_yticks(np.linspace(12000, 32000, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(12, 32, 5)], fontsize=12)

        ax.set_xlabel('')
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=12, rotation=90)

    plt.tight_layout()
    plt.show()

    fig.savefig(f'{figdir}/lineplot_sv_num_diffcov.pdf')


def plot_sv_counts_with_giab(workdir, lowcov_align_datasets, highcov_align_datasets, asm_dataset, figdir):
    sns.set_theme(style="ticks", rc=custom_params, font="Arial", font_scale=1.0)

    datasets = []
    datasets.extend(lowcov_align_datasets)
    datasets.extend(highcov_align_datasets)
    # datasets.extend(asm_datasets)

    colnames = ['caller', 'dataset', 'aligner', 'all_sv_num', 'sv_exbnd_num', 'ins_num', 'del_num', 'inv_num', 'bnd_num']

    svcount_dict = {caller: [0 for i in range(len(lowcov_align_datasets) + len(highcov_align_datasets))] for caller in CALLERS}

    df_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', names=colnames, sep='\t')
    for idx, row in df_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['all_sv_num'])
        dataset_idx = datasets.index(dataset)
        svcount_dict[caller][dataset_idx] = svcount


    df_sv_count = pd.read_csv(f'{workdir}/caller_sv_count_highcov.tsv', names=colnames, sep='\t')
    for idx, row in df_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['all_sv_num'])
        dataset_idx = datasets.index(dataset)
        svcount_dict[caller][dataset_idx] = svcount

    xticks = list(np.arange(len(datasets)))

    fig, ax = plt.subplots(1,1, figsize=(6, 4))

    for caller, lowcov_counts in svcount_dict.items():

        ax.plot(xticks, lowcov_counts, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], markersize=9, lw=2,
                label=TOOLMAP[caller])


    df_sv_count = pd.read_csv(f'{workdir}/{asm_dataset}/pav/pav_sv_counts.tsv', sep='\t', names=['tag', 'count'])
    for idx, row in df_sv_count.iterrows():
        if row['tag'] == 'Total':
            ax.axhline(int(row['count']), lw=2, ls='--', color=TOOLCOLORS['giab-pav'])

    ax.set_xticks(xticks)
    labels = []
    for ele in lowcov_align_datasets:
        if '27X' in ele:
            labels.append(f'Low-{PLATMAP[ele]}')
        else:
            labels.append(PLATMAP[ele])

    for ele in highcov_align_datasets:
        labels.append(f'High-{PLATMAP[ele]}')


    ax.set_xticklabels(labels, fontsize=12, rotation=90)
    ax.legend(bbox_to_anchor=(1.05, 1), borderaxespad=0)

    ax.set_ylim(12000, 32000)
    ax.set_yticks(np.linspace(12000, 32000, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(12, 32, 5)], fontsize=12)
    ax.set_ylabel('Number of SVs (x1000)', fontsize=13)

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{figdir}/sv_count_with_giab.pdf')


def caller_insdel_pcrt_scatter(workdir, datasets_list, aligners):

    sns.set_theme(style="ticks", font="Arial", font_scale=1.0)
    colnames = ['caller', 'dataset', 'aligner', 'all_sv_num', 'sv_exbnd_num', 'ins_num', 'del_num', 'inv_num', 'bnd_num']

    df_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', names=colnames, sep='\t')

    fig, axes = plt.subplots(3, len(aligners), sharex='col', sharey='row', figsize=(14, 8))
    row_idx = 0

    for datasets in datasets_list:
        caller_dict = {}
        for idx, row in df_sv_count.iterrows():
            all_sv_num, ins_num, del_num, aligner, dataset, caller = int(row['all_sv_num']), int(row['ins_num']), int(
                row['del_num']), row['aligner'], row['dataset'], row['caller']

            if dataset not in datasets:
                continue

            dataset_idx = datasets.index(dataset)
            ins_pcrt = ins_num / all_sv_num
            del_pcrt = del_num / all_sv_num

            if caller in caller_dict:
                if aligner in caller_dict[caller]:
                    caller_dict[caller][aligner][dataset_idx][0] = ins_pcrt
                    caller_dict[caller][aligner][dataset_idx][1] = del_pcrt
                else:
                    caller_dict[caller][aligner] = [[0, 0] for i in range(len(datasets))]
                    caller_dict[caller][aligner][dataset_idx][0] = ins_pcrt
                    caller_dict[caller][aligner][dataset_idx][1] = del_pcrt
            else:
                caller_dict[caller] = {}
                caller_dict[caller][aligner] = [[0, 0] for i in range(len(datasets))]
                caller_dict[caller][aligner][dataset_idx][0] = ins_pcrt
                caller_dict[caller][aligner][dataset_idx][1] = del_pcrt


        platform_legends = [Patch(facecolor=PLATCOLORS[plat], edgecolor='black', label=PLATMAP[plat]) for plat in datasets]
        caller_legends = [Line2D([0], [0], linestyle='None', label=TOOLMAP[caller],
                                 marker=TOOLMARKERS[caller], markerfacecolor='None', markeredgecolor='black', markersize=12) for caller in CALLERS]



        for i, caller in enumerate(CALLERS):
            for aligner, pcrt_list in caller_dict[caller].items():
                col_idx = aligners.index(aligner)
                this_ax = axes[row_idx][col_idx]

                for j, pcrt in enumerate(pcrt_list):
                    this_ax.scatter(pcrt[0], pcrt[1], color=PLATCOLORS[datasets[j]], marker=TOOLMARKERS[caller], s=150, edgecolor='black')

                if col_idx == 0:
                    this_ax.set_ylabel('% of DEL', fontsize=13)
                    if row_idx == 0:
                        legend1 = this_ax.legend(handles=platform_legends, loc='lower left')
                        legend2 = this_ax.legend(handles=caller_legends, loc='upper left')

                        this_ax.add_artist(legend1)
                        this_ax.add_artist(legend2)
                    else:
                        this_ax.legend(handles=platform_legends, loc='upper left')


                if row_idx == 0:
                    this_ax.set_title(aligner, fontsize=13)
                if row_idx == 2:
                    this_ax.set_xlabel('% of INS', fontsize=13)

                this_ax.set_xlim(0.35, 0.65)
                this_ax.set_xticks(np.linspace(0.35, 0.65, 4))
                this_ax.set_xticklabels([f'{int(val * 100)}%' for val in np.linspace(0.35, 0.65, 4)], fontsize=12)

                this_ax.set_ylim(0.3, 0.5)
                this_ax.set_yticks(np.linspace(0.3, 0.5, 5))
                this_ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0.3, 0.5, 5)], fontsize=12)

        row_idx += 1

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{iMACALIGNFIGURE}/insdel_pcrt_platforms.pdf')
def boxplot_num_svtypes(workdir, datasets, fig_path):

    df_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', header=[0], sep='\t')

    hue_order = ['pbsv', 'svim', 'cutesv', 'sniffles']
    hue_color = [TOOLCOLORS[caller] for caller in hue_order]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    xticks = np.arange(len(datasets))

    axes[0][0].set_title('DEL', fontsize=13)
    sns.boxplot(data=df_sv_count, x='dataset', y='del_num', hue='caller', hue_order=hue_order, palette=hue_color, ax=axes[0][0])
    axes[0][0].set_ylim(7000, 9500, 6)
    axes[0][0].set_yticks(np.linspace(7000, 9500, 6))
    axes[0][0].set_yticklabels([int(val) for val in np.linspace(7000, 9500, 6)], fontsize=12)

    axes[0][1].set_title('INS', fontsize=13)
    sns.boxplot(data=df_sv_count, x='dataset', y='ins_num', hue='caller', hue_order=hue_order, palette=hue_color,ax=axes[0][1])
    axes[0][1].set_ylim(6000, 14000)
    axes[0][1].set_yticks(np.linspace(6000, 14000, 6))
    axes[0][1].set_yticklabels([int(val) for val in np.linspace(6000, 14000, 6)], fontsize=12)

    axes[1][0].set_title('DUP', fontsize=13)
    sns.boxplot(data=df_sv_count, x='dataset', y='dup_num', hue='caller', hue_order=hue_order, palette=hue_color,ax=axes[1][0])
    axes[1][0].set_ylim(0, 5000)
    axes[1][0].set_yticks(np.linspace(0, 5000, 6))
    axes[1][0].set_yticklabels([int(val) for val in np.linspace(0, 5000, 6)], fontsize=12)

    axes[1][1].set_title('INV', fontsize=13)
    sns.boxplot(data=df_sv_count, x='dataset', y='inv_num', hue='caller', hue_order=hue_order, palette=hue_color, ax=axes[1][1])
    axes[1][1].set_ylim(0, 250)
    axes[1][1].set_yticks(np.linspace(0, 250, 6))
    axes[1][1].set_yticklabels([int(val) for val in np.linspace(0, 250, 6)], fontsize=12)

    for i in range(2):
        for j in range(2):
            ax = axes[i][j]
            ax.set_xticks(xticks)
            ax.set_xticklabels([PLATMAP[plat] for plat in datasets], fontsize=12)
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.legend(title='', fontsize=12)


    plt.tight_layout()
    plt.show()
    fig.savefig(fig_path)

'''