#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/5/5

'''

import matplotlib.pyplot as plt
import pandas as pd
import pysam
import numpy as np
import seaborn as sns
import vcf
from scipy import stats
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

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

def overview_svnum(workdir, datasets):

    ins_pcrt = []
    del_pcrt = []

    sv_count_list = []
    align_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', header=[0], sep='\t')

    for idx, row in align_sv_count.iterrows():
        plat = 'HiFi'
        all_sv_num, ins_num, del_num, aligner, dataset, caller = int(row['all_num']), int(row['ins_num']), int(
            row['del_num']), row['aligner'], row['dataset'], row['caller']
        if 'ont' in dataset:
            plat = 'ONT'

        ins_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], ins_num * 100 / all_sv_num, plat, 'Read-based'))
        del_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], del_num * 100 / all_sv_num, plat, 'Read-based'))
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], all_sv_num, plat, 'Read-based'))

    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        plat = 'HiFi'
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        if 'ont' in dataset:
            plat = 'ONT'
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], svcount, plat, 'Assembly-based'))

        ins_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], ins_num * 100 / svcount, plat, 'Assembly-based'))
        del_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], del_num * 100 / svcount, plat, 'Assembly-based'))

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        plat = 'HiFi'
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        if 'ont' in dataset:
            plat = 'ONT'
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], svcount, plat, 'Assembly-based'))
        ins_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], ins_num * 100 / svcount, plat, 'Assembly-based'))
        del_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], del_num * 100 / svcount, plat, 'Assembly-based'))

    df_ins_pcrt = pd.DataFrame(ins_pcrt, columns=['caller', 'dataset', 'pcrt', 'plat', 'stra'])
    df_del_pcrt = pd.DataFrame(del_pcrt, columns=['caller', 'dataset', 'pcrt', 'plat' ,'stra'])
    df_sv_counts = pd.DataFrame(sv_count_list, columns=['caller', 'dataset', 'count', 'plat', 'stra'])


    fig, axes = plt.subplots(1, 3, sharex='col', sharey='row', figsize=(10, 3))

    xticks = np.arange(len(datasets))
    hue_order = ['Assembly-based', 'Read-based']
    hue_color = [STRACOLORS[ele.split('-')[0]] for ele in hue_order]
    sns.boxplot(data=df_sv_counts, y='dataset', x='count', hue='stra', hue_order=hue_order, palette=hue_color, ax=axes[0])
    axes[0].set_yticks(xticks)
    axes[0].set_yticklabels([PLATMAP[plat] for plat in datasets], fontsize=13)

    axes[0].set_xlim(15000, 45000)
    axes[0].set_xticks(np.linspace(15000, 45000, 4))
    axes[0].set_xticklabels([int(val) for val in np.linspace(15, 45, 4)], fontsize=12)

    axes[0].set_xlabel('Number of SVs (x$10^3$)', fontsize=13)
    axes[0].legend(title='', loc='upper left')

    sns.boxplot(data=df_ins_pcrt, y='dataset', x='pcrt', hue='stra', hue_order=hue_order, palette=hue_color, ax=axes[1])
    sns.boxplot(data=df_del_pcrt, y='dataset', x='pcrt', hue='stra', hue_order=hue_order, palette=hue_color, ax=axes[2])

    axes[1].set_xlim(20, 80)
    axes[1].set_xticks(np.linspace(20, 80, 4))
    axes[1].set_xticklabels([f'{int(val)}%' for val in np.linspace(20, 80, 4)], fontsize=12)
    axes[1].legend('')
    axes[1].set_xlabel('Percent of INS', fontsize=13)

    axes[2].set_xlim(20, 60)
    axes[2].set_xticks(np.linspace(20, 60, 5))
    axes[2].set_xticklabels([f'{int(val)}%' for val in np.linspace(20, 60, 5)], fontsize=12)
    axes[2].legend('')
    axes[2].set_xlabel('Percent of DEL', fontsize=13)

    for ax in axes:
        ax.set_ylabel('')
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')

    fig.tight_layout()

    plt.show()
    fig.savefig(f'{Figure1}/overview_svnums_of_dataset.pdf')

def svsize_distribution(workdir, datasets, aligners):

    assm_svs_small = []
    assm_svs_large = []

    read_svs_small = []
    read_svs_large = []

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for aligner in aligners:
            for read_caller in CALLERS:
                read_calls = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/HG002.{read_caller}.exbnd.bed'
                with open(read_calls, 'r') as f:
                    for line in f:
                        entries = line.strip().split('\t')
                        svsize = abs(int(entries[4]))
                        if svsize <= 1000:
                            read_svs_small.append((svsize, entries[3], read_caller, plat, 'Read-based'))
                        elif svsize <= 10000:
                            read_svs_large.append((svsize, entries[3], read_caller, plat, 'Read-based'))

            if aligner == 'minimap2':
                for asm_caller in ASMCALLERS:
                    asm_calls = f'{workdir}/{plat}/minimap2_{dataset}/filtered/HG002.{asm_caller}.bed'
                    with open(asm_calls, 'r') as f:
                        for line in f:
                            entries = line.strip().split('\t')
                            svsize = abs(int(entries[4]))
                            if svsize <= 1000:
                                assm_svs_small.append((svsize, entries[3], asm_caller, plat, 'Assembly-based'))
                            elif svsize <= 10000:
                                assm_svs_large.append((svsize, entries[3], asm_caller, plat, 'Assembly-based'))

    hue_order = ['HiFi', 'ONT']
    hue_color = [PLATCOLORS[plat] for plat in hue_order]
    plat_legend = [Patch(facecolor=PLATCOLORS[plat], label=plat) for plat in hue_order]

    df_readsvs_small = pd.DataFrame(read_svs_small, columns=['size', 'svtype', 'caller', 'plat', 'stra'])
    df_readsvs_large = pd.DataFrame(read_svs_large, columns=['size', 'svtype', 'caller', 'plat', 'stra'])

    df_assmsvs_small = pd.DataFrame(assm_svs_small, columns=['size', 'svtype', 'caller', 'plat', 'stra'])
    df_assmsvs_large = pd.DataFrame(assm_svs_large, columns=['size', 'svtype', 'caller', 'plat', 'stra'])


    fig, axes = plt.subplots(2, 1, figsize=(5, 3))

    sns.histplot(data=df_assmsvs_small, x='size', hue='plat', bins=100, hue_order=hue_order, palette=hue_color, ax=axes[0])
    # ax.set_title('Assembly-based', fontsize=13)
    axes[0].set_ylim(0, 20000)
    axes[0].set_yticks(np.linspace(0, 20000, 3))
    axes[0].set_yticklabels([f'{int(val)}' for val in np.linspace(0, 20000, 3)])

    sns.histplot(data=df_assmsvs_large, x='size', hue='plat',bins=100,  hue_order=hue_order, palette=hue_color, ax=axes[1])
    # ax.set_title('Assembly-based', fontsize=13)
    axes[1].set_ylim(0, 2000)
    axes[1].set_yticks(np.linspace(0, 2000, 3))
    axes[1].set_yticklabels([f'{int(val)}' for val in np.linspace(0, 2000, 3)])

    for ax in axes:
        ax.set_xlabel('')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylabel('Number of SVs', fontsize=13)
        ax.legend(handles=plat_legend)

    fig.tight_layout()

    fig1, axes = plt.subplots(2, 1, figsize=(5, 3))
    sns.histplot(data=df_readsvs_small, x='size', hue='plat', bins=100,  hue_order=hue_order, palette=hue_color,ax=axes[0])
    # ax.set_title('Assembly-based', fontsize=13)
    axes[0].set_ylim(0, 60000)
    axes[0].set_yticks(np.linspace(0, 60000, 4))
    axes[0].set_yticklabels([f'{int(val)}' for val in np.linspace(0, 60000, 4)])

    sns.histplot(data=df_readsvs_large, x='size', hue='plat', bins=100, hue_order=hue_order, palette=hue_color,ax=axes[1])
    # ax.set_title('Assembly-based', fontsize=13)
    axes[1].set_ylim(0, 6000)
    axes[1].set_yticks(np.linspace(0, 6000, 4))
    axes[1].set_yticklabels([f'{int(val)}' for val in np.linspace(0, 6000, 4)])

    for ax in axes:
        ax.set_xlabel('')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylabel('Number of SVs', fontsize=13)
        ax.legend(handles=plat_legend)


    fig1.tight_layout()
    fig.savefig(f'{Figure1}/assembly_based_svsize.pdf')
    fig1.savefig(f'{Figure1}/read_based_svsize.pdf')
    plt.show()

def svnum_by_regions(workdir, datasets, aligners):
    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    df_sv_counts = pd.read_csv(f'{workdir}/caller_sv_counts_region.tsv', sep='\t', header=0)

    svnums_list = []
    ins_pcrt_list = []
    del_pcrt_list = []

    read_svnums = {caller: {dataset: {aligner: [0 for i in range(len(region_labels))] for aligner in aligners} for dataset in datasets} for caller in CALLERS}
    assm_svnums = {caller: {dataset: {aligner: [0 for i in range(len(region_labels))] for aligner in ['minimap2']} for dataset in datasets} for caller in ASMCALLERS}

    read_svtypes = {
        caller: {dataset: {aligner: [{'INS': 0, 'DEL': 0} for i in range(len(region_labels))] for aligner in aligners} for dataset in
                 datasets} for caller in CALLERS}
    assm_svtypes = {
        caller: {dataset: {aligner: [{'INS': 0, 'DEL': 0} for i in range(len(region_labels))] for aligner in ['minimap2']} for dataset in
                 datasets} for caller in ASMCALLERS}

    for idx, row in df_sv_counts.iterrows():
        caller, dataset, aligner, region, count, svtype = row['caller'], row['dataset'], row['aligner'], row['region'], int(row['count']), row['svtype']
        read_svnums[caller][dataset][aligner][region_labels.index(region)] += count
        if svtype in ['INS', 'DEL']:
            read_svtypes[caller][dataset][aligner][region_labels.index(region)][svtype] += count

    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts_region.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        dataset, region, svcount, aligner, svtype = row['dataset'], row['region'], int(row['count']), row['aligner'], row['svtype']
        assm_svnums['pav'][dataset][aligner][region_labels.index(region)] += svcount
        if svtype in ['INS', 'DEL']:
            assm_svtypes['pav'][dataset][aligner][region_labels.index(region)][svtype] += svcount

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts_region.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        dataset, region, svcount, aligner, svtype = row['dataset'], row['region'], int(row['count']), row['aligner'], row['svtype']
        assm_svnums['svimasm'][dataset][aligner][region_labels.index(region)] += svcount
        if svtype in ['INS', 'DEL']:
            assm_svtypes['svimasm'][dataset][aligner][region_labels.index(region)][svtype] += svcount

    for caller, dataset_aligner in read_svnums.items():
        for dataset, aligner_count in dataset_aligner.items():
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            for aligner, count_list in aligner_count.items():
                sums = np.sum(count_list)
                for idx, count in enumerate(count_list):
                    svnums_list.append((caller, region_labels[idx], count / sums * 100, plat, 'Read-based'))

                    ins_pcrt = read_svtypes[caller][dataset][aligner][idx]['INS'] / sums * 100
                    del_pcrt = read_svtypes[caller][dataset][aligner][idx]['DEL'] / sums * 100

                    ins_pcrt_list.append((caller, region_labels[idx], ins_pcrt, plat, 'Read-based'))
                    del_pcrt_list.append((caller, region_labels[idx], del_pcrt, plat, 'Read-based'))



    for caller, dataset_aligner in assm_svnums.items():
        for dataset, aligner_count in dataset_aligner.items():
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'
            for aligner, count_list in aligner_count.items():
                sums = np.sum(count_list)
                for idx, count in enumerate(count_list):
                    svnums_list.append((caller, region_labels[idx], count / sums * 100, plat, 'Assembly-based'))
                    ins_pcrt = assm_svtypes[caller][dataset][aligner][idx]['INS'] / sums * 100
                    del_pcrt = assm_svtypes[caller][dataset][aligner][idx]['DEL'] / sums * 100

                    ins_pcrt_list.append((caller, region_labels[idx], ins_pcrt, plat, 'Assembly-based'))
                    del_pcrt_list.append((caller, region_labels[idx], del_pcrt, plat, 'Assembly-based'))

    df_svnums = pd.DataFrame(svnums_list, columns=['caller', 'region', 'count', 'plat', 'stra'])
    # df_inspcrt = pd.DataFrame(ins_pcrt_list, columns=['caller', 'region', 'pcrt', 'plat', 'stra'])
    # df_delpcrt = pd.DataFrame(del_pcrt_list, columns=['caller', 'region', 'pcrt', 'plat', 'stra'])

    hue_order = ['Assembly-based', 'Read-based']
    hue_color = [STRACOLORS[ele.split('-')[0]] for ele in hue_order]
    legends = [Patch(facecolor=STRACOLORS[stra.split('-')[0]], label=stra) for stra in hue_order]

    fig, axes = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(5, 3))
    sns.barplot(data=df_svnums, y='region', x='count', hue='stra', hue_order=hue_order, palette=hue_color, ax=axes)

    axes.set_xlabel('Percent of SVs', fontsize=13)
    axes.set_xlim(0, 100)
    axes.set_xticks(np.linspace(0, 100, 5))
    axes.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
    axes.set_ylabel('')
    axes.legend(handles=legends)
    axes.set_yticklabels(axes.get_yticklabels(), fontsize=13)

    # sns.boxplot(data=df_inspcrt, y='region', x='pcrt', hue='stra', hue_order=hue_order, palette=hue_color, ax=axes[1])
    # axes[1].set_xlabel('Percent of INS', fontsize=13)
    # sns.boxplot(data=df_delpcrt, y='region', x='pcrt', hue='stra', hue_order=hue_order, palette=hue_color, ax=axes[2])
    # axes[2].set_xlabel('Percent of DEL', fontsize=13)
    #
    # for idx in [1, 2]:
    #     axes[idx].set_xlim(0, 40)
    #     axes[idx].set_xticks(np.linspace(0, 40, 5))
    #     axes[idx].set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 40, 5)], fontsize=12)
    #
    # for ax in axes:
    #     ax.set_ylabel('')
    #     ax.legend(title='')

    fig.tight_layout()
    plt.show()

    fig.savefig(f'{Figure1}/overview_svs_region_hotspots.pdf')


def sv_regions_pieplot(workdir):
    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    assembly_counts = {region: [] for region in region_labels}
    read_counts = {region: [] for region in region_labels}

    df_sv_counts = pd.read_csv(f'{workdir}/caller_sv_counts_region.tsv', sep='\t', header=0)

    for idx, row in df_sv_counts.iterrows():
        caller, dataset, aligner, region, count, svtype = row['caller'], row['dataset'], row['aligner'], row[
            'region'], int(row['count']), row['svtype']
        read_counts[region].append(count)


    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts_region.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        dataset, region, svcount, assembler, svtype = row['dataset'], row['region'], int(row['count']), row['assembler'], \
                                                    row['svtype']
        assembly_counts[region].append(svcount)

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts_region.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        dataset, region, svcount, assembler, svtype = row['dataset'], row['region'], int(row['count']), row['assembler'], \
                                                    row['svtype']
        assembly_counts[region].append(svcount)


    avg_read_counts = []
    avg_assm_counts = []

    for region in region_labels:
        avg_read_counts.append(np.mean(read_counts[region]))
        avg_assm_counts.append(np.mean(assembly_counts[region]))

    avg_read_outer = [avg_read_counts[0] + avg_read_counts[1] + avg_read_counts[2], avg_read_counts[3]]
    avg_assm_outer = [avg_assm_counts[0] + avg_assm_counts[1] + avg_assm_counts[2], avg_assm_counts[3]]

    size = 0.3
    inner_colors = [REGIONCOLORS[region] for region in region_labels]
    outer_colors = ['#9ac0cd', REGIONCOLORS['Unique']]

    fig, axes = plt.subplots(1, 2, figsize=(6, 4))

    axes[0].pie(avg_read_outer, autopct='%1.1f%%', radius=1, startangle=90, colors=outer_colors, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))
    axes[0].pie(avg_read_counts, autopct='%1.1f%%', radius=1 - size, colors=inner_colors, startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))

    axes[1].pie(avg_assm_outer, autopct='%1.1f%%', radius=1, startangle=90, colors=outer_colors, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))
    axes[1].pie(avg_assm_counts, autopct='%1.1f%%', radius=1 - size, colors=inner_colors, startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))

    plt.tight_layout()

    fig.savefig(f'{Figure1}/sv_regions_avg_pieplot.pdf')
    plt.show()

def read_based_svtypes(workdir, datasets, aligners):

    svtype_colors = {'INS': '#007bba', 'DEL': '#ff7800', 'DUP': '#00a330'}

    df_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', header=[0], sep='\t')

    caller_sv_counts = {'INS': [0 for i in range(len(datasets))],
                               'DEL': [0 for i in range(len(datasets))],
                               'INV': [0 for i in range(len(datasets))],
                               'DUP': [0 for i in range(len(datasets))],
                               'Others': [0 for i in range(len(datasets))]}

    for idx, row in df_sv_count.iterrows():
        all_sv_num, ins_num, del_num, inv_num, dup_num, aligner, dataset, caller = int(row['all_num']), \
                                                                                   int(row['ins_num']), int(
            row['del_num']), int(row['inv_num']), int(row['dup_num']), row['aligner'], \
                                                                                   row['dataset'], row['caller']

        dataset_idx = datasets.index(dataset)

        caller_sv_counts['INS'][dataset_idx] += ins_num
        caller_sv_counts['DEL'][dataset_idx] += del_num
        caller_sv_counts['INV'][dataset_idx] += inv_num
        caller_sv_counts['DUP'][dataset_idx] += dup_num
        caller_sv_counts['Others'][dataset_idx] += all_sv_num - ins_num - del_num - inv_num - dup_num

    svtype_legends = [Patch(facecolor=SVTYPECOLORS[svtype], edgecolor='black', label=svtype) for svtype in
                      ['INS', 'DEL', 'INV', 'DUP', 'Others']]

    barwidth = 0.5
    xticks = np.arange(len(datasets))
    # r2 = [rr + barwidth + 0.06 for rr in r1]
    # rs = [r1, r2]

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))


    ax.bar(xticks, caller_sv_counts['INS'], color=SVTYPECOLORS['INS'], edgecolor='black', width=barwidth)
    ax.bar(xticks, caller_sv_counts['DEL'], color=SVTYPECOLORS['DEL'], bottom=caller_sv_counts['INS'], edgecolor='black',
           width=barwidth)
    ax.bar(xticks, caller_sv_counts['DUP'], color=SVTYPECOLORS['DUP'],
           bottom=[i + j for i, j in zip(caller_sv_counts['INS'], caller_sv_counts['DEL'])], edgecolor='black',
           width=barwidth)
    ax.bar(xticks, caller_sv_counts['INV'], color=SVTYPECOLORS['INV'],
           bottom=[i + j + k for i, j, k in zip(caller_sv_counts['INS'], caller_sv_counts['DEL'], caller_sv_counts['DUP'])],
           edgecolor='black', width=barwidth)
    ax.bar(xticks, caller_sv_counts['Others'], color=SVTYPECOLORS['Others'], bottom=[i + j + k + m for i, j, k, m in
                                                                               zip(caller_sv_counts['INS'],
                                                                                   caller_sv_counts['DEL'],
                                                                                   caller_sv_counts['INV'],
                                                                                   caller_sv_counts['DUP'])],
           edgecolor='black', width=barwidth)

    ax.legend(handles=svtype_legends, ncol=2)
    ax.set_xticks(xticks)
    ax.set_xticklabels([PLATMAP[ele] for ele in datasets], rotation=90, fontsize=12)

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.grid(axis='y', ls='--', color='grey', lw=2)

    # ax.set_ylim(0, 400000)
    # ax.set_yticks(np.linspace(0, 400000, 5))
    # ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 40, 5)], fontsize=12)
    # ax.set_ylabel('Number of SVs (x$10^4$)', fontsize=13)


    plt.tight_layout()
    plt.show()
    fig.savefig(f'{Figure1}/read_based_svtypes.pdf')


def size_distribution_at_region(workdir, datasets):

    asm_methods = ['svimasm', 'pav']
    aligner = 'minimap2'

    hifi_read_calls = []
    hifi_asm_calls = []

    region = 'Simple Repeats'

    for dataset in datasets:
        if 'ont' in dataset:
            for read_caller in CALLERS:
                read_calls = f'{workdir}/ONT/{aligner}_{dataset}/filtered/HG002.{read_caller}.exbnd.bed'
                with open(read_calls, 'r') as f:
                    for line in f:
                        entries = line.strip().split('\t')
                        if entries[5] == region:
                            hifi_read_calls.append((abs(int(entries[4])), entries[3], read_caller))

            for asm_caller in asm_methods:
                asm_calls = f'{workdir}/ONT/minimap2_{dataset}/filtered/HG002.{asm_caller}.bed'
                with open(asm_calls, 'r') as f:
                    for line in f:
                        entries = line.strip().split('\t')
                        if entries[5] == region:
                            hifi_asm_calls.append((abs(int(entries[4])), entries[3], asm_caller))

    df_hifi_reads_calls = pd.DataFrame(hifi_read_calls, columns=['svlen', 'svtype', 'caller'])
    df_hifi_asm_calls = pd.DataFrame(hifi_asm_calls, columns=['svlen', 'svtype', 'caller'])

    fig, hifi_ax = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(10, 4))

    sns.histplot(data=df_hifi_asm_calls, x='svlen', hue='caller', log_scale=True, bins=50, ax=hifi_ax[0])
    sns.histplot(data=df_hifi_reads_calls, x='svlen', hue='caller', log_scale=True, bins=50, ax=hifi_ax[1])

    for ax in hifi_ax:
        ax.set_xlabel('')

    plt.tight_layout()
    plt.show()

