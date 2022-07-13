#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/3/17

'''

import matplotlib.pyplot as plt
import pandas as pd
import pysam
import numpy as np
import seaborn as sns
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

def print_sv_num_stats(workdir, datasets, aligners):
    assemblers = ['flye', 'hifiasm', 'shasta']

    read_caller_svnums = {caller: {dataset: {aligner: 0 for aligner in aligners} for dataset in datasets} for caller in CALLERS}
    assm_caller_svnums = {caller: {dataset: {ele: 0 for ele in assemblers} for dataset in datasets} for caller in ASMCALLERS}

    align_caller_counts = {caller: {aligner: [] for aligner in aligners} for caller in CALLERS}
    asm_caller_counts = {caller: [] for caller in ASMCALLERS}

    repeats = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    read_region_svnum = {region: {dataset: {aligner: [0 for i in range(len(CALLERS))] for aligner in aligners} for dataset in datasets} for region in repeats}
    assm_region_svnum = {region: {dataset: {ele: [0 for i in range(len(ASMCALLERS))] for ele in assemblers} for dataset in datasets} for region in repeats}

    read_sv_counts = []
    assm_sv_counts = []
    read_insdel_counts = {'ins':[], 'del': []}
    assm_insdel_counts = {'ins':[], 'del': []}


    read_plat_svnums = {'HiFi': [], 'ONT': []}
    assm_plat_svnums = {'HiFi': [], 'ONT': []}

    align_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', header=[0], sep='\t')

    for idx, row in align_sv_count.iterrows():
        caller, dataset, aligner, svcount, ins_num, del_num = row['caller'], row['dataset'], row['aligner'], int(row['all_num']), int(row['ins_num']), int(row['del_num'])
        align_caller_counts[caller][aligner].append(svcount)
        read_sv_counts.append(svcount)
        read_insdel_counts['ins'].append(ins_num)
        read_insdel_counts['del'].append(del_num)
        read_caller_svnums[caller][dataset][aligner] = svcount

        if 'hifi' in dataset:
            read_plat_svnums['HiFi'].append(svcount)
        elif 'ont' in dataset:
            read_plat_svnums['ONT'].append(svcount)

    df_read_region_svnum = pd.read_csv(f'{workdir}/caller_sv_counts_region.tsv', header=[0], sep='\t')
    for idx, row in df_read_region_svnum.iterrows():
        caller, dataset, aligner, region, svtype, count = row['caller'], row['dataset'], row['aligner'], row['region'], row['svtype'], row['count']
        caller_idx = CALLERS.index(caller)
        read_region_svnum[region][dataset][aligner][caller_idx] += count

    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        asm_caller_counts[caller].append(svcount)
        assm_caller_svnums['pav'][dataset][assembler] = svcount
        assm_sv_counts.append(svcount)
        assm_insdel_counts['ins'].append(ins_num)
        assm_insdel_counts['del'].append(del_num)

        if 'hifi' in dataset:
            assm_plat_svnums['HiFi'].append(svcount)
        elif 'ont' in dataset:
            assm_plat_svnums['ONT'].append(svcount)

    df_pav_region_svnum = pd.read_csv(f'{workdir}/pav_sv_counts_region.tsv', header=[0], sep='\t')
    for idx, row in df_pav_region_svnum.iterrows():
        dataset, assembler, region, svcount = row['dataset'], row['assembler'], row['region'], int(row['count'])
        assm_region_svnum[region][dataset][assembler][ASMCALLERS.index('pav')] += svcount

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        asm_caller_counts[caller].append(svcount)
        assm_caller_svnums['svimasm'][dataset][assembler] = svcount
        assm_sv_counts.append(svcount)
        assm_insdel_counts['ins'].append(ins_num)
        assm_insdel_counts['del'].append(del_num)

        if 'hifi' in dataset:
            assm_plat_svnums['HiFi'].append(svcount)
        elif 'ont' in dataset:
            assm_plat_svnums['ONT'].append(svcount)

    df_svimasm_region_svnum = pd.read_csv(f'{workdir}/svimasm_sv_counts_region.tsv', header=[0], sep='\t')
    for idx, row in df_svimasm_region_svnum.iterrows():
        dataset, assembler, region, svcount = row['dataset'], row['assembler'], row['region'], int(row['count'])
        assm_region_svnum[region][dataset][assembler][ASMCALLERS.index('svimasm')] += svcount

    # for asm_caller, count in asm_caller_counts.items():
    #     print(f'{asm_caller} \n\tstd: {np.std(count)} \n\tmean: {np.mean(count)} \n\t {np.std(count)/np.mean(count)}')

    # for read_caller, count_aligner_dict in align_caller_counts.items():
    #     aligner_avg_count = []
    #     for aligner, count in count_aligner_dict.items():
            # print(f'{read_caller}-{aligner}: {np.std(count)}')
            # aligner_avg_count.append(np.mean(count))
        # print(f'{read_caller} \n\tstd: {np.std(aligner_avg_count)} \n\tmean: {np.mean(aligner_avg_count)} \n\t{np.std(aligner_avg_count)/np.mean(aligner_avg_count)}')

    # print(f'Read-based: \n\tmedian: {np.median(read_sv_counts)} \n\tmean: {np.mean(read_sv_counts)}')
    # for svtype, counts in read_insdel_counts.items():
    #     print(f'\t{svtype} median: {np.median(counts)}')
    #
    # print(f'Assembly-based: \n\tmedian: {np.median(assm_sv_counts)} \n\tmean: {np.mean(assm_sv_counts)}')
    # for svtype, counts in assm_insdel_counts.items():
    #     print(f'\t{svtype} median: {np.median(counts)}')


    # read_region_pcrt = {region: [] for region in repeats}
    # for region, dataset_dict in read_region_svnum.items():
    #     for dataset, aligner_counts in dataset_dict.items():
    #         for aligner, counts in aligner_counts.items():
    #             for i, count in enumerate(counts):
    #                 pcrt = count / read_caller_svnums[CALLERS[i]][dataset][aligner] * 100
    #                 read_region_pcrt[region].append(pcrt)
    #
    # for region, vals in read_region_pcrt.items():
    #     print(f'{region}: {np.mean(vals)}')
    #
    # assm_region_pcrt = {region: [] for region in repeats}
    # for region, dataset_dict in assm_region_svnum.items():
    #     for dataset, aligner_counts in dataset_dict.items():
    #         for aligner, counts in aligner_counts.items():
    #             for i, count in enumerate(counts):
    #                 pcrt = count / assm_caller_svnums[ASMCALLERS[i]][dataset][aligner] * 100
    #                 assm_region_pcrt[region].append(pcrt)
    #
    # for region, vals in assm_region_pcrt.items():
    #     print(f'{region}: {np.mean(vals)}')

    print('Read-based calls =====')
    for plat, vals in read_plat_svnums.items():
        print(f'{plat}: {np.median(vals)}')

    print('Assembly-based calls =====')
    for plat, vals in assm_plat_svnums.items():
        print(f'{plat}: {np.median(vals)}')

def sv_num_stats_of_regions(workdir, datasets, aligners):

    read_caller_counts = {caller: {aligner: [] for aligner in aligners} for caller in CALLERS}
    asm_caller_counts = {caller: [] for caller in ASMCALLERS}

    read_counts = {'CMRG': [], 'PASS': []}
    assm_counts = {'CMRG': [], 'PASS': []}

    cmrg_nums = []
    highconf_nums = []
    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        caller, dataset, aligner, cmrg_count, highconf_num = row['caller'], row['dataset'], row['aligner'], int(row['cmrg_num']), int(row['highconf_num'])
        assm_counts['CMRG'].append(cmrg_count)
        assm_counts['PASS'].append(highconf_num)
        highconf_nums.append((highconf_num, 'Assembly', dataset))
        cmrg_nums.append((cmrg_count, 'Assembly', dataset))

        asm_caller_counts[caller].append(highconf_num)

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        caller, dataset, aligner, cmrg_count, highconf_num = row['caller'], row['dataset'], row['aligner'], int(row['cmrg_num']), int(row['highconf_num'])
        assm_counts['CMRG'].append(cmrg_count)
        assm_counts['PASS'].append(highconf_num)

        highconf_nums.append((highconf_num, 'Assembly', dataset))
        cmrg_nums.append((cmrg_count, 'Assembly', dataset))

        asm_caller_counts[caller].append(highconf_num)

    align_sv_count = pd.read_csv(f'{workdir}/caller_sv_in_cmrg_passed.tsv', header=[0], sep='\t')
    for idx, row in align_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['count'])
        read_counts[row['region_class']].append(svcount)

        if row['region_class'] == 'PASS':
            highconf_nums.append((svcount, 'Read', dataset))
            read_caller_counts[caller][aligner].append(svcount)

        elif row['region_class'] == 'CMRG':
            cmrg_nums.append((svcount, 'Read', dataset))

    # for asm_caller, count in asm_caller_counts.items():
    #     print(f'{asm_caller}: {np.mean(count)}')
    #
    # for read_caller, count_aligner_dict in read_caller_counts.items():
    #     aligner_avg_count = []
    #     for aligner, count in count_aligner_dict.items():
    #         # print(f'{read_caller}-{aligner}: {np.std(count)}')
    #         aligner_avg_count.append(np.mean(count))
    #     print(f'{read_caller}: {np.mean(aligner_avg_count)}')

    print('Read-based calls ======')
    for key, vals in read_counts.items():
        print(f'{key}: {np.median(vals)}')

    print('Assembly-based calls ======')
    for key, vals in assm_counts.items():
        print(f'{key}: {np.median(vals)}')


    df_cmrg = pd.DataFrame(cmrg_nums, columns=['count', 'stra', 'dataset'])
    df_highconf = pd.DataFrame(highconf_nums, columns=['count', 'stra', 'dataset'])

    fig, axes = plt.subplots(1, 2, figsize=(7, 3))
    sns.boxplot(data=df_cmrg, x='dataset', y='count', hue='stra', ax=axes[0])
    axes[0].set_title('CMRG')
    sns.boxplot(data=df_highconf, x='dataset', y='count', hue='stra', ax=axes[1])
    axes[1].set_title('High-confident')

    for ax in axes:
        ax.set_ylabel('')
        ax.set_xlabel('')

    fig.tight_layout()
    plt.show()


def print_regions_stats():

    df_highconf = pd.read_csv('/Users/apple/Evaluation/HG002/HG002_SVs_Tier1_v0.6.pass.bed', sep='\t', names=['chr', 'start', 'end', 'svtype', 'sample'])
    df_cmrgs = pd.read_csv('/Users/apple/Evaluation/HG002/CMRGs/HG002_GRCh37_CMRG_SV_v1.00.bed', sep='\t', names=['chr', 'start', 'end'])

    highconf_size = []

    for idx, row in df_highconf.iterrows():
        size = abs(int(row['end']) - int(row['start']))
        highconf_size.append(size)

    cmrg_size = []
    for idx, row in df_cmrgs.iterrows():
        size = abs(int(row['end']) - int(row['start']))
        cmrg_size.append(size)


    fig, axes = plt.subplots(1, 2, figsize=(7, 3))
    sns.histplot(highconf_size, log_scale=True, ax=axes[0])
    sns.histplot(cmrg_size, log_scale=True, ax=axes[1])

    fig.tight_layout()
    plt.show()

def sv_num_by_aligner(workdir, datasets, aligners, callers, fig_path):

    align_callers = ['pbsv', 'svim', 'cutesv', 'sniffles']

    minimap2_sv_count = {caller: [0 for i in range(len(datasets))] for caller in callers}

    sv_count_dict = {aligner: {caller: [0 for i in range(len(datasets))] for caller in align_callers} for aligner in ['ngmlr', 'lra', 'winnowmap']}

    align_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', header=[0], sep='\t')

    for idx, row in align_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['all_num'])
        dataset_idx = datasets.index(dataset)
        if aligner in sv_count_dict:
            sv_count_dict[aligner][caller][dataset_idx] = svcount
        else:
            minimap2_sv_count[caller][dataset_idx] = svcount


    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['total'])
        dataset_idx = datasets.index(dataset)
        minimap2_sv_count[caller][dataset_idx] = svcount

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['total'])
        dataset_idx = datasets.index(dataset)
        minimap2_sv_count[caller][dataset_idx] = svcount


    xticks = np.arange(len(datasets))
    fig, axes = plt.subplots(1, len(aligners), sharex='col', sharey='row', figsize=(12, 4))


    for aligner_idx, aligner in enumerate(aligners):
        ax = axes[aligner_idx]
        ax.set_title(aligner, fontsize=13)
        if aligner_idx == 0:
            ax.set_ylabel('Number of SVs (x1000)', fontsize=13)
            for caller, counts in minimap2_sv_count.items():
                ax.plot(xticks, counts, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], markersize=9, lw=2, label=TOOLMAP[caller])

        else:
            for caller, counts in sv_count_dict[aligner].items():
                ax.plot(xticks, counts, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], markersize=9, lw=2, label=TOOLMAP[caller])

        # if aligner_idx == 3:
        #     ax.legend(bbox_to_anchor=(1.05, 1), borderaxespad=0)

        ax.legend()
        ax.set_ylim(15000, 35000)
        ax.set_yticks(np.linspace(10000, 35000, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(15, 35, 5)])

        ax.set_xticks(xticks)
        ax.set_xticklabels([PLATMAP[plat] for plat in datasets], fontsize=12, rotation=90)

    plt.tight_layout()
    plt.show()
    fig.savefig(fig_path)

def insdel_pcrt_by_aligner(workdir, datasets, aligners, callers, fig_path):

    minimap2_del_counts = [[0 for i in range(len(datasets))] for i in range(len(callers))]
    minimap2_ins_counts = [[0 for i in range(len(datasets))] for i in range(len(callers))]

    del_counts = {aligner: [[0 for i in range(len(datasets))] for i in range(len(CALLERS))] for aligner in ['ngmlr', 'lra', 'winnowmap']}
    ins_counts = {aligner: [[0 for i in range(len(datasets))] for i in range(len(CALLERS))] for aligner in ['ngmlr', 'lra', 'winnowmap']}

    df_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', header=[0], sep='\t')

    for idx, row in df_sv_count.iterrows():
        all_sv_num, ins_num, del_num, aligner, dataset, caller = int(row['all_num']), int(row['ins_num']), int(
            row['del_num']), row['aligner'], row['dataset'], row['caller']

        caller_idx = CALLERS.index(caller)
        dataset_idx = datasets.index(dataset)

        if aligner != 'minimap2':
            del_counts[aligner][caller_idx][dataset_idx] += del_num * 100 / all_sv_num
            ins_counts[aligner][caller_idx][dataset_idx] += ins_num * 100 / all_sv_num
        else:
            minimap2_del_counts[caller_idx][dataset_idx] += del_num * 100 / all_sv_num
            minimap2_ins_counts[caller_idx][dataset_idx] += ins_num * 100 / all_sv_num

    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        caller, dataset, aligner, svcount, ins_num, del_num = row['caller'], row['dataset'], row['aligner'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        caller_idx = callers.index(caller)
        dataset_idx = datasets.index(dataset)
        minimap2_del_counts[caller_idx][dataset_idx] += del_num * 100 / svcount
        minimap2_ins_counts[caller_idx][dataset_idx] += ins_num * 100 / svcount

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        caller, dataset, aligner, svcount, ins_num, del_num = row['caller'], row['dataset'], row['aligner'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        caller_idx = callers.index(caller)
        dataset_idx = datasets.index(dataset)
        minimap2_del_counts[caller_idx][dataset_idx] += del_num * 100 / svcount
        minimap2_ins_counts[caller_idx][dataset_idx] += ins_num * 100 / svcount


    fig, axes = plt.subplots(2, len(aligners), sharex='col', sharey='row', figsize=(12, 6))
    xticks = np.arange(len(datasets))
    for col_idx, aligner in enumerate(aligners):
        del_ax = axes[0][col_idx]
        del_ax.set_title(aligner, fontsize=13)
        ins_ax = axes[1][col_idx]
        if col_idx == 0:
            del_ax.set_ylabel('% of DEL', fontsize=13)
            ins_ax.set_ylabel('% of INS', fontsize=13)

        del_ax.set_ylim(30, 70)
        del_ax.set_yticks(np.linspace(30, 70, 5))
        del_ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(30, 70, 5)])

        ins_ax.set_ylim(30, 70)
        ins_ax.set_yticks(np.linspace(30, 70, 5))
        ins_ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(30, 70, 5)])

        if aligner == 'minimap2':
            for caller_idx, caller in enumerate(callers):
                del_ax.plot(xticks, minimap2_del_counts[caller_idx], color=TOOLCOLORS[caller],
                            marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])
                ins_ax.plot(xticks, minimap2_ins_counts[caller_idx], color=TOOLCOLORS[caller],
                            marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])
        else:
            this_aligner_del_counts = del_counts[aligner]
            this_aligner_ins_counts = ins_counts[aligner]
            for caller_idx, caller in enumerate(CALLERS):
                del_ax.plot(xticks, this_aligner_del_counts[caller_idx], color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])
                ins_ax.plot(xticks, this_aligner_ins_counts[caller_idx], color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])

        ins_ax.set_xticks(xticks)
        ins_ax.set_xticklabels([PLATMAP[plat] for plat in datasets], rotation=90, fontsize=13)
        del_ax.legend()
        ins_ax.legend()

    plt.tight_layout()
    plt.show()
    fig.savefig(fig_path)


def insdel_pcrt_by_regions(workdir, datasets, aligners, callers):

    sv_count_dict = {aligner: {caller: [0 for i in range(len(datasets))] for caller in callers} for aligner in aligners}

    read_sv_count = pd.read_csv(f'{workdir}/caller_sv_count.tsv', header=[0], sep='\t')

    for idx, row in read_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['all_num'])
        dataset_idx = datasets.index(dataset)
        sv_count_dict[aligner][caller][dataset_idx] += svcount

    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['total'])
        dataset_idx = datasets.index(dataset)
        sv_count_dict[aligner][caller][dataset_idx] += svcount

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['total'])
        dataset_idx = datasets.index(dataset)
        sv_count_dict[aligner][caller][dataset_idx] += svcount

    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    for label in region_labels:

        ins_count_by_regions = {aligner: {caller: [0 for i in range(len(datasets))] for caller in callers} for aligner in aligners}
        del_count_by_regions = {aligner: {caller: [0 for i in range(len(datasets))] for caller in callers} for aligner in aligners}

        read_sv_counts = pd.read_csv(f'{workdir}/caller_sv_counts_region.tsv', sep='\t', header=0)
        for idx, row in read_sv_counts.iterrows():
            caller, dataset, aligner, region, count = row['caller'], row['dataset'], row['aligner'], row['region'], int(row['count'])
            dataset_idx = datasets.index(dataset)
            if label == region:
                if row['svtype'] == 'INS':
                    ins_count_by_regions[aligner][caller][dataset_idx] += count / sv_count_dict[aligner][caller][dataset_idx]
                if row['svtype'] == 'DEL':
                    del_count_by_regions[aligner][caller][dataset_idx] += count / sv_count_dict[aligner][caller][dataset_idx]

        pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts_region.tsv', header=[0], sep='\t')
        for idx, row in pav_sv_count.iterrows():
            dataset, region, svcount = row['dataset'], row['region'], int(row['count'])
            dataset_idx = datasets.index(dataset)
            if region == label:
                if row['svtype'] == 'INS':
                    ins_count_by_regions['minimap2']['pav'][dataset_idx] += svcount / sv_count_dict['minimap2']['pav'][dataset_idx]
                if row['svtype'] == 'DEL':
                    del_count_by_regions['minimap2']['pav'][dataset_idx] += svcount / sv_count_dict['minimap2']['pav'][dataset_idx]

        svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts_region.tsv', header=[0], sep='\t')
        for idx, row in svimasm_sv_count.iterrows():
            dataset, region, svcount = row['dataset'], row['region'], int(row['count'])
            dataset_idx = datasets.index(dataset)
            if region == label:
                if row['svtype'] == 'INS':
                    ins_count_by_regions['minimap2']['svimasm'][dataset_idx] += svcount / sv_count_dict['minimap2']['svimasm'][dataset_idx]
                if row['svtype'] == 'DEL':
                    del_count_by_regions['minimap2']['svimasm'][dataset_idx] += svcount / sv_count_dict['minimap2']['svimasm'][dataset_idx]

        fig, axes = plt.subplots(2, len(aligners), sharex='col', sharey='row', figsize=(12, 6))
        xticks = np.arange(len(datasets))
        for col_idx, aligner in enumerate(aligners):
            del_ax = axes[0][col_idx]
            del_ax.set_title(aligner, fontsize=13)
            ins_ax = axes[1][col_idx]
            if col_idx == 0:
                del_ax.set_ylabel('% of DEL', fontsize=13)
                ins_ax.set_ylabel('% of INS', fontsize=13)

            if aligner == 'minimap2':
                for caller_idx, caller in enumerate(callers):
                    del_ax.plot(xticks, del_count_by_regions[aligner][caller], color=TOOLCOLORS[caller],
                                marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])
                    ins_ax.plot(xticks, ins_count_by_regions[aligner][caller], color=TOOLCOLORS[caller],
                                marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])
            else:
                for caller_idx, caller in enumerate(CALLERS):
                    del_ax.plot(xticks, del_count_by_regions[aligner][caller], color=TOOLCOLORS[caller],
                                marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])
                    ins_ax.plot(xticks, ins_count_by_regions[aligner][caller], color=TOOLCOLORS[caller],
                                marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])

            ins_ax.set_xticks(xticks)
            ins_ax.set_xticklabels([PLATMAP[plat] for plat in datasets], rotation=90, fontsize=13)
            del_ax.legend()
            ins_ax.legend()

        plt.tight_layout()
        plt.show()
        fig.savefig(f'{iMACCOMFIGURE}/{label}_insdel_pcrt.pdf')

def plot_highconf_cmrg_sv_num(workdir, datasets, aligners, callers):


    highconf_minimap2_num = {caller: [0 for i in range(len(datasets))] for caller in callers}

    highconf_count_dict = {aligner: {caller: [0 for i in range(len(datasets))] for caller in CALLERS} for aligner in
                     ['ngmlr', 'lra', 'winnowmap']}

    cmrg_minimap2_num = {caller: [0 for i in range(len(datasets))] for caller in callers}
    cmrg_count_dict = {aligner: {caller: [0 for i in range(len(datasets))] for caller in CALLERS} for aligner in
                           ['ngmlr', 'lra', 'winnowmap']}

    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        caller, dataset, aligner, cmrg_count, highconf_num = row['caller'], row['dataset'], row['aligner'], int(row['cmrg_num']), int(row['highconf_num'])
        dataset_idx = datasets.index(dataset)
        highconf_minimap2_num[caller][dataset_idx] = highconf_num
        cmrg_minimap2_num[caller][dataset_idx] = cmrg_count

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        caller, dataset, aligner, cmrg_count, highconf_num = row['caller'], row['dataset'], row['aligner'], int(row['cmrg_num']), int(row['highconf_num'])
        dataset_idx = datasets.index(dataset)
        highconf_minimap2_num[caller][dataset_idx] = highconf_num
        cmrg_minimap2_num[caller][dataset_idx] = cmrg_count

    align_sv_count = pd.read_csv(f'{workdir}/caller_sv_in_cmrg_passed.tsv', header=[0], sep='\t')
    for idx, row in align_sv_count.iterrows():
        caller, dataset, aligner, svcount = row['caller'], row['dataset'], row['aligner'], int(row['count'])
        dataset_idx = datasets.index(dataset)
        if row['region_class'] == 'PASS':
            if aligner in highconf_count_dict:
                highconf_count_dict[aligner][caller][dataset_idx] = svcount
            else:
                highconf_minimap2_num[caller][dataset_idx] = svcount
        else:
            if aligner in cmrg_count_dict:
                cmrg_count_dict[aligner][caller][dataset_idx] = svcount
            else:
                cmrg_minimap2_num[caller][dataset_idx] = svcount

    xticks = np.arange(len(datasets))
    fig, axes = plt.subplots(1, len(aligners), sharex='col', sharey='row', figsize=(12, 4))

    for aligner_idx, aligner in enumerate(aligners):
        ax = axes[aligner_idx]
        ax.set_title(aligner, fontsize=13)
        if aligner_idx == 0:
            ax.set_ylabel('Number of SVs (x1000)', fontsize=13)
            for caller, counts in highconf_minimap2_num.items():
                ax.plot(xticks, counts, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], markersize=9, lw=2,
                        label=TOOLMAP[caller])

        else:
            for caller, counts in highconf_count_dict[aligner].items():
                ax.plot(xticks, counts, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], markersize=9, lw=2,
                        label=TOOLMAP[caller])

        ax.legend()
        ax.set_ylim(6000, 16000)
        ax.set_yticks(np.linspace(6000, 16000, 6))
        ax.set_yticklabels([int(val) for val in np.linspace(6, 16, 6)])

        ax.set_xticks(xticks)
        ax.set_xticklabels([PLATMAP[plat] for plat in datasets], fontsize=12, rotation=90)

    fig1, axes = plt.subplots(1, len(aligners), sharex='col', sharey='row', figsize=(12, 4))

    for aligner_idx, aligner in enumerate(aligners):
        ax = axes[aligner_idx]
        ax.set_title(aligner, fontsize=13)
        if aligner_idx == 0:
            ax.set_ylabel('Number of SVs (x100)', fontsize=13)
            for caller, counts in cmrg_minimap2_num.items():
                ax.plot(xticks, counts, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], markersize=9, lw=2,
                        label=TOOLMAP[caller])

        else:
            for caller, counts in cmrg_count_dict[aligner].items():
                ax.plot(xticks, counts, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], markersize=9, lw=2,
                        label=TOOLMAP[caller])

        ax.legend()
        ax.set_ylim(100, 500)
        ax.set_yticks(np.linspace(100, 500, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(1, 5, 5)])

        ax.set_xticks(xticks)
        ax.set_xticklabels([PLATMAP[plat] for plat in datasets], fontsize=12, rotation=90)

    fig.tight_layout()
    fig1.tight_layout()
    plt.show()
    fig.savefig(f'{iMACCOMFIGURE}/highconf_region_sv_num.pdf')
    fig1.savefig(f'{iMACCOMFIGURE}/cmrg_region_sv_num.pdf')

