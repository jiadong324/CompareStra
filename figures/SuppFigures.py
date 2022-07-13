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

def svnum_supps(workdir, datasets, aligners):

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

        ins_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], ins_num * 100 / all_sv_num, plat, aligner, 'Read-based'))
        del_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], del_num * 100 / all_sv_num, plat, aligner, 'Read-based'))
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], all_sv_num, plat, 'Read-based'))

    pav_sv_count = pd.read_csv(f'{workdir}/pav_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in pav_sv_count.iterrows():
        plat = 'HiFi'
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        if 'ont' in dataset:
            plat = 'ONT'
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], svcount, plat, 'Assembly-based'))

        ins_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], ins_num * 100 / svcount, plat, assembler, 'Assembly-based'))
        del_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], del_num * 100 / svcount, plat, assembler, 'Assembly-based'))

    svimasm_sv_count = pd.read_csv(f'{workdir}/svimasm_sv_counts.tsv', header=[0], sep='\t')
    for idx, row in svimasm_sv_count.iterrows():
        plat = 'HiFi'
        caller, dataset, assembler, svcount, ins_num, del_num = row['caller'], row['dataset'], row['assembler'], int(row['total']), int(row['ins_num']), int(row['del_num'])
        if 'ont' in dataset:
            plat = 'ONT'
        sv_count_list.append((TOOLMAP[caller], PLATMAP[dataset], svcount, plat, 'Assembly-based'))
        ins_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], ins_num * 100 / svcount, plat, assembler, 'Assembly-based'))
        del_pcrt.append((TOOLMAP[caller], PLATMAP[dataset], del_num * 100 / svcount, plat, assembler, 'Assembly-based'))

    df_ins_pcrt = pd.DataFrame(ins_pcrt, columns=['caller', 'dataset', 'pcrt', 'plat', 'aligner', 'stra'])
    # df_del_pcrt = pd.DataFrame(del_pcrt, columns=['caller', 'dataset', 'pcrt', 'plat', 'aligner', 'stra'])
    # df_sv_counts = pd.DataFrame(sv_count_list, columns=['caller', 'dataset', 'count', 'plat', 'stra'])

    aligner_colors = [ALIGNERCOLOR[aligner] for aligner in aligners]
    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))

    sns.stripplot(data=df_ins_pcrt[df_ins_pcrt['stra'] == 'Read-based'], x='dataset', y='pcrt', hue='aligner', hue_order=aligners,
                  palette=aligner_colors, size=7, ax=axes[0])

    sns.stripplot(data=df_ins_pcrt[df_ins_pcrt['stra'] == 'Assembly-based'], x='dataset', y='pcrt', hue='aligner', size=7, ax=axes[1])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylim(30, 90)
            ax.set_yticks(np.linspace(30, 90, 4))
            ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(30, 90, 4)], fontsize=12)

            ax.set_ylabel('Percent of INS', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_xlabel('')
        ax.set_xticks(np.arange(len(datasets)))
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=13, rotation=90)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.grid(axis='y', ls='--', color='grey', lw=2)

        ax.legend(title='')

    fig.tight_layout()

    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))

    sns.stripplot(data=df_ins_pcrt[df_ins_pcrt['stra'] == 'Read-based'], x='dataset', y='pcrt', hue='caller',
                  hue_order=[TOOLMAP[caller] for caller in CALLERS], palette=[TOOLCOLORS[caller] for caller in CALLERS],
                  size=7, ax=axes[0])

    sns.stripplot(data=df_ins_pcrt[df_ins_pcrt['stra'] == 'Assembly-based'], x='dataset', y='pcrt', hue='caller',
                  hue_order=[TOOLMAP[caller] for caller in ASMCALLERS], palette=[TOOLCOLORS[caller] for caller in ASMCALLERS],
                  size=7, ax=axes[1])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylim(30, 90)
            ax.set_yticks(np.linspace(30, 90, 4))
            ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(30, 90, 4)], fontsize=12)

            ax.set_ylabel('Percent of INS', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_xlabel('')
        ax.set_xticks(np.arange(len(datasets)))
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=13, rotation=90)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.grid(axis='y', ls='--', color='grey', lw=2)

        ax.legend(title='')

    fig1.tight_layout()

    plt.show()
    fig.savefig(f'{SUPPFIG}/ins_pcrt_of_dataset.pdf')
    fig1.savefig(f'{SUPPFIG}/ins_pcrt_of_dataset_by_callers.pdf')

def histplot_bpstd(workdir, aligner):
    bpstd_list = []

    for asm_caller in ASMCALLERS:

        matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
        for idx, row in matched_info.iterrows():
            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
            if svtype in ['INS', 'DEL'] and int(row['SUPP']) == 6:
                bpstd_list.append((start_std, TOOLMAP[asm_caller], 'Assembly-based'))

    for caller in CALLERS:
        matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
        for idx, row in matched_info.iterrows():
            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])),row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
            if svtype in ['INS', 'DEL'] and int(row['SUPP']) == 6:
                bpstd_list.append((start_std, TOOLMAP[caller], 'Read-based'))

    df_bpstd = pd.DataFrame(bpstd_list, columns=['bpstd', 'caller', 'stra'])
    fig, ax = plt.subplots(1, 1, figsize=(5, 3))

    stra_legends = [Patch(label='Assembly', color=STRACOLORS['Assembly']), Patch(label='Read', color=STRACOLORS['Read'])]

    sns.histplot(data=df_bpstd, x='bpstd', hue='stra', log_scale=[False, True], bins=40,
                 hue_order=['Assembly-based', 'Read-based'], palette=[STRACOLORS['Assembly'], STRACOLORS['Read']], ax=ax)

    # sns.ecdfplot(data=df_bpstd, x='bpstd', hue='stra', hue_order=['Assembly-based', 'Read-based'], palette=[STRACOLORS['Assembly'], STRACOLORS['Read']], ax=ax)
    ax.set_xlim(None, 800)
    ax.set_xticks(np.linspace(0, 800, 5))
    ax.set_xticklabels([int(val) for val in np.linspace(0, 800, 5)], fontsize=12)
    ax.set_xlabel('Breakpoint std.(bp)', fontsize=13)
    ax.legend(handles=stra_legends)

    ax.set_ylabel('# of SVs', fontsize=13)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    plt.show()
    fig.savefig(f'{SUPPFIG}/stra_all_dataset_repro_bpstd_distr.pdf')


