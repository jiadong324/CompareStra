#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/5/19

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


def stra_unique_overview(workdir, aligners, datasets, callers):
    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}

    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    stra_unique_region = []
    stra_unique_svtypes = []

    read_unique_svsize = []
    assm_unique_svsize = []

    stra_uniques = []

    assm_unique_region_dict = {'HiFi': {region: [] for region in regions}, 'ONT': {region: [] for region in regions}}
    read_unique_region_dict = {'HiFi': {region: [] for region in regions}, 'ONT': {region: [] for region in regions}}

    assm_unique_nums = []
    read_unique_nums = []

    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for col_idx, aligner in enumerate(aligners):

            for assembler in plat_assemblers[plat]:
                for caller in CALLERS:
                    for asm_caller in ASMCALLERS:
                        unique_info = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/comstra/{caller}-{aligner}.{asm_caller}-minimap2-{assembler}.{dataset}.uniques.info.tsv'
                        df_unique_info = pd.read_csv(unique_info, sep='\t', header=[0])
                        assm_uniques_by_regions = {region: 0 for region in regions}
                        assm_uniques_svtypes = {}
                        uniques_total = [0, 0]

                        read_uniques_by_regions = {region: 0 for region in regions}
                        read_uniques_svtypes = {}

                        for idx, row in df_unique_info.iterrows():
                            svtype, caller, svlen, svregion, pcrt = row['TYPE_MATCH'], row['CALLER'], abs(int(row['SVLEN'])), row['REGION_TYPE'], float(row['PCRT'])

                            if caller in ASMCALLERS:
                                assm_unique_svsize.append((aligner, assembler, plat, svlen, caller))
                                uniques_total[0] += 1
                                if pcrt >= 50:
                                    assm_uniques_by_regions[svregion] += 1

                                if svtype in assm_uniques_svtypes:
                                    assm_uniques_svtypes[svtype] += 1
                                else:
                                    assm_uniques_svtypes[svtype] = 1

                            elif caller in CALLERS:
                                if svlen >= 50:
                                    read_unique_svsize.append((aligner, assembler, plat, svlen, caller))
                                uniques_total[1] += 1
                                if pcrt >= 50:
                                    read_uniques_by_regions[svregion] += 1
                                if svtype in read_uniques_svtypes:
                                    read_uniques_svtypes[svtype] += 1
                                else:
                                    read_uniques_svtypes[svtype] = 1

                        stra_uniques.append((plat, dataset, aligner, assembler, uniques_total[0], asm_caller, 'Assembly-Uniques'))
                        assm_unique_nums.append(uniques_total[0])
                        stra_uniques.append((plat, dataset, aligner, assembler, uniques_total[1], caller, 'Read-Uniques'))
                        read_unique_nums.append(uniques_total[1])

                        for region, count in assm_uniques_by_regions.items():
                            stra_unique_region.append((plat, region, count / uniques_total[0] * 100, aligner, assembler, 'Assembly-Uniques', asm_caller))
                            assm_unique_region_dict[plat][region].append(count)
                            # assm_unique_region_pcrt[dataset][region].append(count / uniques_total[0] * 100)

                        for svtype, count in assm_uniques_svtypes.items():
                            stra_unique_svtypes.append((plat, svtype, count / uniques_total[0] * 100, aligner, assembler,'Assembly-Uniques', asm_caller))


                        for region, count in read_uniques_by_regions.items():
                            stra_unique_region.append((plat, region, count / uniques_total[1] * 100, aligner, assembler, 'Read-Uniques', caller))
                            read_unique_region_dict[plat][region].append(count)
                            # read_unique_region_pcrt[dataset][region].append(count / uniques_total[1] * 100)

                        for svtype, count in read_uniques_svtypes.items():
                            stra_unique_svtypes.append((plat, svtype, count / uniques_total[1] * 100, aligner, assembler, 'Read-Uniques', caller))

    print(f'Avg. assembly unique counts: {np.mean(assm_unique_nums)}')
    print(f'Avg. read unique counts: {np.mean(read_unique_nums)}')


    df_stra_uniques = pd.DataFrame(stra_uniques, columns=['plat', 'dataset', 'aligner', 'assembler', 'count', 'caller', 'stra'])
    fig, ax = plt.subplots(1, 1, figsize=(3, 4))

    plat_order = ['HiFi', 'ONT']
    plat_color = [PLATCOLORS[plat] for plat in plat_order]
    plat_legends = [Patch(facecolor=PLATCOLORS[plat], label=plat) for plat in plat_order]

    assemblers = ['shasta', 'flye', 'hifiasm']
    assm_colors = [ASSMBLERCOLOR[ele] for ele in assemblers]

    sns.stripplot(data=df_stra_uniques, x='stra', y='count', hue='plat', hue_order=plat_order, palette=plat_color, ax=ax)
    box = sns.boxplot(data=df_stra_uniques, y='count', x='stra', ax=ax)
    for i in [0, 1]:
        mybox = box.artists[i]
        mybox.set_facecolor('white')
        mybox.set_edgecolor('black')

    ax.set_xlabel('')
    ax.set_xticks(np.arange(2))
    ax.set_xticklabels(['Assembly', 'Read'], fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.legend()

    ax.set_ylim(2000, 30000)
    ax.set_yticks(np.linspace(2000, 30000, 5))
    ax.set_yticklabels([f'{int(val)}' for val in np.linspace(2, 30, 5)], fontsize=12)
    ax.set_ylabel('# of SVs', fontsize=13)

    fig.tight_layout()
    fig.savefig(f'{Figure5}/stra_unique_boxplot.pdf')

    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))
    aligner_color = [ALIGNERCOLOR[ele] for ele in aligners]

    sns.boxplot(data=df_stra_uniques[df_stra_uniques['stra'] == 'Assembly-Uniques'], x='assembler', y='count',
                hue_order=aligners, palette=aligner_color, hue='aligner', ax=axes[0])
    axes[0].set_title('Assembly-Uniques')

    sns.boxplot(data=df_stra_uniques[df_stra_uniques['stra'] == 'Read-Uniques'], x='assembler', y='count',
                hue_order=aligners, palette=aligner_color, hue='aligner', ax=axes[1])
    axes[1].set_title('Read-Uniques')

    for i, ax in enumerate(axes):
        ax.set_xlabel('')
        ax.set_xticks(np.arange(3))
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=13)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.legend()
        if i == 0:
            ax.set_ylim(2000, 30000)
            ax.set_yticks(np.linspace(2000, 30000, 5))
            ax.set_yticklabels([f'{int(val)}' for val in np.linspace(2, 30, 5)], fontsize=12)
            ax.set_ylabel('# of SVs', fontsize=13)
        else:
            ax.set_ylabel('')

    fig1.tight_layout()
    fig1.savefig(f'{Figure5}/strategy_uniques_aligner_assembler_count.pdf')

    # df_read_uniques_svsize = pd.DataFrame(read_unique_svsize, columns=['aligner', 'assembler', 'plat', 'svlen', 'stra'])
    # df_assm_uniques_svsize = pd.DataFrame(assm_unique_svsize, columns=['aligner', 'assembler', 'plat', 'svlen', 'stra'])

    # plat_order = ['HiFi', 'ONT']
    # plat_color = [PLATCOLORS[ele] for ele in plat_order]
    # plat_legends = [Patch(facecolor=PLATCOLORS[ele], label=ele) for ele in plat_order]

    # fig2, axes = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(5, 4))
    # sns.histplot(data=df_read_uniques_svsize[df_read_uniques_svsize['svlen'] < 1000], x='svlen', hue='plat',
    #              hue_order=plat_order, palette=plat_color, bins=50, ax=axes[0])
    # axes[0].set_ylabel('Read unique size', fontsize=13)
    #
    # sns.histplot(data=df_assm_uniques_svsize[df_assm_uniques_svsize['svlen'] < 1000], x='svlen', hue='plat',
    #              hue_order=plat_order, palette=plat_color, bins=50, ax=axes[1])
    # axes[1].set_ylabel('Assembly unique size', fontsize=13)
    #
    # for ax in axes:
    #     ax.set_xlabel('')
    #     ax.set_ylabel('Number of SVs')
    #     ax.set_ylim(0, 500000)
    #     ax.set_yticks(np.linspace(0, 500000, 3))
    #     ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 500, 3)], fontsize=12)
    #
    #     ax.spines['top'].set_visible(False)
    #     ax.spines['right'].set_visible(False)
    #     ax.set_xlabel('')
    #
    #     ax.legend(handles=plat_legends)
    #
    # fig2.tight_layout()
    #
    # fig2.savefig(f'{Figure5}/stra_unique_svsize.pdf')

    fig3, axes = plt.subplots(1, 2, figsize=(6, 4))
    size= 0.3
    colors = [REGIONCOLORS[region] for region in regions]
    for fig_idx, plat in enumerate(['HiFi', 'ONT']):
        this_ax = axes[fig_idx]

        read_region_avg = []
        assm_region_avg = []
        for region in regions:
            read_region_avg.append(np.mean(read_unique_region_dict[plat][region]))
            assm_region_avg.append(np.mean(assm_unique_region_dict[plat][region]))

        print(read_region_avg)
        print(assm_region_avg)
        this_ax.pie(assm_region_avg, autopct='%1.1f%%', radius=1, colors=colors,
                    startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))

        this_ax.pie(read_region_avg, autopct='%1.1f%%', radius=1 - size, colors=colors,
                    startangle=90, pctdistance=0.8, wedgeprops=dict(width=size, edgecolor='w'))

    fig3.tight_layout()
    fig3.savefig(f'{Figure5}/stra_unique_region_pieplot.pdf')

    plt.show()

def plot_merged_stra_uniques(datasets):

    xticklabels = ['SUPP=1', 'SUPP=2', 'SUPP=3', 'SUPP=4', 'SUPP=5', 'SUPP=6', 'SUPP=7', 'SUPP=8', 'SUPP=9', 'SUPP=10']

    count_legend = [Line2D([0], [0], color=PLATCOLORS['HiFi'], label=f'Count (HiFi)', lw=2, marker='o'),
                    Line2D([0], [0], color=PLATCOLORS['ONT'], label=f'Count (ONT)', lw=2, marker='X')]

    for stra in ['read', 'assembly']:
        count_dict = {'HiFi': [0 for i in range(len(xticklabels))], 'ONT': [0 for i in range(len(xticklabels))]}
        pcrt_dict = {'HiFi': [0 for i in range(len(xticklabels))], 'ONT': [0 for i in range(len(xticklabels))]}

        for dataset in datasets:
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'

            print(f'{plat} =====')
            output_dir = f'{iMACDIR}/{plat}/minimap2_{dataset}/filtered/comstra'

            unique_supp, total = get_survivor_supp(f'{output_dir}/{stra}_based.uniques.jasmine.merged.vcf')

            print(f'{stra} total: {total}')

            missed = 0

            for supp, count in unique_supp.items():
                supp_idx = xticklabels.index(f'SUPP={supp}')
                missed += count
                count_dict[plat][supp_idx] += count
                pcrt_dict[plat][supp_idx] += count / total * 100

            print(f'% of {stra}-based supp=8: {missed}, {missed / total}')

        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
        barwidth = 0.3
        r1 = np.arange(len(xticklabels))
        r2 = [r + barwidth for r in r1]
        xticks = [r + barwidth / 2 for r in r1]

        ax.bar(r1, pcrt_dict['HiFi'], width=barwidth, edgecolor='white', facecolor=PLATCOLORS['HiFi'], label='Percent (HiFi)')
        ax.bar(r2, pcrt_dict['ONT'], width=barwidth, edgecolor='white', facecolor=PLATCOLORS['ONT'], label='Percent (ONT)')

        ax.set_xlabel('# of supporting strategy unique sets', fontsize=13)
        ax.set_xticks(xticks)
        ax.set_xticklabels([int(val) for val in np.linspace(1, 10, 10)], fontsize=13)

        ax.set_ylabel(f'% of {stra} unique SVs', fontsize=13)
        ax.legend(loc='upper right')
        if stra == 'read':
            ax.set_ylim(0, 80)
            ax.set_yticks(np.linspace(0, 80, 5))
            ax.set_yticklabels([int(val) for val in np.linspace(0, 80, 5)], fontsize=12)
        else:
            ax.set_ylim(0, 40)
            ax.set_yticks(np.linspace(0, 40, 5))
            ax.set_yticklabels([int(val) for val in np.linspace(0, 40, 5)], fontsize=12)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # ax1 = ax.twinx()
        # ax1.plot(r1, count_dict['HiFi'], lw=2, marker='o', color=PLATCOLORS['HiFi'])
        # ax1.plot(r2, count_dict['ONT'], lw=2, marker='X',  color=PLATCOLORS['ONT'])
        #
        # ax1.set_ylabel(f'# of {stra} unique SVs (x100)', fontsize=13)
        # if stra == 'assembly':
        #     ax1.set_ylim(0, 7500)
        #     ax1.set_yticks(np.linspace(0, 7500, 4))
        #     ax1.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 75, 4)], fontsize=12)
        # else:
        #     ax1.set_ylim(0, 20000)
        #     ax1.set_yticks(np.linspace(0, 20000, 5))
        #     ax1.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 200, 5)], fontsize=12)
        #
        # ax1.legend(handles=count_legend, loc='upper left')

        fig.tight_layout()
        plt.show()
        fig.savefig(f'{Figure5}/merged_{stra}_uniques_supp.pdf')


def overview_unique_annot(workdir, datasets):
    columns = ['chrom', 'start', 'end', 'svtype', 'svlen', 'region', 'rptype', 'pcrt', 'suppvec', 'supp', 'mapq', 'sigs', 'nsigs']
    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    annot_info_list = []
    region_dict = {'HiFi': {'read': {region: 0 for region in regions}, 'assembly': {region: 0 for region in regions}},
                   'ONT': {'read': {region: 0 for region in regions}, 'assembly': {region: 0 for region in regions}}}
    region_list = []
    total_counts = {'HiFi': {'read': 0, 'assembly': 0},
                    'ONT': {'read': 0, 'assembly': 0}}

    read_svtype_dict = {'INS': [0, 0], 'DEL': [0, 0]}
    assembly_svtype_dict = {'INS': [0, 0], 'DEL': [0, 0]}

    stra_map = {'read': 'Read', 'assembly': 'Assembly'}
    plats = ['HiFi', 'ONT']

    insdel_info_list = []

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        plat_idx = plats.index(plat)
        for stra in ['read', 'assembly']:
            df_uniques = pd.read_csv(f'{workdir}/{plat}/minimap2_{dataset}/filtered/comstra/{stra}_based.uniques.annot.sigs.tsv', sep='\t', names=columns)
            svtypes_dict = {'INS': 0, 'DEL': 0}
            totals = 0

            for idx, row in df_uniques.iterrows():
                this_annot = [row['svtype'], row['svlen'], row['region'], row['pcrt'], row['mapq'], row['sigs'], row['nsigs'], plat, stra_map[stra]]
                if int(row['supp']) == 10:
                    totals += 1
                    annot_info_list.append(this_annot)
                    region_dict[plat][stra][row['region']] += 1
                    total_counts[plat][stra] += 1

                    if row['svtype'] == 'INS':
                        svtypes_dict['INS'] += 1
                        insdel_info_list.append(('INS', stra_map[stra], int(row['svlen'] + 1)))
                        if stra == 'read':
                            read_svtype_dict[row['svtype']][plat_idx] += 1
                        else:
                            assembly_svtype_dict[row['svtype']][plat_idx] += 1

                    if row['svtype'] == 'DEL':
                        svtypes_dict['DEL'] += 1
                        insdel_info_list.append(('DEL', stra_map[stra], int(row['svlen'] + 1)))
                        if stra == 'read':
                            read_svtype_dict[row['svtype']][plat_idx] += 1
                        else:
                            assembly_svtype_dict[row['svtype']][plat_idx] += 1

            print(f'{plat} {stra} {totals}')
            print('\t #INS:', svtypes_dict['INS'])
            print('\t #DEL:', svtypes_dict['DEL'])

    legends = [Line2D([0], [0], color='black', ls='-',label='HiFi', lw=2),
               Line2D([0], [0], color='black', ls='--', label='ONT', lw=2),
               Line2D([0], [0], color=STRACOLORS['Read'], label='Read-based', lw=2),
               Line2D([0], [0], color=STRACOLORS['Assembly'], label='Assembly-based', lw=2)]

    xticks = np.arange(len(regions))

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    for plat, stra_region_count in region_dict.items():
        for stra, region_count in stra_region_count.items():
            pcrt_list = []
            for region, count in region_count.items():
                pcrt_list.append(count / total_counts[plat][stra] * 100)
                region_list.append((plat, f'{stra_map[stra]}-based', region, count / total_counts[plat][stra] * 100))
                print(f'{plat} {stra} {region} {count} ({count / total_counts[plat][stra] * 100})')

            ax.plot(pcrt_list, xticks, ls=PLATLS[plat], lw=2, marker='o', color=STRACOLORS[stra_map[stra]])

    ax.set_xlabel('Percent of unique SVs', fontsize=13)
    ax.set_xlim(0, 100)
    ax.set_xticks(np.linspace(0, 100, 5))
    ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
    ax.set_ylabel('')
    ax.set_yticks(np.arange(len(regions)))
    ax.set_yticklabels(regions, fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.legend(handles=legends)

    fig.tight_layout()

    legends = [Line2D([0], [0], color='black', ls='--', lw=2, label='ONT'),
               Line2D([0], [0], color='black', lw=2, label='HiFi'),
               Line2D([0], [0], color=STRACOLORS['Assembly'], lw=2, label='Assembly-based'),
               Line2D([0], [0], color=STRACOLORS['Read'], lw=2, label='Read-based'),]

    df_annot_info = pd.DataFrame(annot_info_list, columns=['svtype', 'svlen', 'region', 'pcrt', 'mapq', 'sigs', 'nsigs', 'plat', 'stra'])
    fig1, ax = plt.subplots(1, 1, figsize=(5, 4))
    hue_order = ['Read', 'Assembly']
    hue_color = [STRACOLORS[val] for val in hue_order]
    sns.ecdfplot(data=df_annot_info[df_annot_info['plat'] == 'HiFi'], x='mapq', hue='stra', hue_order=hue_order, palette=hue_color, ax=ax, lw=2)
    ax.set_ylabel('Percent of unique SVs', fontsize=13)
    ax.set_ylim(0, 1)
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)], fontsize=12)
    sns.ecdfplot(data=df_annot_info[df_annot_info['plat'] == 'ONT'], ls='--', x='mapq', hue='stra', hue_order=hue_order, palette=hue_color, ax=ax, lw=2)
    ax.axvline(20, 0, 0.4, color='#ce9d35', lw=2)
    ax.set_xlabel('Mapping quality', fontsize=13)
    ax.legend(handles=legends)
    fig1.tight_layout()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    svtype_legends = [Patch(facecolor=SVTYPECOLORS['INS'], edgecolor='white', label='INS'),
                      Patch(facecolor=SVTYPECOLORS['DEL'], edgecolor='white', label='DEL')]

    fig2, ax = plt.subplots(1, 1, figsize=(4, 2))
    xticks = np.arange(len(plats))
    ax.barh(xticks, read_svtype_dict['INS'], height=0.5, color=SVTYPECOLORS['INS'])
    ax.barh(xticks, read_svtype_dict['DEL'], height=0.5, left=read_svtype_dict['INS'], color=SVTYPECOLORS['DEL'])
    ax.set_xlim(0, 400)
    ax.set_xticks(np.linspace(0, 400, 5))
    ax.set_xticklabels([int(val) for val in np.linspace(0, 400, 5)], fontsize=12)
    ax.set_xlabel('Number of SVs', fontsize=13)
    ax.set_yticks(xticks)
    ax.set_yticklabels(plats, fontsize=13)
    ax.legend(handles=svtype_legends)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig2.tight_layout()

    fig3, ax = plt.subplots(1, 1, figsize=(4, 2))
    xticks = np.arange(len(plats))
    ax.barh(xticks, assembly_svtype_dict['INS'], height=0.5, color=SVTYPECOLORS['INS'])
    ax.barh(xticks, assembly_svtype_dict['DEL'], height=0.5, left=assembly_svtype_dict['INS'], color=SVTYPECOLORS['DEL'])
    ax.set_xlim(0, 3000)
    ax.set_xticks(np.linspace(0, 3000, 4))
    ax.set_xticklabels([int(val) for val in np.linspace(0, 3000, 4)], fontsize=12)
    ax.set_xlabel('Number of SVs', fontsize=13)
    ax.set_yticks(xticks)
    ax.set_yticklabels(plats, fontsize=13)
    ax.legend(handles=svtype_legends)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig3.tight_layout()

    df_insdel = pd.DataFrame(insdel_info_list, columns=['svtype', 'stra', 'svlen'])
    svtype_order = ['INS', 'DEL']
    svtype_color = [SVTYPECOLORS['INS'], SVTYPECOLORS['DEL']]

    fig4, axes = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(6, 4))
    sns.histplot(data=df_insdel[df_insdel['stra'] == 'Read'], x='svlen', log_scale=True, bins=50, hue='svtype', hue_order=svtype_order, palette=svtype_color, ax=axes[0])
    axes[0].set_ylabel('Read', fontsize=13)
    axes[0].set_ylim(0, 60)
    axes[0].set_yticks(np.linspace(0, 60, 4))
    axes[0].set_yticklabels([int(val) for val in np.linspace(0, 60, 4)], fontsize=12)

    sns.histplot(data=df_insdel[df_insdel['stra'] == 'Assembly'], x='svlen', log_scale=True, bins=50, hue='svtype', hue_order=svtype_order, palette=svtype_color, ax=axes[1])
    axes[1].set_ylabel('Assembly', fontsize=13)
    axes[1].set_ylim(0, 300)
    axes[1].set_yticks(np.linspace(0, 300, 4))
    axes[1].set_yticklabels([int(val) for val in np.linspace(0, 300, 4)], fontsize=12)

    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(handles=svtype_legends)
        ax.set_xlabel('')

    fig4.tight_layout()
    plt.show()

    fig.savefig(f'{Figure5}/strategy_uniques_regions.pdf')
    fig1.savefig(f'{Figure5}/strategy_uniques_mapq_ecdf.pdf')
    fig2.savefig(f'{Figure5}/read_uniques_insdel_count.pdf')
    fig3.savefig(f'{Figure5}/assembly_uniques_insdel_count.pdf')
    fig4.savefig(f'{Figure5}/strategy_uniques_insdel_size.pdf')

def plot_assm_annot_results(workdir, datasets):

    columns = ['chrom', 'start', 'end', 'svtype', 'svlen', 'region', 'rptype','pcrt', 'suppvec', 'supp', 'mapq', 'sigs', 'nsigs']
    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    annot_info_list = []
    total_counts = {'HiFi': 0, 'ONT': 0}
    region_dict = {'HiFi': {region: 0 for region in regions},
                   'ONT': {region: 0 for region in regions}}

    svtypes_dict = {}

    high_mapqs = {'HiFi': [], 'ONT': []}
    high_mapqs_has_sigs = {'HiFi': [], 'ONT': []}
    high_mapqs_no_sigs = {'HiFi': [], 'ONT': []}

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        sig_annot = f'{workdir}/{plat}/minimap2_{dataset}/filtered/comstra/assembly_based.uniques.annot.sigs.tsv'
        df_sig_annot = pd.read_csv(sig_annot, sep='\t', names=columns)

        for idx, row in df_sig_annot.iterrows():
            this_annot = [row['svtype'], row['svlen'], row['region'], row['pcrt'], row['mapq'], row['sigs'], row['nsigs'], plat]
            if int(row['supp']) == 10:
                # annot_info_list.append(this_annot)
                region_dict[plat][row['region']] += 1
                total_counts[plat] += 1

                if row['svtype'] in svtypes_dict:
                    svtypes_dict[row['svtype']] += 1
                else:
                    svtypes_dict[row['svtype']] = 1

                if int(row['mapq']) >= 20:
                    annot_info_list.append(this_annot)
                    high_mapqs[plat].append(this_annot)
                    if row['sigs'] >= 5:
                        high_mapqs_has_sigs[plat].append(this_annot)
                    else:
                        high_mapqs_no_sigs[plat].append(this_annot)

    for plat, info_list in high_mapqs.items():
        print(f'{plat}: {len(info_list)}')
        print(f'{len(high_mapqs_has_sigs[plat])} ({len(high_mapqs_has_sigs[plat]) / len(info_list) * 100})')
        print(f'{len(info_list)} ({len(info_list) / total_counts[plat] * 100})')

    # region_list = []
    # for plat, counts_dict in region_dict.items():
    #     print(f'{plat} ======')
    #     for region, count in counts_dict.items():
    #         region_list.append((plat, region, count / total_counts[plat] * 100))
    #         print(f'{region} {count} ({count / total_counts[plat] * 100})')
    #
    # print(svtypes_dict)

    legends = [Line2D([0], [0], color='black', ls='-', lw=2, label='HiFi'),
               Line2D([0], [0], color='black', ls='--', lw=2, label='ONT')]

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    xticks = np.arange(len(regions))

    for plat, info_list in high_mapqs_no_sigs.items():
        total_nums = len(info_list)
        # print(f'{plat} =====')
        region_dict = {region: 0 for region in regions}
        for ele in info_list:
            region_dict[ele[2]] += 1

        pcrt_list = []
        for region, counts in region_dict.items():
            # print(f'{region} : {counts}, {counts / total_nums * 100}')
            pcrt_list.append(counts / total_nums * 100)

        ax.plot(pcrt_list, xticks, ls=PLATLS[plat], marker='o', lw=2, color=STRACOLORS['Assembly'])

    ax.set_yticks(xticks)
    ax.set_yticklabels(regions, fontsize=13)

    ax.set_xlim(0, 80)
    ax.set_xticks(np.linspace(0, 80, 5))
    ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 80, 5)], fontsize=12)
    ax.set_xlabel('Percent of assembly unique SVs', fontsize=13)

    ax.legend(handles=legends)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    plt.show()
    fig.savefig(f'{Figure5}/assembly_unique_nosigs_regions.pdf')


def plot_read_annot_results(workdir, datasets):

    columns = ['chrom', 'start', 'end', 'svtype', 'svlen', 'region', 'rptype', 'pcrt', 'suppvec', 'supp', 'mapq', 'sigs', 'nsigs']
    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    annot_info_list = []
    total_counts = {'HiFi': 0, 'ONT': 0}
    region_dict = {'HiFi': {region: 0 for region in regions},
                   'ONT': {region: 0 for region in regions}}

    svtypes_dict = {}

    high_mapqs = {'HiFi': [], 'ONT': []}
    high_mapqs_has_sigs = {'HiFi': [], 'ONT': []}

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        sig_annot = f'{workdir}/{plat}/minimap2_{dataset}/filtered/comstra/read_based.uniques.annot.sigs.tsv'
        df_sig_annot = pd.read_csv(sig_annot, sep='\t', names=columns)

        for idx, row in df_sig_annot.iterrows():
            this_annot = [row['svtype'], row['svlen'], row['region'], row['pcrt'], row['mapq'], row['sigs'], row['nsigs'], plat]
            if int(row['supp']) == 8:
                # annot_info_list.append(this_annot)
                region_dict[plat][row['region']] += 1
                total_counts[plat] += 1

                if row['svtype'] in svtypes_dict:
                    svtypes_dict[row['svtype']] += 1
                else:
                    svtypes_dict[row['svtype']] = 1

                if int(row['mapq']) >= 20:
                    annot_info_list.append(this_annot)
                    high_mapqs[plat].append(this_annot)

    for plat, info_list in high_mapqs.items():
        print(f'{plat}: {len(info_list)} ({len(info_list) / total_counts[plat] * 100})')

    region_list = []
    for plat, counts_dict in region_dict.items():
        print(f'{plat} ======')
        for region, count in counts_dict.items():
            region_list.append((plat, region, count / total_counts[plat] * 100))
            print(f'{region} {count} ({count / total_counts[plat] * 100})')

    print(svtypes_dict)

    # df_region_counts = pd.DataFrame(region_list, columns=['plat', 'region', 'count'])
    # fig1, ax = plt.subplots(1, 1, figsize=(5, 4))
    # sns.barplot(data=df_region_counts, x='region', y='count', hue='plat', ax=ax)
    # ax.set_ylabel('% of SVs')
    # ax.set_ylim(0, 100)
    # ax.set_yticks(np.linspace(0, 100, 5))
    # ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
    #
    # ax.set_xlabel('')
    # ax.set_xticklabels(ax.get_xticklabels(), fontsize=13, rotation=90)
    # ax.legend(title='')
    # fig1.tight_layout()
    #
    # plt.show()