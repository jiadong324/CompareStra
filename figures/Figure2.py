#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/5/6

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



def print_platform_unique_stats(workdir, aligner):

    unique_svs = []
    xticklabels = []
    read_caller_counts = {'hifi': [], 'ont': []}
    asm_caller_counts = {'hifi': [], 'ont': []}

    for caller in CALLERS:
        df_hifi_unique = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.hifi.uniques.info.tsv', sep='\t', header=[0])
        df_ont_unique = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.ont.uniques.info.tsv', sep='\t', header=[0])

        unique_svs.append((len(df_hifi_unique), TOOLMAP[caller], 'HiFi-Uniques'))
        unique_svs.append((len(df_ont_unique), TOOLMAP[caller], 'ONT-Uniques'))
        xticklabels.append(TOOLMAP[caller])

        read_caller_counts['hifi'].append(len(df_hifi_unique))
        read_caller_counts['ont'].append(len(df_ont_unique))

    for asm_caller in ASMCALLERS:
        df_hifi_unique = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{aligner}.hifi.uniques.info.tsv', sep='\t', header=[0])
        df_ont_unique = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{aligner}.ont.uniques.info.tsv', sep='\t', header=[0])

        unique_svs.append((len(df_hifi_unique), TOOLMAP[asm_caller], 'HiFi-Uniques'))
        unique_svs.append((len(df_ont_unique), TOOLMAP[asm_caller], 'ONT-Uniques'))
        xticklabels.append(TOOLMAP[asm_caller])

        asm_caller_counts['hifi'].append(len(df_hifi_unique))
        asm_caller_counts['ont'].append(len(df_ont_unique))

    print('Read-based calling === ')
    for plat, counts in read_caller_counts.items():
        print(f'{plat} median: {np.median(counts)} \tstd: {np.std(counts)}')

    print('Assembly-based calling === ')
    for plat, counts in asm_caller_counts.items():
        print(f'{plat} median: {np.median(counts)} \tstd: {np.std(counts)}')

def overview_dataset_repro(workdir, aligner):

    labels = [1, 2, 3, 4, 5, 6]
    supp_pcrt = []

    xticks = np.arange(len(labels))

    fig, axes = plt.subplots(1, 3, sharex='col', sharey='row', figsize=(7, 3))



    for col_idx, caller in enumerate(ASMCALLERS):

        for assembler in ASSEMBLER:
            merged_vcf = f'{workdir}/assm_dataset_repro/{caller}.{assembler}.jasmine.merged.vcf'
            assm_pcrt = []
            supp_dict, merged_total = get_survivor_supp(merged_vcf)

            for supp in labels:
                count = supp_dict[supp]
                supp_pcrt.append((count / merged_total, f'SUPP={supp}', 'Assembly-based'))
                assm_pcrt.append(count / merged_total)
                if supp == 6:
                    print(f'{caller}, {assembler}, {count / merged_total * 100}')

            axes[col_idx].plot(xticks, assm_pcrt, color=TOOLCOLORS[caller], marker=ASSEMBLERMK[assembler],
                         lw=2, label=assembler)


    for caller in CALLERS:
        merged_vcf = f'{workdir}/read_dataset_repro/{caller}.{aligner}.jasmine.merged.vcf'
        supp_dict, merged_total = get_survivor_supp(merged_vcf)
        read_pcrt = []
        for supp in labels:
            count = 0
            if supp in supp_dict:
                count = supp_dict[supp]
            supp_pcrt.append((count / merged_total, f'SUPP={supp}', 'Read-based'))
            read_pcrt.append(count / merged_total)
            if supp == 6:
                print(f'{caller}, {count / merged_total * 100}')

        axes[2].plot(xticks, read_pcrt, color=TOOLCOLORS[caller], marker='o', lw=2, label=TOOLMAP[caller])

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('Percent of strategy', fontsize=13)
        else:
            ax.set_ylabel('')
        if i == 2:
            ax.legend()
        elif i == 1:
            assembler_legends = [Line2D([0], [0], marker=mk, color=TOOLCOLORS['svimasm'], lw=2, markersize=6, label=asm) for asm, mk
                                 in ASSEMBLERMK.items()]
            ax.legend(handles=assembler_legends)
        else:
            assembler_legends = [Line2D([0], [0], marker=mk, color=TOOLCOLORS['pav'], lw=2, markersize=6, label=asm) for asm, mk
                in ASSEMBLERMK.items()]
            ax.legend(handles=assembler_legends)

        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, fontsize=13)


        ax.set_ylim(0, .6)
        ax.set_yticks(np.linspace(0, .6, 4))
        ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, .6, 4)], fontsize=12)


        ax.set_xlabel('')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig.tight_layout()

    hue_order = ['Assembly-based', 'Read-based']
    hue_color = [STRACOLORS[ele.split('-')[0]] for ele in hue_order]

    df_supp_pcrt = pd.DataFrame(supp_pcrt, columns=['pcrt', 'supp', 'strategy'])
    fig1, ax = plt.subplots(1, 1, figsize=(5, 3))
    bars = sns.barplot(data=df_supp_pcrt, x='supp', y='pcrt', hue='strategy', hue_order=hue_order,
                palette=hue_color, alpha=0.6, edgecolor='black', capsize=0.15, ci='sd', ax=ax)

    for bar in bars.patches:
        x = bar.get_x()
        width = bar.get_width()
        center = x + width / 2.

        bar.set_x(center - 0.3 / 2.)
        bar.set_width(0.3)
        bar.set_lw(2)

    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=13)
    ax.set_xlabel('')

    ax.set_ylim(0, .6)
    ax.set_yticks(np.linspace(0, .6, 4))
    ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, .6, 4)], fontsize=12)
    ax.set_ylabel('Percent of strategy', fontsize=13)
    ax.legend(title='')

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig1.tight_layout()

    plt.show()

    fig.savefig(f'{Figure2}/all_datasets_repro_by_callers.pdf')
    fig1.savefig(f'{Figure2}/all_datasets_repro_overview.pdf')


def overview_dataset_repro_highconf(workdir, aligner):

    labels = [1, 2, 3, 4, 5, 6]
    supp_pcrt = []

    for caller in ASMCALLERS:
        merged_vcf = f'{workdir}/assm_dataset_repro/highconf_regions/{caller}.{aligner}.jasmine.merged.vcf'
        supp_dict, merged_total = get_survivor_supp(merged_vcf)
        for supp in labels:
            count = supp_dict[supp]
            supp_pcrt.append((count / merged_total, f'SUPP={supp}', 'Assembly-based'))


    for caller in CALLERS:
        merged_vcf = f'{workdir}/read_dataset_repro/highconf_regions/{caller}.{aligner}.jasmine.merged.vcf'

        supp_dict, merged_total = get_survivor_supp(merged_vcf)
        for supp in labels:
            count = 0
            if supp in supp_dict:
                count = supp_dict[supp]
            supp_pcrt.append((count / merged_total, f'SUPP={supp}', 'Read-based'))

    df_supp_pcrt = pd.DataFrame(supp_pcrt, columns=['pcrt', 'supp', 'strategy'])
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    sns.boxplot(data=df_supp_pcrt, x='supp', y='pcrt', hue='strategy',ax=ax)
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels([f'SUPP={supp}' for supp in labels], fontsize=12, rotation=90)

    ax.set_ylim(0, 1)
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)])

    ax.legend(title='')
    ax.set_xlabel('')
    ax.set_ylabel('% of SVs', fontsize=13)

    fig.tight_layout()

    flags = ['hifi', 'ont']
    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))
    labels = [1, 2, 3]

    for col_idx, flag in enumerate(flags):
        this_ax = axes[col_idx]
        supp_pcrt = []
        for caller in ASMCALLERS:
            merged_vcf = f'{workdir}/assm_dataset_repro/{flag}/highconf_regions/{caller}.{aligner}.jasmine.merged.vcf'
            supp_dict, merged_total = get_survivor_supp(merged_vcf)
            for supp in labels:
                count = supp_dict[supp]
                supp_pcrt.append((count / merged_total, f'SUPP={supp}', 'Assembly-based'))

        for caller in CALLERS:
            merged_vcf = f'{workdir}/read_dataset_repro/{flag}/highconf_regions/{caller}.{aligner}.jasmine.merged.vcf'
            supp_dict, merged_total = get_survivor_supp(merged_vcf)

            for supp in labels:
                count = 0
                if supp in supp_dict:
                    count = supp_dict[supp]
                supp_pcrt.append((count / merged_total, f'SUPP={supp}', 'Read-based'))

        df_supp_pcrt = pd.DataFrame(supp_pcrt, columns=['pcrt', 'supp', 'strategy'])

        sns.boxplot(data=df_supp_pcrt, x='supp', y='pcrt', hue='strategy', ax=this_ax)

        if col_idx == 0:
            this_ax.set_ylabel('% of SVs', fontsize=13)
        else:
            this_ax.set_ylabel('')

        this_ax.set_ylim(0, 1)
        this_ax.set_yticks(np.linspace(0, 1, 5))
        this_ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)])
        this_ax.set_xlabel('')
        this_ax.legend(title='')

        this_ax.set_xticks(np.arange(len(labels)))
        this_ax.set_xticklabels([f'SUPP={ele}' for ele in labels], fontsize=12, rotation=90)

    fig1.tight_layout()
    plt.show()

    fig.savefig(f'{Figure2}/all_datasets_highconf_repro_overview.pdf')
    fig1.savefig(f'{Figure2}/all_datasets_highconf_repro_platform_overview.pdf')

def overview_datasets_unique(workdir, callers, aligner):

    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    platform_uniques = {'hifi': [[0 for i in range(len(callers))] for j in range(len(region_labels))],
                             'ont': [[0 for i in range(len(callers))] for j in range(len(region_labels))]}

    total_uniques = [0 for i in range(len(callers))]
    total_unique_list = []

    read_platform_uniques = {'hifi': [0 for i in range(len(CALLERS))], 'ont': [0 for i in range(len(CALLERS))]}

    assm_platform_uniques = {'hifi': [0 for i in range(len(ASMCALLERS))], 'ont': [0 for i in range(len(ASMCALLERS))]}

    for i, caller in enumerate(CALLERS):
        caller_idx = callers.index(caller)
        for plat in ['hifi', 'ont']:
            df_unique = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.{plat}.uniques.info.tsv', sep='\t', header=[0])
            unique_regions = {}
            total_uniques[caller_idx] += len(df_unique)
            total_unique_list.append((len(df_unique), 'Read-based', PLATMAP[plat], caller))
            for idx, row in df_unique.iterrows():
                if row['REGION_TYPE'] in unique_regions:
                    unique_regions[row['REGION_TYPE']] += 1
                else:
                    unique_regions[row['REGION_TYPE']] = 1

            for region, count in unique_regions.items():
                region_idx = region_labels.index(region)
                platform_uniques[plat][region_idx][caller_idx] += count
                read_platform_uniques[plat][i] += count


    for j, asm_caller in enumerate(ASMCALLERS):
        caller_idx = callers.index(asm_caller)
        for assembler in ASSEMBLER:
            for plat in ['hifi', 'ont']:
                df_unique = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{assembler}.{plat}.uniques.info.tsv', sep='\t', header=[0])
                unique_regions = {}
                total_uniques[caller_idx] += len(df_unique)
                total_unique_list.append((len(df_unique), 'Assembly-based', PLATMAP[plat], asm_caller))
                for idx, row in df_unique.iterrows():
                    if row['REGION_TYPE'] in unique_regions:
                        unique_regions[row['REGION_TYPE']] += 1
                    else:
                        unique_regions[row['REGION_TYPE']] = 1

                for region, count in unique_regions.items():
                    region_idx = region_labels.index(region)
                    platform_uniques[plat][region_idx][caller_idx] += count

                    assm_platform_uniques[plat][j] += count


    read_region_pcrt_plat = {'hifi': [], 'ont': []}
    assm_region_pcrt_plat = {'hifi': [], 'ont': []}
    region_pcrts = []
    platform_unique_pcrt = {'hifi': [[0 for i in range(len(callers))] for j in range(len(region_labels))],
                             'ont': [[0 for i in range(len(callers))] for j in range(len(region_labels))]}

    for plat, region_caller in platform_uniques.items():
        read_plat_counts = read_platform_uniques[plat]
        assm_plat_counts = assm_platform_uniques[plat]
        for i, caller_counts in enumerate(region_caller):
            for j, count in enumerate(caller_counts):
                if j < len(read_plat_counts):
                    print(f'{plat} {region_labels[i]}: {count / read_plat_counts[j]}')
                    region_pcrts.append((PLATMAP[plat], region_labels[i], count / read_plat_counts[j] * 100, 'Read-based'))
                    platform_unique_pcrt[plat][i][j] += count / read_plat_counts[j] * 100

                    if region_labels[i] == 'Simple Repeats':
                        read_region_pcrt_plat[plat].append(count / read_plat_counts[j] * 100)
                else:
                    print(f'{plat} {region_labels[i]}: {count / assm_plat_counts[j - len(read_plat_counts)]}')
                    region_pcrts.append((PLATMAP[plat], region_labels[i], count / assm_plat_counts[j - len(read_plat_counts)] * 100, 'Assembly-based'))
                    platform_unique_pcrt[plat][i][j] += count / assm_plat_counts[j - len(read_plat_counts)] * 100

                    if region_labels[i] == 'Simple Repeats':
                        assm_region_pcrt_plat[plat].append(count / assm_plat_counts[j - len(read_plat_counts)] * 100)


    print('Read-based =======')

    for plat, counts in read_platform_uniques.items():
        print(f'{plat}: {np.median(counts)}')
        print(f'Simple Repeats: {np.median(read_region_pcrt_plat[plat])}')


    print('Assembly-based =======')
    for plat, counts in assm_platform_uniques.items():
        print(f'{plat}: {np.median(counts)}')
        print(f'Simple Repeats: {np.median(assm_region_pcrt_plat[plat])}')


    df_total_uniques = pd.DataFrame(total_unique_list, columns=['count', 'stra', 'plat', 'caller'])
    fig, ax = plt.subplots(1, 1, figsize=(5, 2))
    plat_order= ['HiFi', 'ONT']
    plat_color = [PLATCOLORS[plat] for plat in plat_order]
    plat_legends = [Patch(facecolor=PLATCOLORS[plat], label=plat) for plat in plat_order]

    sns.boxplot(data=df_total_uniques, ax=ax, y='stra', x='count', hue='plat', hue_order=plat_order, palette=plat_color)


    ax.legend(handles=plat_legends)
    ax.set_xlim(0, 40000)
    ax.set_xticks(np.linspace(0, 40000, 5))
    ax.set_xticklabels([int(val) for val in np.linspace(0, 40, 5)], fontsize=13)
    ax.set_xlabel('Number of dataset unique SVs (x$10^3$)', fontsize=13)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)
    ax.set_ylabel('')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


    fig.tight_layout()

    barwidth = 0.3
    r1 = np.arange(len(callers))
    r2 = [r + barwidth + 0.05 for r in r1]

    plat_xticks = {'hifi': r1, 'ont': r2}

    fig1, ax = plt.subplots(1, 1, figsize=(6, 4))

    legends = [Patch(facecolor=REGIONCOLORS[region], edgecolor='black', label=region) for region in region_labels]
    legends.append(Patch(facecolor='white', edgecolor='black', hatch='///', label='ONT-Uniques'))
    legends.append(Patch(facecolor='white', edgecolor='black', label='HiFi-Uniques'))

    # for plat, unique_counts in platform_uniques.items():
    #     hatch = '///'
    #     if plat == 'hifi':
    #         hatch = ''
    #     for i, region in enumerate(region_labels):
    #         this_data = unique_counts[i]
    #         if i == 0:
    #             ax.bar(plat_xticks[plat], this_data, width=barwidth, color=REGIONCOLORS[region], edgecolor='black', hatch=hatch)
    #         else:
    #             bottoms = []
    #             for j in range(0, i):
    #                 bottoms.append(unique_counts[j])
    #             bottom_sum = [sum(x) for x in zip(*bottoms)]
    #             ax.bar(plat_xticks[plat], this_data, width=barwidth, bottom=bottom_sum, color=REGIONCOLORS[region], edgecolor='black', hatch=hatch)

    for plat, unique_counts in platform_unique_pcrt.items():
        hatch = '///'
        if plat == 'hifi':
            hatch = ''
        for i, region in enumerate(region_labels):
            this_data = unique_counts[i]
            if i == 0:
                ax.bar(plat_xticks[plat], this_data, width=barwidth, color=REGIONCOLORS[region], edgecolor='black', hatch=hatch)
            else:
                bottoms = []
                for j in range(0, i):
                    bottoms.append(unique_counts[j])
                bottom_sum = [sum(x) for x in zip(*bottoms)]
                ax.bar(plat_xticks[plat], this_data, width=barwidth, bottom=bottom_sum, color=REGIONCOLORS[region], edgecolor='black', hatch=hatch)

    ax.set_ylim(0, 100)
    ax.set_yticks(np.linspace(0, 100, 5))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
    ax.set_ylabel('Percent of SVs', fontsize=13)

    ax.legend(handles=legends, ncol=2)
    ax.set_xticks([r + (barwidth + 0.05) / 2 for r in r1])
    ax.set_xticklabels([TOOLMAP[caller] for caller in callers], fontsize=13, rotation=90)

    # ax1 = ax.twinx()
    # ax1.plot([r + (barwidth + 0.05) / 2 for r in r1], total_uniques, marker='X', lw=2, color='red')
    # ax1.set_ylim(0, 20000)
    # ax1.set_yticks(np.linspace(0, 20000, 5))
    # ax1.set_yticklabels([int(val) for val in np.linspace(0, 20, 5)], fontsize=12)
    # ax1.set_ylabel('Number of SVs (x$10^3$)', fontsize=13)

    fig1.tight_layout()

    df_region_pcrt = pd.DataFrame(region_pcrts, columns=['plat', 'region', 'pcrt', 'stra'])
    fig2, ax = plt.subplots(1, 1, figsize=(5, 3))
    hue_order = ['Assembly-based', 'Read-based']
    hue_color = [STRACOLORS[ele.split('-')[0]] for ele in hue_order]

    legends = [Patch(facecolor=STRACOLORS['Assembly'], label='Assembly-based'),
               Patch(facecolor=STRACOLORS['Read'], label='Read-based')]

    sns.boxplot(data=df_region_pcrt, y='region', x='pcrt', hue='stra', hue_order=hue_order, palette=hue_color, ax=ax)

    ax.set_xlim(0, 100)
    ax.set_xticks(np.linspace(0, 100, 5))
    ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
    ax.set_xlabel('Percent of dataset unique SVs', fontsize=13)

    ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylabel('')
    ax.legend(handles=legends, loc='lower right')

    fig2.tight_layout()
    plt.show()

    fig.savefig(f'{Figure2}/datasets_uniques_overview.pdf')
    fig1.savefig(f'{Figure2}/datasets_uniques_stackedplot_of_callers.pdf')
    fig2.savefig(f'{Figure2}/datasets_uniques_regions.pdf')

def datasets_unique_svsize(workdir, aligner):

    # unique_size_dict = {'Simple Repeats':[], 'Repeat Masked':[], 'Segment Dup':[], 'Unique':[]}
    unique_svtypes_in_regions = {'Simple Repeats':[], 'Repeat Masked':[], 'Segment Dup':[], 'Unique':[]}
    unique_svtypes = []

    unique_del_size = []
    unique_ins_size = []

    # plat_unique_size = {'hifi': unique_size_dict, 'ont': unique_size_dict}

    xticklabels = []
    for caller in CALLERS:
        xticklabels.append(TOOLMAP[caller])
        for plat in ['hifi', 'ont']:
            df_unique = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.{plat}.uniques.info.tsv', sep='\t', header=[0])
            svtypes_dict = {'Simple Repeats':{}, 'Repeat Masked':{}, 'Segment Dup': {}, 'Unique':{}}
            this_unique_types = {}
            for idx, row in df_unique.iterrows():
                svtype = row['TYPE_MATCH']
                # unique_size_dict[row['REGION_TYPE']].append((abs(int(row['SVLEN'])) + 1, TOOLMAP[caller], f'{PLATMAP[plat]}-Uniques'))
                svlen = abs(int(row['SVLEN']))

                if svlen < 50:
                    continue

                if svtype == 'DEL':
                    unique_del_size.append((svlen, TOOLMAP[caller], f'{PLATMAP[plat]}', 'Read-based'))
                elif svtype == 'INS':
                    unique_ins_size.append((svlen, TOOLMAP[caller], f'{PLATMAP[plat]}', 'Read-based'))

                if svtype in this_unique_types:
                    this_unique_types[svtype] += 1
                else:
                    this_unique_types[svtype] = 1


                if svtype in svtypes_dict[row['REGION_TYPE']]:
                    svtypes_dict[row['REGION_TYPE']][svtype] += 1
                else:
                    svtypes_dict[row['REGION_TYPE']][svtype] = 1

            for region, count_dict in svtypes_dict.items():
                for svtype, count in count_dict.items():
                    unique_svtypes_in_regions[region].append((svtype, count, TOOLMAP[caller], f'{PLATMAP[plat]}'))

            for svtype, count in this_unique_types.items():
                unique_svtypes.append((svtype, count / len(df_unique) * 100, f'{PLATMAP[plat]}', TOOLMAP[caller]))


    for asm_caller in ASMCALLERS:
        xticklabels.append(TOOLMAP[asm_caller])
        for assembler in ASSEMBLER:
            for plat in ['hifi', 'ont']:
                df_unique = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{assembler}.{plat}.uniques.info.tsv', sep='\t', header=[0])
                svtypes_dict = {'Simple Repeats': {}, 'Repeat Masked': {}, 'Segment Dup': {}, 'Unique': {}}
                this_unique_types = {}
                for idx, row in df_unique.iterrows():
                    svtype = row['TYPE_MATCH']
                    svlen = abs(int(row['SVLEN']))
                    if svlen < 50:
                        continue

                    # unique_size_dict[row['REGION_TYPE']].append((abs(int(row['SVLEN'])) + 1, TOOLMAP[asm_caller], f'{PLATMAP[plat]}-Uniques'))
                    if svtype == 'DEL':
                        unique_del_size.append((svlen, TOOLMAP[asm_caller], f'{PLATMAP[plat]}', 'Assembly-based'))
                    elif svtype == 'INS':

                        unique_ins_size.append((svlen, TOOLMAP[asm_caller], f'{PLATMAP[plat]}', 'Assembly-based'))

                    # if svtype in this_unique_types:
                    #     this_unique_types[svtype] += 1
                    # else:
                    #     this_unique_types[svtype] = 1

                    if svtype in svtypes_dict[row['REGION_TYPE']]:
                        svtypes_dict[row['REGION_TYPE']][svtype] += 1
                    else:
                        svtypes_dict[row['REGION_TYPE']][svtype] = 1

                for region, count_dict in svtypes_dict.items():
                    for svtype, count in count_dict.items():
                        unique_svtypes_in_regions[region].append((svtype, count, TOOLMAP[asm_caller], f'{PLATMAP[plat]}'))

                # for svtype, count in this_unique_types.items():
                #     unique_svtypes.append((svtype, count / len(df_unique) * 100, f'{PLATMAP[plat]}', TOOLMAP[asm_caller]))

    fig, axes = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(6, 4))
    df_unqie_ins = pd.DataFrame(unique_ins_size, columns=['svlen', 'caller', 'plat', 'stra'])
    df_unqie_del = pd.DataFrame(unique_del_size, columns=['svlen', 'caller', 'plat', 'stra'])
    stra_hue_order = ['Read-based', 'Assembly-based']
    stra_hue_color = [STRACOLORS[stra.split('-')[0]] for stra in stra_hue_order]

    legends = [Patch(facecolor=STRACOLORS[stra.split('-')[0]], label=stra) for stra in stra_hue_order]
    sns.histplot(data=df_unqie_ins, x='svlen', hue='stra', hue_order=stra_hue_order, palette=stra_hue_color,
                 log_scale=[False, True], bins=100, ax=axes[0])

    # ax1 = axes[0].inset_axes([.35, .35, .4, .45])
    # sns.histplot(data=df_unqie_ins[df_unqie_ins['svlen'] < 200], x='svlen', hue='plat', hue_order=plat_hue_order, log_scale=[False, True], palette=plat_hue_color, bins=50, ax=ax1)
    # ax1.set_ylabel('')
    # # ax1.set_yscale('log')
    # ax1.set_xlabel('')
    # ax1.legend([])
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    # ax1.tick_params(which='minor', length=4, width=1)
    #
    # ax1.set_xticks(np.linspace(0, 200, 3))
    # ax1.set_xticklabels([0, 100, 200])

    sns.histplot(data=df_unqie_del, x='svlen', hue='stra', hue_order=stra_hue_order, palette=stra_hue_color,
                 log_scale=[False, True], bins=100, ax=axes[1])

    # ax1 = axes[1].inset_axes([.35, .35, .4, .45])
    # sns.histplot(data=df_unqie_del[df_unqie_del['svlen'] < 200], x='svlen', hue='plat', hue_order=plat_hue_order, log_scale=[False, True], palette=plat_hue_color, bins=50, ax=ax1)
    # ax1.set_ylabel('')
    # ax1.set_xlabel('')
    # ax1.legend([])
    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)
    # ax1.set_xticks(np.linspace(0, 200, 3))
    # ax1.set_xticklabels([0, 100, 200])
    # ax1.tick_params(which='minor', length=4, width=1)

    # sns.ecdfplot(data=df_unqie_size, x='svlen', hue='plat', ax=axes[1])
    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('INS (x$10^3$)', fontsize=13)
        else:
            ax.set_ylabel('DEL (x$10^3$)', fontsize=13)

        ax.set_xlabel('')
        ax.legend(handles=legends)
        # ax.set_ylim(0, 8000)
        # ax.set_yticks(np.linspace(0, 8000, 5))
        # ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 8, 5)], fontsize=12)
        # ax.set_xlim(0, 100000)
        ax.set_xticks(np.linspace(0, 100000, 6))
        ax.set_xticklabels([int(val) for val in np.linspace(0, 100, 6)], fontsize=12)
        ax.set_xlabel('SV size (kbp)', fontsize=13)


        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{Figure2}/platform_uniques_svlen.pdf')

    plt.show()

    # hue_order = ['DEL', 'INS', 'DUP', 'INV', 'BND', 'DUP:TANDEM']
    # hue_color = [SVTYPECOLORS[svtype] for svtype in hue_order]

    # fig1, axes =  plt.subplots(2, 1, sharex='col', sharey='row', figsize=(6, 4))
    # df_unique_types = pd.DataFrame(unique_svtypes, columns=['svtype', 'pcrt', 'plat', 'caller'])
    #
    # sns.barplot(data=df_unique_types[df_unique_types['plat'] == 'HiFi'], x='caller', y='pcrt', hue='svtype', hue_order=hue_order, palette=hue_color, ax=axes[0])
    # sns.barplot(data=df_unique_types[df_unique_types['plat'] == 'ONT'], x='caller', y='pcrt', hue='svtype',
    #             hue_order=hue_order, palette=hue_color, ax=axes[1])
    #
    # for i, ax in enumerate(axes):
    #     ax.set_xlabel('')
    #     if i == 0:
    #         ax.set_ylabel('Percent of SVs', fontsize=13)
    #     else:
    #         ax.set_ylabel('')
    #     ax.set_ylim(0, 80)
    #     ax.set_yticks(np.linspace(0, 80, 5))
    #     ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 80, 5)], fontsize=12)
    #     ax.spines['top'].set_visible(False)
    #     ax.spines['right'].set_visible(False)
    #
    #     ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=13)
    #
    #     ax.legend(title='', ncol=3)
    #
    # fig1.tight_layout()
    # fig1.savefig(f'{Figure2}/datasets_uniques_svtypes.pdf')

    # region_labels = ['Simple Repeats']
    # for col_idx, region in enumerate(region_labels):
    #     fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    #     df_regions = pd.DataFrame(unique_size_dict[region], columns=['svlen', 'caller', 'plat'])
    #     sns.histplot(data=df_regions, x='svlen', hue='plat', ax=ax, bins=50)
    #     sns.ecdfplot(data=df_regions, x='svlen', hue='plat', ax=ax)
    #     ax.set_xlabel('')
    #     ax.set_yticks(np.linspace(0, 8000, 5))
    #     ax.set_yticklabels([int(val) for val in np.linspace(0, 8, 5)])
    #     ax.set_ylabel('# of SVs (x1000)')
    #     fig.tight_layout()
        # fig.savefig(f'{Figure2}/platform_uniques_simrep_svlen.pdf')

    # hue_order = ['DEL', 'INS', 'DUP', 'INV', 'BND', 'DUP:TANDEM']
    # hue_color = [SVTYPECOLORS[svtype] for svtype in hue_order]
    #
    # fig1, ax1 = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 3))
    # df_svtypes_regions = pd.DataFrame(unique_svtypes_in_regions['Simple Repeats'], columns=['svtype', 'count', 'caller', 'plat'])
    # sns.barplot(data=df_svtypes_regions[df_svtypes_regions['plat'] == 'HiFi-Uniques'], x='caller', y='count', hue='svtype', hue_order=hue_order, palette=hue_color, ax=ax1[0])
    # ax1[0].set_xlabel('')
    # ax1[0].set_ylabel('# of SVs (x1000)')
    # ax1[0].set_title('HiFi-Unique')
    # ax1[0].set_ylim(0, 5000)
    # ax1[0].set_yticks(np.linspace(0, 5000, 6))
    # ax1[0].set_yticklabels([int(val) for val in np.linspace(0, 5, 6)])
    # ax1[0].legend(ncol=2)
    #
    # df_svtypes_regions = pd.DataFrame(unique_svtypes_in_regions['Simple Repeats'], columns=['svtype', 'count', 'caller', 'plat'])
    # sns.barplot(data=df_svtypes_regions[df_svtypes_regions['plat'] == 'ONT-Uniques'], x='caller', y='count', hue='svtype', hue_order=hue_order, palette=hue_color, ax=ax1[1])
    # ax1[1].set_xlabel('')
    # ax1[1].set_ylabel('')
    # ax1[1].set_title('ONT-Unique')
    # ax1[1].legend('')
    # fig1.tight_layout()

    # fig1.savefig(f'{Figure2}/platform_uniques_simrep_svtypes.pdf')
    # fig2.savefig(f'{Figure2}/ont_uniques_simrep_svtypes.pdf')


def read_dataset_unique_svtype(workdir, aligner):

    svtypes = ['INS', 'DEL', 'DUP', 'INV']

    plt.figure(figsize=(5, 4))

    for i, plat in enumerate(['hifi', 'ont']):
        this_plat_dict = {'group': CALLERS,
                            'INS': [0 for i in range(len(CALLERS))],
                            'DEL': [0 for i in range(len(CALLERS))],
                            'INV': [0 for i in range(len(CALLERS))],
                            'DUP': [0 for i in range(len(CALLERS))]}

        this_ax = plt.subplot(1, 2, i + 1, polar=True)
        this_ax.set_title(PLATMAP[plat], fontsize=13)

        for j, caller in enumerate(CALLERS):

            df_unique = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.{plat}.uniques.info.tsv', sep='\t', header=[0])
            this_unique_types = {}
            for idx, row in df_unique.iterrows():
                svtype = row['TYPE_MATCH']
                # unique_size_dict[row['REGION_TYPE']].append((abs(int(row['SVLEN'])) + 1, TOOLMAP[caller], f'{PLATMAP[plat]}-Uniques'))

                if svtype in this_unique_types:
                    this_unique_types[svtype] += 1
                else:
                    this_unique_types[svtype] = 1

            for svtype, count in this_unique_types.items():
                if svtype in svtypes:
                    this_pcrt = count / len(df_unique) * 100
                    this_plat_dict[svtype][j] += this_pcrt

        radar_plot(pd.DataFrame(this_plat_dict), this_ax, CALLERS)

    plt.tight_layout()

    plt.savefig(f'{Figure2}/read_datasets_uniques_svtypes_radar.pdf')
    plt.show()


def assm_dataset_unique_svtype(workdir):

    unique_svtype_pcrt = []
    for plat in ['hifi', 'ont']:
        for j, caller in enumerate(ASMCALLERS):
            for assembler in ASSEMBLER:
                df_unique = pd.read_csv(f'{workdir}/assm_dataset_repro/{caller}.{assembler}.{plat}.uniques.info.tsv', sep='\t', header=[0])
                this_unique_types = {'INS': 0, 'DEL': 0}
                total_uniques = len(df_unique)
                for idx, row in df_unique.iterrows():
                    svtype = row['TYPE_MATCH']
                    # unique_size_dict[row['REGION_TYPE']].append((abs(int(row['SVLEN'])) + 1, TOOLMAP[caller], f'{PLATMAP[plat]}-Uniques'))

                    if svtype in this_unique_types:
                        this_unique_types[svtype] += 1
                    else:
                        this_unique_types[svtype] = 1
                unique_svtype_pcrt.append(('INS', PLATMAP[plat], this_unique_types['INS'] / total_uniques * 100, assembler, caller))
                unique_svtype_pcrt.append(('DEL', PLATMAP[plat], this_unique_types['DEL'] / total_uniques * 100, assembler, caller))

    df_unique_svtype = pd.DataFrame(unique_svtype_pcrt, columns=['svtype', 'plat', 'pcrt', 'assembler', 'caller'])
    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 3))

    hue_order = ['INS', 'DEL']
    hue_color = [SVTYPECOLORS[svtype] for svtype in hue_order]

    sns.barplot(data=df_unique_svtype[df_unique_svtype['plat'] == 'HiFi'], x='caller', y='pcrt',
                hue='svtype', hue_order=hue_order, palette=hue_color, capsize=0.15, ci='sd', ax=axes[0])
    axes[0].set_title('HiFi', fontsize=13)
    sns.barplot(data=df_unique_svtype[df_unique_svtype['plat'] == 'ONT'], x='caller', y='pcrt',
                hue='svtype', hue_order=hue_order, palette=hue_color, capsize=0.15, ci='sd', ax=axes[1])
    axes[1].set_title('ONT', fontsize=13)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of SVs', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_xlabel('')
        ax.set_xticks(np.arange(2))
        ax.set_xticklabels(['PAV', 'SVIM-asm'], fontsize=13)
        ax.legend(title='')
        ax.set_ylim(0, 80)
        ax.set_yticks(np.linspace(0, 80, 5))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 80, 5)], fontsize=12)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig.tight_layout()

    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 3))


    hifi_unique_svtypes = df_unique_svtype[df_unique_svtype['plat'] == 'HiFi']
    ont_unique_svtypes = df_unique_svtype[df_unique_svtype['plat'] == 'ONT']

    sns.stripplot(data=hifi_unique_svtypes[hifi_unique_svtypes['svtype'] == 'INS'], x='caller', y='pcrt',
                hue='assembler', size=8, ax=axes[0])
    axes[0].set_title('HiFi', fontsize=13)

    sns.stripplot(data=ont_unique_svtypes[ont_unique_svtypes['svtype'] == 'INS'], x='caller', y='pcrt',
                hue='assembler', size=8, ax=axes[1])

    axes[1].set_title('ONT', fontsize=13)

    for i, ax in enumerate(axes):

        if i == 0:
            ax.set_ylabel('% of INS', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_xlabel('')
        ax.set_xticks(np.arange(2))
        ax.set_xticklabels(['PAV', 'SVIM-asm'], fontsize=13)
        ax.legend('')
        ax.set_ylim(40, 80)
        ax.set_yticks(np.linspace(40, 80, 3))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(40, 80, 3)], fontsize=12)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig1.tight_layout()

    fig.savefig(f'{Figure2}/assm_datasets_uniques_svtypes.pdf')
    fig1.savefig(f'{Figure2}/assm_datasets_uniques_ins.pdf')
    plt.show()

def radar_plot(df, ax, groups):
    categories = list(df)[1:]
    N = len(categories)

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Initialise the spider plot
    # plt.figure(figsize=(3, 3))
    # ax = plt.subplot(111, polar=True)

    # If you want the first axis to be on top:
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("NW")
    # Draw one axe per variable + add labels
    plt.xticks(angles[:-1], categories, fontsize=12)

    # Draw ylabels
    ax.set_rlabel_position(180)
    plt.yticks([0, 40, 80], ["0", "40%", "80%"], size=10)
    plt.ylim(0, 80)

    # ------- PART 2: Add plots

    # Plot each individual = each line of the data
    # I don't make a loop, because plotting more than 3 groups makes the chart unreadable

    for i, group_name in enumerate(groups):
        values = df.loc[i].drop('group').values.flatten().tolist()
        values += values[:1]
        ls = '-'
        if group_name in ASMCALLERS:
            ls = '--'

        ax.plot(angles, values, linewidth=2, color=TOOLCOLORS[group_name], ls=ls, label=TOOLMAP[group_name])
        # ax.fill(angles, values, TOOLCOLORS[group_name], alpha=0.1)
        ax.scatter(angles, values, s=20, color=TOOLCOLORS[group_name])



def dataset_repro_at_regions(workdir, aligner):

    flags = ['all', 'hifi', 'ont']

    for col_idx, flag in enumerate(flags):

        fig, this_ax = plt.subplots(1, 1, figsize=(6, 4))
        labels = [1, 2, 3, 4, 5, 6]

        if flag == 'hifi' or flag == 'ont':
            fig, this_ax = plt.subplots(1, 1, figsize=(4, 4))
            labels = [1, 2, 3]

        this_ax.set_ylabel('% of SVs', fontsize=13)

        this_ax.set_ylim(0, 1)
        this_ax.set_yticks(np.linspace(0, 1, 5))
        this_ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)])

        this_ax.set_xticks(np.arange(len(labels)))
        this_ax.set_xticklabels([f'SUPP={ele}' for ele in labels], fontsize=12, rotation=90)

        for caller in ASMCALLERS:
            merged_vcf = f'{workdir}/assm_dataset_repro/highconf_regions/{caller}.{aligner}.jasmine.merged.vcf'
            if flag == 'hifi' or flag == 'ont':
                merged_vcf = f'{workdir}/assm_dataset_repro/{flag}/highconf_regions/{caller}.{aligner}.jasmine.merged.vcf'
            supp_dict, merged_total = get_survivor_supp(merged_vcf)

            data = []
            for supp in labels:
                count = supp_dict[supp]
                data.append(count / merged_total)

            this_ax.plot(np.arange(len(labels)), data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, label=TOOLMAP[caller])
            this_ax.legend()

        for caller in CALLERS:
            merged_vcf = f'{workdir}/read_dataset_repro/highconf_regions/{caller}.{aligner}.jasmine.merged.vcf'
            if flag != 'all':
                merged_vcf = f'{workdir}/read_dataset_repro/{flag}/highconf_regions/{caller}.{aligner}.jasmine.merged.vcf'

            supp_dict, merged_total = get_survivor_supp(merged_vcf)

            data = []
            for supp in labels:
                count = 0
                if supp in supp_dict:
                    count = supp_dict[supp]
                data.append(count / merged_total)

            this_ax.plot(np.arange(len(labels)), data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, label=TOOLMAP[caller])
            this_ax.legend()

        plt.tight_layout()
        plt.show()
        fig.savefig(f'{Figure2}/{flag}_datasets_repro_of_callers_at_highconf.pdf')


'''
def dataset_repro_of_callers(workdir, aligner):

    flags = ['hifi', 'ont']

    fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    labels = [1, 2, 3, 4, 5, 6]
    xticks = np.arange(len(labels))
    for caller in ASMCALLERS:
        merged_vcf = f'{workdir}/assm_dataset_repro/{caller}.{aligner}.jasmine.merged.vcf'
        supp_dict, merged_total = get_survivor_supp(merged_vcf)

        data = []
        for supp in labels:
            count = supp_dict[supp]
            data.append(count / merged_total)

        ax.plot([r - 0.2 for r in xticks], data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, label=TOOLMAP[caller])
        ax.legend()

    for caller in CALLERS:
        merged_vcf = f'{workdir}/read_dataset_repro/{caller}.{aligner}.jasmine.merged.vcf'
        supp_dict, merged_total = get_survivor_supp(merged_vcf)

        data = []
        for supp in labels:
            count = 0
            if supp in supp_dict:
                count = supp_dict[supp]
            data.append(count / merged_total)

        ax.plot([r + 0.2 for r in xticks], data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, label=TOOLMAP[caller])
        ax.legend()

    ax.set_ylabel('Percent of SVs', fontsize=13)
    ax.set_ylim(0, 1)
    ax.set_yticks(np.linspace(0, 1, 5))
    ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)], fontsize=12)

    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels([f'SUPP={ele}' for ele in labels], fontsize=13, rotation=90)
    fig.tight_layout()
    fig.savefig(f'{Figure2}/all_datasets_repro_of_callers.pdf')

    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(5, 3))
    for col_idx, flag in enumerate(flags):
        this_ax = axes[col_idx]

        labels = [1, 2, 3]

        if col_idx == 0:
            this_ax.set_ylabel('Percent of SVs', fontsize=13)
        else:
            this_ax.set_ylabel('')

        this_ax.set_ylim(0, 1)
        this_ax.set_yticks(np.linspace(0, 1, 5))
        this_ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)], fontsize=12)

        this_ax.set_xticks(np.arange(len(labels)))
        this_ax.set_xticklabels([f'SUPP={ele}' for ele in labels], fontsize=13, rotation=90)


        for caller in ASMCALLERS:
            merged_vcf = f'{workdir}/assm_dataset_repro/{caller}.{aligner}.jasmine.merged.vcf'
            if flag == 'hifi' or flag == 'ont':
                merged_vcf = f'{workdir}/assm_dataset_repro/{flag}/{caller}.{aligner}.jasmine.merged.vcf'
            supp_dict, merged_total = get_survivor_supp(merged_vcf)

            data = []
            for supp in labels:
                count = supp_dict[supp]
                data.append(count / merged_total)

            this_ax.plot(np.arange(len(labels)), data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller],
                         lw=2, label=TOOLMAP[caller])
            this_ax.legend()

        for caller in CALLERS:
            merged_vcf = f'{workdir}/read_dataset_repro/{caller}.{aligner}.jasmine.merged.vcf'
            if flag != 'all':
                merged_vcf = f'{workdir}/read_dataset_repro/{flag}/{caller}.{aligner}.jasmine.merged.vcf'

            supp_dict, merged_total = get_survivor_supp(merged_vcf)

            data = []
            for supp in labels:
                count = 0
                if supp in supp_dict:
                    count = supp_dict[supp]
                data.append(count / merged_total)

            this_ax.plot(np.arange(len(labels)), data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller],
                         lw=2, label=TOOLMAP[caller])
            this_ax.legend()


    fig1.tight_layout()
    plt.show()
    fig1.savefig(f'{Figure2}/datasets_repro_of_callers_at_platforms.pdf')
'''