#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/5/27

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
plt.rcParams["xtick.minor.width"] = 2
plt.rcParams["xtick.labelsize"] = 11

plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.minor.width"] = 2
plt.rcParams["ytick.labelsize"] = 11

plt.rcParams["axes.linewidth"] = 2


def overview_aligner_repro(workdir, datasets):
    total_counts = []

    datasets_dict = {'HiFi': ['hifi_10kb', 'hifi_15kb', 'hifi_18kb'],
                     'ONT': ['ont_9kb', 'ont_19kb', 'ont_30kb']}

    total_count_of_supp = []

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))

    for caller in CALLERS:

        plat_supp_pcrt1 = {'HiFi': [0, 0, 0], 'ONT': [0, 0 ,0]}
        plat_supp_pcrt2 = {'HiFi': [0, 0, 0], 'ONT': [0, 0, 0]}

        for dataset_idx, dataset in enumerate(datasets):
            plat = 'HiFi'
            if 'ont' in dataset:
                plat = 'ONT'

            this_dataset_index = datasets_dict[plat].index(dataset)

            merged_vcf = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.{dataset}.jasmine.merged.vcf'
            supp_dict, merged_total = get_survivor_supp(merged_vcf)
            total_counts.append((merged_total, plat, dataset, TOOLMAP[caller]))

            for supp, count in supp_dict.items():
                # pcrt_list.append((count / merged_total * 100, plat, caller, f'SUPP={supp}'))
                # total_count_of_supp.append((f'SUPP={supp}', count, plat, caller))
                if supp == 1:
                    plat_supp_pcrt1[plat][this_dataset_index] += count / merged_total * 100
                    total_count_of_supp.append((f'SUPP={supp}', count, plat, caller))
                elif supp == 4:
                    plat_supp_pcrt2[plat][this_dataset_index] += count / merged_total * 100
                    total_count_of_supp.append((f'SUPP={supp}', count, plat, caller))


        ax.scatter(plat_supp_pcrt1['HiFi'], plat_supp_pcrt1['ONT'], color='#416e9f', edgecolor='black', marker=TOOLMARKERS[caller], s=120)
        ax.scatter(plat_supp_pcrt2['HiFi'], plat_supp_pcrt2['ONT'], color='#599261', edgecolor='black', marker=TOOLMARKERS[caller], s=120)

    ax.set_xlabel('HiFi datasets', fontsize=13)
    ax.set_xlim(20, 60, 5)
    ax.set_xticks(np.linspace(20, 60, 5))
    ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(20, 60, 5)], fontsize=12)

    ax.set_ylabel('ONT datasets', fontsize=13)
    ax.set_ylim(20, 60, 5)
    ax.set_yticks(np.linspace(20, 60, 5))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(20, 60, 5)], fontsize=12)


    all_legends = [Patch(facecolor='#416e9f', label='SUPP=1'), Patch(facecolor='#599261', label='SUPP=4')]
    caller_legends = [Line2D([0], [0], marker=TOOLMARKERS[caller], mfc='white', mec='black', color='w',
                             markersize=10, label=TOOLMAP[caller]) for caller in CALLERS]

    all_legends.extend(caller_legends)
    ax.legend(handles=all_legends)

    ax.plot([20, 60], [20, 60], color='#c96150', ls='--')
    fig.tight_layout()
    fig.savefig(f'{Figure6}/aligner_repro_pcrt.pdf')

    plt.show()

def aligner_concordant_bpstd(workdir, datasets, aligners):

    shift_labels = ['0', '0,10', '10,50', '>50']
    svsize_list = []

    svsize_plat = {'HiFi': 0, 'ONT': 0}
    leftbp_pcrt_list = []

    for platform_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for caller in CALLERS:
            left_shift_label = {shift: 0 for shift in shift_labels}
            match_info_out = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-concordant.tsv'
            df_matched = pd.read_csv(match_info_out, sep='\t', header=[0])
            this_total = 0

            for idx, row in df_matched.iterrows():
                start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                if int(row['SUPP']) == 4:
                    this_total += 1
                    svsize_list.append((svlen + 1, f'{plat}-Unique'))
                    if start_std == 0:
                        left_shift_label['0'] += 1
                    elif start_std <= 10:
                        left_shift_label['0,10'] += 1
                    elif start_std <= 50:
                        left_shift_label['10,50'] += 1
                    else:
                        left_shift_label['>50'] += 1

                    if svlen > 10000:
                        svsize_plat[plat] += 1

            for shift in shift_labels:
                leftbp_pcrt_list.append((plat, PLATMAP[dataset], TOOLMAP[caller], left_shift_label[shift] / this_total * 100, shift))


    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 4))
    df_leftbp = pd.DataFrame(leftbp_pcrt_list, columns=['plat', 'dataset', 'caller', 'pcrt', 'shift'])
    plat_order = ['HiFi', 'ONT']
    plat_color = [PLATCOLORS[plat] for plat in plat_order]
    plat_legends = [Patch(facecolor=PLATCOLORS[plat], label=plat) for plat in plat_order]
    sns.stripplot(data=df_leftbp, x='shift', y='pcrt', hue='dataset', palette="Dark2", size=8, ax=axes[0])

    axes[0].set_xticks(np.arange(len(shift_labels)))
    axes[0].set_xticklabels(['0', '0~10', '10~50', '>50'], fontsize=13)
    axes[0].legend(title='')

    # shift_order = ['0', '0,10', '10,50', '>50']
    # shift_color = [SHIFTCOLORDICT[shift] for shift in shift_order]
    # shift_legend = [Patch(facecolor=SHIFTCOLORDICT['0'], label='0'), Patch(facecolor=SHIFTCOLORDICT['0,10'], label='0~10'),
    #                 Patch(facecolor=SHIFTCOLORDICT['10,50'], label='10~50'), Patch(facecolor=SHIFTCOLORDICT['>50'], label='>50')]

    caller_color = [TOOLCOLORS[caller] for caller in CALLERS]
    # caller_legend = [Patch(facecolor=TOOLCOLORS[caller], label=TOOLMAP[caller]) for caller in CALLERS]


    # sns.barplot(data=df_leftbp, x='shift', y='pcrt', hue='caller', hue_order=shift_order, palette=shift_color, capsize=.05, ax=axes[1])
    sns.stripplot(data=df_leftbp, x='shift', y='pcrt', hue='caller', palette=caller_color, size=8, ax=axes[1])
    # sns.stripplot(data=df_leftbp, x='shift', y='pcrt', hue='aligner', palette=aligner_color, size=8, ax=axes[1])
    axes[1].set_xticklabels(['0', '0~10', '10~50', '>50'], fontsize=13)
    axes[1].legend(title='')

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('Percent of concordant SVs', fontsize=13)
            ax.set_ylim(0, 100)
            ax.set_yticks(np.linspace(0, 100, 5))
            ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)

        else:
            ax.set_ylabel('')

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)
        ax.set_xlabel('')

    fig.tight_layout()
    plt.show()
    fig.savefig(f'{Figure6}/aligner_repro_bpstd.pdf')



def aligner_unique_features(workdir, datasets, aligners):
    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    aligner_unique_counts = []
    aligner_unique_region_pcrt = []

    aligner_unique_info = []

    aligner_uniques = {aligner: [] for aligner in aligners}
    aligner_uniques_region_pcrt = {aligner: {region: [] for region in regions} for aligner in aligners}

    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for j, caller in enumerate(CALLERS):

            unique_info_out = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-unique.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])
            aligner_unique_dict = {aligner: 0 for aligner in aligners}
            aligner_unique_region_dict = {aligner: [0 for i in range(len(regions))] for aligner in aligners}

            for idx, row in df_unique.iterrows():
                svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['CALLER'], row['REGION_TYPE']
                region_idx = regions.index(sv_region)
                aligner_unique_dict[aligner] += 1
                aligner_unique_region_dict[aligner][region_idx] += 1
                if svlen >= 50:
                    aligner_unique_info.append((svlen, aligner, sv_region, svtype, plat))

            for aligner, count in aligner_unique_dict.items():
                aligner_unique_counts.append((count, plat, PLATMAP[dataset], caller, aligner))
                aligner_uniques[aligner].append(count)
                for k, val in enumerate(aligner_unique_region_dict[aligner]):
                    aligner_unique_region_pcrt.append((val / count * 100, regions[k], plat, caller, aligner))
                    aligner_uniques_region_pcrt[aligner][regions[k]].append(val)


    df_unique_counts = pd.DataFrame(aligner_unique_counts, columns=['count', 'plat', 'dataset', 'caller', 'aligner'])
    # df_unique_pcrt_region = pd.DataFrame(aligner_unique_region_pcrt, columns=['count', 'region', 'plat', 'caller', 'aligner'])

    for aligner, counts in aligner_uniques.items():
        print(f'{aligner} median: {np.median(counts)} avg: {np.mean(counts)}')

    df_aligner_unique = pd.DataFrame(aligner_unique_info, columns=['svlen', 'aligner', 'region', 'svtype', 'plat'])
    fig, axes = plt.subplots(2, 1, sharex='col', sharey='row', figsize=(5, 4))

    aligner_color = [ALIGNERCOLOR[aligner] for aligner in aligners]
    aligner_legends = [Patch(facecolor=ALIGNERCOLOR[aligner], label=aligner) for aligner in aligners]

    sns.histplot(data=df_aligner_unique[df_aligner_unique['plat'] == 'HiFi'], x='svlen', hue='aligner', bins=50,
                 hue_order=aligners, palette=aligner_color, log_scale=True, ax=axes[0])
    axes[0].set_ylabel('HiFi (x$10^3$)', fontsize=13)

    sns.histplot(data=df_aligner_unique[df_aligner_unique['plat'] == 'ONT'], x='svlen', hue='aligner', bins=50,
                 hue_order=aligners, palette=aligner_color, log_scale=True, ax=axes[1])

    axes[1].set_ylabel('ONT (x$10^3$)', fontsize=13)

    for ax in axes:
        ax.legend(handles=aligner_legends)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel('')
        ax.set_ylim(0, 4000)
        ax.set_yticks(np.linspace(0, 4000, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 4, 5)], fontsize=12)


    fig.tight_layout()
    fig.savefig(f'{Figure6}/aligner_unique_svsize.pdf')

    # colors = [REGIONCOLORS[region] for region in regions]
    # for fig_idx, aligner in enumerate(aligners):
    #     fig, ax = plt.subplots(1, 1, figsize=(3, 2))
    #     this_aligner_avg = []
    #     for region, counts in aligner_uniques_region_pcrt[aligner].items():
    #         this_aligner_avg.append(np.mean(counts))
    #     ax.pie(this_aligner_avg, autopct='%1.1f%%', colors=colors, startangle=90, pctdistance=0.7,
    #            wedgeprops=dict(width=0.7, edgecolor='w'))
    #
    #     fig.tight_layout()
    #     fig.savefig(f'{Figure6}/{aligner}_unique_region.pdf')

    # plat_order = ['HiFi', 'ONT']
    # plat_color = [PLATCOLORS[plat] for plat in plat_order]
    # plat_legends = [Patch(facecolor=PLATCOLORS[plat], label=plat) for plat in plat_order]

    fig1, ax = plt.subplots(1, 1, figsize=(4, 4))
    df_unique_region = pd.DataFrame(aligner_unique_region_pcrt, columns=['pcrt', 'region', 'plat', 'caller', 'aligner'])
    sns.barplot(data=df_unique_region, y='region', x='pcrt', hue='aligner', hue_order=aligners, palette=aligner_color, ax=ax)
    ax.set_xlabel('Percent of SVs', fontsize=13)
    ax.set_xlim(0, 100)
    ax.set_xticks(np.linspace(0, 100, 5))
    ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)

    ax.set_ylabel('')
    ax.set_yticks(np.arange(len(regions)))
    ax.set_yticklabels(regions, fontsize=13)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.grid(axis='x', ls='--', color='grey', lw=2)
    ax.legend(handles=aligner_legends)

    fig1.tight_layout()
    fig1.savefig(f'{Figure6}/aligner_unique_region.pdf')

    fig2, ax = plt.subplots(1, 1, figsize=(4, 4))
    sns.stripplot(data=df_unique_counts, x='aligner', y='count', hue='dataset', palette="Dark2", size=8, ax=ax)
    ax.set_xlabel('')
    ax.set_xticks(np.arange(len(aligners)))
    ax.set_xticklabels(aligners, fontsize=13)

    ax.set_ylabel('Number of SVs', fontsize=13)
    ax.set_ylim(500, 7500)
    ax.set_yticks(np.linspace(500, 7500, 5))
    ax.set_yticklabels([f'{int(val)}' for val in np.linspace(500, 7500, 5)], fontsize=12)
    ax.legend(title='')

    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig2.tight_layout()
    fig2.savefig(f'{Figure6}/aligner_unique_counts.pdf')
    plt.show()


def aligner_unique_svtypes(workdir, datasets, aligners):

    aligner_unique_svtypes_pcrt = []
    aligner_unique_svtypes_region = []
    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']
    svtypes = ['INS', 'DUP', 'DEL']

    aligner_pcrt = {aligner: {svtype: [] for svtype in svtypes} for aligner in aligners}
    caller_unique_svtypes = {'HiFi': [], 'ONT': []}

    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for j, caller in enumerate(CALLERS):
            unique_info_out = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-unique.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])
            aligner_unique_dict = {aligner: 0 for aligner in aligners}
            aligner_unique_region_dict = {aligner: {region: 0 for region in regions} for aligner in aligners}

            aligner_unique_svtypes_dict = {aligner: [0 for i in range(len(svtypes))] for aligner in aligners}
            svtypes_region_dict = {aligner: {region: [0 for i in range(len(svtypes))] for region in regions} for aligner in aligners}
            unique_svtypes = {aligner: {ele: 0 for ele in svtypes} for aligner in aligners}

            for idx, row in df_unique.iterrows():
                svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['CALLER'], row['REGION_TYPE']
                aligner_unique_dict[aligner] += 1
                aligner_unique_region_dict[aligner][sv_region] += 1

                if svtype in svtypes:

                    aligner_unique_svtypes_dict[aligner][svtypes.index(svtype)] += 1
                    svtypes_region_dict[aligner][sv_region][svtypes.index(svtype)] += 1
                    unique_svtypes[aligner][svtype] += 1

            for aligner, svtype_count in unique_svtypes.items():
                for svtype, count in svtype_count.items():
                    caller_unique_svtypes[plat].append((caller, aligner, svtype, count / aligner_unique_dict[aligner] * 100))

            # print(aligner_unique_svtypes_dict)
            for aligner, count in aligner_unique_dict.items():
                for m, val in enumerate(aligner_unique_svtypes_dict[aligner]):
                    aligner_unique_svtypes_pcrt.append((val / count * 100, svtypes[m], plat, caller, aligner))
                    aligner_pcrt[aligner][svtypes[m]].append(val / count * 100)
                for region, region_svtypes in svtypes_region_dict[aligner].items():
                    for k, val in enumerate(region_svtypes):
                        aligner_unique_svtypes_region.append((val, svtypes[k], region, aligner, plat, caller))

    # df_unique_pcrt_svtypes = pd.DataFrame(aligner_unique_svtypes_pcrt, columns=['pcrt', 'svtype', 'plat', 'caller', 'aligner'])
    # df_unique_svtypes_region = pd.DataFrame(aligner_unique_svtypes_region, columns=['pcrt', 'svtype', 'region', 'aligner', 'plat', 'caller'])


    # sns.boxplot(data=df_unique_svtypes_region[df_unique_svtypes_region['region'] == 'Simple Repeats'], x='aligner', y='pcrt', hue='svtype', ax=ax)

    radar_dict = {'group': aligners, 'INS':[], 'DEL': [], 'DUP': []}

    for aligner in aligners:
        for svtype in svtypes:
            this_svtype_avg = np.mean(aligner_pcrt[aligner][svtype])
            radar_dict[svtype].append(this_svtype_avg)

    df_radar = pd.DataFrame(radar_dict)

    categories = list(df_radar)[1:]
    N = len(categories)

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    plt.figure(figsize=(5, 4))
    # Initialise the spider plot
    ax = plt.subplot(111, polar=True)

    # If you want the first axis to be on top:
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    # Draw one axe per variable + add labels
    plt.xticks(angles[:-1], categories, fontsize=13)

    # Draw ylabels
    ax.set_rlabel_position(315)
    plt.yticks([10, 30, 50], ["", "", ""], color="grey", size=12)
    plt.ylim(0, 50)

    for i, aligner in enumerate(aligners):
        values = df_radar.loc[i].drop('group').values.flatten().tolist()
        values += values[:1]
        ax.plot(angles, values, linewidth=2, linestyle='solid', label=aligner, color=ALIGNERCOLOR[aligner])
        ax.fill(angles, values, ALIGNERCOLOR[aligner], alpha=0.1)


    plt.tight_layout()
    plt.savefig(f'{Figure6}/aligner_unique_svtypes_radar.pdf')

    caller_color = [TOOLCOLORS[ele] for ele in CALLERS]
    svtype_color = [SVTYPECOLORS[ele] for ele in svtypes]

    for plat, plat_uniques in caller_unique_svtypes.items():
        fig, axes = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(8, 3))
        df_caller_unique_svtypes = pd.DataFrame(plat_uniques, columns=['caller', 'aligner', 'svtype', 'pcrt'])
        for col_idx, aligner in enumerate(aligners):
            this_ax = axes[col_idx]
            this_ax.set_title(aligner, fontsize=13)
            sns.boxplot(data=df_caller_unique_svtypes, x='caller', y='pcrt', hue='svtype', palette=svtype_color, ax=this_ax)
            this_ax.set_xticks(np.arange(len(CALLERS)))
            this_ax.set_xticklabels([TOOLMAP[ele] for ele in CALLERS], fontsize=13, rotation=90)

            this_ax.legend('')
            this_ax.set_xlabel('')
            if col_idx == 0:
                this_ax.set_ylabel('% of SVs', fontsize=13)
            else:
                this_ax.set_ylabel('')

            this_ax.spines['left'].set_visible(False)
            this_ax.spines['right'].set_visible(False)
            this_ax.grid(axis='y', ls='--', color='grey', lw=2)

            this_ax.set_ylim(0, 100)
            this_ax.set_yticks(np.linspace(0, 100, 5))
            this_ax.set_yticklabels([int(ele) for ele in np.linspace(0, 100, 5)], fontsize=12)

        fig.tight_layout()
        fig.savefig(f'{Figure6}/{plat}_caller_unique_svtypes_byaligner.pdf')

    plt.show()


def aligner_unique_subset_svtypes(workdir, datasets, aligners):

    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']
    svtypes = ['INS', 'DUP', 'DEL']

    svtypes_pcrt = {aligner: {region: {svtype: [] for svtype in svtypes} for region in regions} for aligner in aligners}

    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for j, caller in enumerate(CALLERS):
            unique_info_out = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-unique.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])

            this_total = {aligner: 0 for aligner in aligners}
            svtypes_count = {aligner: {region: [0, 0, 0] for region in regions} for aligner in aligners}
            for idx, row in df_unique.iterrows():
                svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['CALLER'], row['REGION_TYPE']

                if svlen >= 100 and svlen <= 1000:
                    this_total[aligner] += 1
                    if svtype in svtypes:
                        svtype_idx = svtypes.index(svtype)
                        svtypes_count[aligner][sv_region][svtype_idx] += 1

            for aligner, region_counts in svtypes_count.items():
                for region, counts in region_counts.items():
                    for i, val in enumerate(counts):
                        svtype = svtypes[i]
                        svtypes_pcrt[aligner][region][svtype].append(val / this_total[aligner] * 100)

    fig, axes = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(6, 4))
    xticks = np.arange(len(svtypes))
    bar_width = 0.5

    for col_idx, aligner in enumerate(aligners):
        ax = axes[col_idx]
        for j, region in enumerate(regions):
            this_pcrt_avg = []
            for svtype in svtypes:
                this_pcrt_avg.append(np.mean(svtypes_pcrt[aligner][region][svtype]))

            if j == 0:
                ax.bar(xticks, this_pcrt_avg, color=REGIONCOLORS[region], width=bar_width, label=region)

            else:
                bottoms = []
                for k in range(0, j):
                    this_avg = []
                    for svtype in svtypes:
                        this_avg.append(np.mean(svtypes_pcrt[aligner][regions[k]][svtype]))
                    bottoms.append(this_avg)
                bottom_sum = [sum(x) for x in zip(*bottoms)]
                ax.bar(xticks, this_pcrt_avg, bottom=bottom_sum, width=bar_width, color=REGIONCOLORS[region], label=region)

        ax.set_xticks(xticks)
        ax.set_xticklabels(svtypes, fontsize=13)

        ax.set_title(aligner, fontsize=13)

        if col_idx == 0:
            ax.set_ylabel('Percent of unique SVs', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_ylim(0, 60)
        ax.set_yticks(np.linspace(0, 60, 4))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 60, 4)], fontsize=12)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)

        # ax.legend()

    fig.tight_layout()
    plt.show()

    fig.savefig(f'{Figure6}/aligners_unique_svtypes_byregions.pdf')
'''
def aligner_unique_subset_svtypes(workdir, datasets, aligners):

    aligner_unique_svtypes_pcrt = []
    aligner_unique_svtypes_region = []
    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']
    svtypes = ['INS', 'DUP', 'DEL']

    aligner_pcrt = {aligner: {svtype: [] for svtype in svtypes} for aligner in aligners}

    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for j, caller in enumerate(CALLERS):
            unique_info_out = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-unique.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])
            aligner_unique_dict = {aligner: 0 for aligner in aligners}
            aligner_unique_region_dict = {aligner: {region: 0 for region in regions} for aligner in aligners}

            aligner_unique_svtypes_dict = {aligner: [0 for i in range(len(svtypes))] for aligner in aligners}
            svtypes_region_dict = {aligner: {region: [0 for i in range(len(svtypes))] for region in regions} for aligner in aligners}

            for idx, row in df_unique.iterrows():
                svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['CALLER'], row['REGION_TYPE']

                if svlen >= 100 and svlen <= 1000:
                    aligner_unique_dict[aligner] += 1
                    aligner_unique_region_dict[aligner][sv_region] += 1

                    if svtype in svtypes:
                        aligner_unique_svtypes_dict[aligner][svtypes.index(svtype)] += 1
                        svtypes_region_dict[aligner][sv_region][svtypes.index(svtype)] += 1

            for aligner, count in aligner_unique_dict.items():
                for m, val in enumerate(aligner_unique_svtypes_dict[aligner]):
                    aligner_unique_svtypes_pcrt.append((val / count * 100, svtypes[m], plat, caller, aligner))
                    aligner_pcrt[aligner][svtypes[m]].append(val / count * 100)
                for region, region_svtypes in svtypes_region_dict[aligner].items():
                    for k, val in enumerate(region_svtypes):
                        aligner_unique_svtypes_region.append((val, svtypes[k], region, aligner, plat, caller))

    radar_dict = {'group': aligners, 'INS':[], 'DEL': [], 'DUP': []}

    for aligner in aligners:
        for svtype in svtypes:
            this_svtype_avg = np.mean(aligner_pcrt[aligner][svtype])
            radar_dict[svtype].append(this_svtype_avg)

    df_radar = pd.DataFrame(radar_dict)

    categories = list(df_radar)[1:]
    N = len(categories)

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    fig = plt.figure(figsize=(5, 4))
    # Initialise the spider plot
    ax = plt.subplot(111, polar=True)

    # If you want the first axis to be on top:
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    # Draw one axe per variable + add labels
    plt.xticks(angles[:-1], categories, fontsize=13)

    # Draw ylabels
    ax.set_rlabel_position(45)
    plt.yticks([10, 30, 50, 70], ["10%", "30%", "50%", "70%"], color="grey", size=12)
    plt.ylim(0, 70)

    # ------- PART 2: Add plots

    # Plot each individual = each line of the data
    # I don't make a loop, because plotting more than 3 groups makes the chart unreadable

    # Ind1
    for i, aligner in enumerate(aligners):
        values = df_radar.loc[i].drop('group').values.flatten().tolist()
        values += values[:1]
        ax.plot(angles, values, linewidth=2, linestyle='solid', label=aligner, color=ALIGNERCOLOR[aligner])
        ax.fill(angles, values, ALIGNERCOLOR[aligner], alpha=0.1)


    # Add legend
    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
    plt.tight_layout()

    plt.savefig(f'{Figure6}/aligner_unique_subset_svtypes_radar.pdf')
    plt.show()

'''

