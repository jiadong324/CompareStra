#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/1/31
'''
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import vcf
import os, math

import upsetplot
from matplotlib.lines import Line2D
import pysam
from statistics import mean
from matplotlib.patches import Patch


from helpers.Constant import *
from helpers.Reader import *
from helpers.Annot import *
# from helpers.Repro import *


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


def get_aligner_compare_info(aligners, merged_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds):
    matched_list = []
    unique_list = []

    repeats = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    svs_by_regions = {rep: [0 for i in range(len(aligners))] for rep in repeats}

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')

            info_tokens = entries[7].split(';')
            info_dict = {}

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            merged_type = info_dict['SVTYPE']
            merged_id = entries[2]
            supp_vec = info_dict['SUPP_VEC']

            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            if supp > 1:
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)

                region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_aligner = aligners[supp_vec.index('1')]
                region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)
                unique_list.append((merged_id, merged_type, unique_aligner, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))

                if supp_vec == '01':
                    svs_by_regions[region_label][2] += 1
                if supp_vec == '10':
                    svs_by_regions[region_label][0] += 1

    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN',
                                       'START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(match_info_out, header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', '#CHROM', 'POS', 'END', 'SVLEN',
                                       'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(unique_info_out, sep='\t', header=True, index=False)


def compare_between_aligner(workdir, aligners, datasets):

    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for caller in CALLERS:
            outdir = f'{workdir}/{plat}/{dataset}_aligner_repro'
            if not os.path.exists(outdir):
                os.mkdir(outdir)

            tmp_file = f'{outdir}/{caller}.{dataset}.txt'
            tmp_file_writer = open(tmp_file, 'w')

            for aligner in aligners:
                vcf_file = f'{workdir}/{plat}/{aligner}_{dataset}/filtered/HG002.{caller}.filtered.vcf'
                print(f'{vcf_file}', file=tmp_file_writer)
            tmp_file_writer.close()

            merged_out_vcf = f'{outdir}/{caller}.{dataset}.jasmine.merged.vcf'
            print(f'Producing {caller} {dataset} merged calls ...')
            cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'
            os.system(cmd)

            match_info_out = f'{outdir}/{caller}.aligner-concordant.tsv'
            unique_info_out = f'{outdir}/{caller}.aligner-unique.tsv'

            print(f'Annotating {caller} {dataset} ...')
            get_aligner_compare_info(aligners, merged_out_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds)

            os.remove(tmp_file)


def print_aligner_unionset_stats(workdir, datasets, aligners):

    caller_counts = {caller: [] for caller in CALLERS}
    aligner_uniques = {aligner: [] for aligner in aligners}
    pcrt_list = []
    totals =  {'HiFi': [], 'ONT': []}
    # totals_of_callers = {caller: [] for aligner in aligners}

    totals_of_callers = []

    for platform_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for caller in CALLERS:
            merged_vcf = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.{dataset}.jasmine.merged.vcf'
            suppd_dict, merged_total = get_survivor_supp(merged_vcf)
            totals_of_callers.append(merged_total)
            sums = 0
            for supp, count in suppd_dict.items():
                if supp > 1:
                    sums += count
                    # pcrt_list.append(count / merged_total * 100)
            caller_counts[caller].append(sums)

                # totals[plat].append(merged_total)
    # print(len(pcrt_list))
    # print(np.median(pcrt_list))
    # print(np.median(totals_of_callers))

    counts_median = []
    for caller, counts in caller_counts.items():
        counts_median.append(np.median(counts))
        print(f'{caller} median: {np.median(counts)}')
    print(np.median(counts_median))


def aligner_repro(workdir, aligners, datasets):
    caller_counts = {caller: [] for caller in CALLERS}
    pcrt_list = []

    labels = [1, 2, 3, 4]

    total_pcrt = {'HiFi': {caller: [0 for i in range(len(labels))] for caller in CALLERS},
                  'ONT': {caller: [0 for i in range(len(labels))] for caller in CALLERS}}

    total_counts = []

    stacked_total = {'HiFi': [[0 for j in range(len(CALLERS))] for i in range(len(labels))],
                     'ONT': [[0 for j in range(len(CALLERS))] for i in range(len(labels))]}

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
            # upset_mem = []
            # data = []
            # for suppvec, count in suppvec_dict.items():
            #     this_mem = []
            #     for i, val in enumerate(suppvec):
            #         if val == '1':
            #             this_mem.append(aligners[i])
            #     upset_mem.append(this_mem)
            #     data.append(count)
            #
            # upset_data = upsetplot.from_memberships(upset_mem, data=data)
            # upsetplot.UpSet(upset_data, show_percentages=True)
            # plt.show()
            sums = 0
            for supp, count in supp_dict.items():
                # pcrt_list.append((count / merged_total * 100, plat, caller, f'SUPP={supp}'))
                # total_count_of_supp.append((f'SUPP={supp}', count, plat, caller))
                if supp == 1:
                    plat_supp_pcrt1[plat][this_dataset_index] += count / merged_total * 100
                    total_count_of_supp.append((f'SUPP={supp}', count, plat, caller))
                elif supp == 4:
                    plat_supp_pcrt2[plat][this_dataset_index] += count / merged_total * 100
                    total_count_of_supp.append((f'SUPP={supp}', count, plat, caller))


        ax.scatter(plat_supp_pcrt1['HiFi'], plat_supp_pcrt1['ONT'], color='#416e9f', edgecolor='black', marker=TOOLMARKERS[caller], s=100)
        ax.scatter(plat_supp_pcrt2['HiFi'], plat_supp_pcrt2['ONT'], color='#599261', edgecolor='black', marker=TOOLMARKERS[caller], s=100)

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
    # barwidth = 0.3
    # r1 = np.arange(len(CALLERS))
    # r2 = [r + barwidth + 0.05 for r in r1]
    # xticks = [r + (barwidth + 0.05)/ 2 for r in r1]
    #
    # plat_ticks = {'HiFi':r1, 'ONT':r2}
    # for plat, count_list in stacked_total.items():
    #     for i, count in enumerate(count_list):
    #
    #         if i == 0:
    #             ax.bar(plat_ticks[plat], count, hatch=PLATHATCH[plat], label=f'SUPP={labels[i]}')
    #         else:
    #             bottoms = []
    #             for j in range(0, i):
    #                 bottoms.append(count_list[j])
    #             bottom_sum = [sum(x) for x in zip(*bottoms)]
    #             ax.bar(plat_ticks[plat], count, bottom=bottom_sum, hatch=PLATHATCH[plat], label=f'SUPP={labels[i]}')


    # df_totals_supp = pd.DataFrame(total_count_of_supp, columns=['supp', 'count', 'plat', 'caller'])
    # fig1, axes = plt.subplots(1, 2, figsize=(7, 4))
    # for fig_idx, supp in enumerate(['SUPP=1', 'SUPP=4']):
    #     ax = axes[fig_idx]
    #     sns.barplot(data=df_totals_supp[df_totals_supp['supp'] == supp], x='caller', y='count', hue='plat', ax=ax)
    #     ax.set_xticklabels([TOOLMAP[caller] for caller in CALLERS], fontsize=13, rotation=90)
    #
    #     ax.set_ylabel('Number of SVs (x$10^3$)', fontsize=13)
    #     # ax.set_ylim(0, 18000)
    #     # ax.set_yticks(np.linspace(0, 18000, 4))
    #     # ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 18, 4)], fontsize=12)
    #     ax.legend(title='')
    #     ax.set_xlabel('')
    #     # ax.legend(handles=legends)
    #
    # fig1.tight_layout()
    # fig1.savefig(f'{Figure6}/aligner_repro_counts.pdf')

    plt.show()

def aligner_concordant_features(workdir, datasets, aligners):

    shift_labels = ['0,10', '10,50', '>50']

    left_bpshift_pcrt = {'HiFi': {shift: [] for shift in shift_labels},
                         'ONT': {shift: [] for shift in shift_labels}}

    right_bpshift_pcrt = {'HiFi': {shift: [] for shift in shift_labels},
                         'ONT': {shift: [] for shift in shift_labels}}

    left_shift_callers = {caller: {shift: [] for shift in shift_labels} for caller in CALLERS}
    right_shift_callers = {caller: {shift: [] for shift in shift_labels} for caller in CALLERS}

    svsize_list = []

    svsize_plat = {'HiFi': 0, 'ONT': 0}
    leftbp_pcrt_list = []
    rightbp_pcrt_list = []

    for platform_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'
        for caller in CALLERS:
            left_shift_label = {shift: 0 for shift in shift_labels}
            right_shift_label = {shift: 0 for shift in shift_labels}

            match_info_out = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-concordant.tsv'
            df_matched = pd.read_csv(match_info_out, sep='\t', header=[0])
            this_total = 0
            for idx, row in df_matched.iterrows():
                start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                if int(row['SUPP']) == 4:
                    this_total += 1
                    svsize_list.append((svlen + 1, f'{plat}-Unique'))
                    if start_std <= 10:
                        left_shift_label['0,10'] += 1
                    elif start_std <= 50:
                        left_shift_label['10,50'] += 1
                    else:
                        left_shift_label['>50'] += 1

                    # if end_std <= 10:
                    #     right_shift_label['0,10'] += 1
                    # elif end_std <= 50:
                    #     right_shift_label['10,50'] += 1
                    # elif end_std <= 100:
                    #     right_shift_label['50,100'] += 1
                    # else:
                    #     right_shift_label['>100'] += 1

                    if svlen > 10000:
                        svsize_plat[plat] += 1

            for shift in shift_labels:
                left_bpshift_pcrt[plat][shift].append(left_shift_label[shift] / this_total * 100)
                right_bpshift_pcrt[plat][shift].append(right_shift_label[shift] / this_total * 100)

                leftbp_pcrt_list.append((plat, TOOLMAP[caller], left_shift_label[shift] / this_total * 100, shift))
                # rightbp_pcrt_list.append((plat, caller, left_shift_label[shift] / this_total * 100, shift))

                # left_shift_callers[caller][shift].append(left_shift_label[shift] / len(df_matched) * 100)
                # right_shift_callers[caller][shift].append(right_shift_label[shift] / len(df_matched) * 100)

    print(svsize_plat)
    # for plat, shift_dict in left_bpshift_pcrt.items():
    #     print(f'{plat} =======')
    #     for shift, pcrt_list in shift_dict.items():
    #         print(f'{shift} mean: {np.mean(pcrt_list)}')
    #
    # for plat, shift_dict in right_bpshift_pcrt.items():
    #     print(f'{plat} =======')
    #     for shift, pcrt_list in shift_dict.items():
    #         print(f'{shift} mean: {np.mean(pcrt_list)}')
    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(7, 4))
    df_leftbp = pd.DataFrame(leftbp_pcrt_list, columns=['plat', 'caller', 'pcrt', 'shift'])
    plat_order = ['HiFi', 'ONT']
    plat_color = [PLATCOLORS[plat] for plat in plat_order]
    plat_legends = [Patch(facecolor=PLATCOLORS[plat], label=plat) for plat in plat_order]
    sns.barplot(data=df_leftbp, x='shift', y='pcrt', hue='plat', hue_order=plat_order, palette=plat_color, ax=axes[0])

    axes[0].set_xticks(np.arange(len(shift_labels)))
    axes[0].set_xticklabels(['0~10', '10~50', '>50'], fontsize=13)
    axes[0].legend(handles=plat_legends)

    sns.barplot(data=df_leftbp, x='caller', y='pcrt', hue='shift', ax=axes[1])
    axes[1].set_xticklabels(axes[1].get_xticklabels(), fontsize=13)
    axes[1].legend(title='')
    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('Percent of concordant SVs', fontsize=13)
            ax.set_ylim(0, 100)
            ax.set_yticks(np.linspace(0, 100, 5))
            ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)

        else:
            ax.set_ylabel('')

        ax.set_xlabel('')

    # for caller, shift_dict in left_shift_callers.items():
    #     print(f'{caller} ======')
    #     for shift, pcrt_list in shift_dict.items():
    #         print(f'{shift} mean: {np.mean(pcrt_list)}')
    #
    # for caller, shift_dict in right_shift_callers.items():
    #     print(f'{caller} ======')
    #     for shift, pcrt_list in shift_dict.items():
    #         print(f'{shift} mean: {np.mean(pcrt_list)}')

    # df_svsize = pd.DataFrame(svsize_list, columns=['svlen', 'plat'])
    # fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    # sns.histplot(data=df_svsize, x='svlen', log_scale=True, hue='plat', bins=50, ax=ax)


    plt.tight_layout()
    plt.show()

def plot_aligner_unique_features(workdir, datasets, aligners):
    regions = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    uniques_count = {caller: [0 for i in range(len(datasets))] for caller in CALLERS}
    uniques_count1 = {aligner: [0 for i in range(len(datasets))] for aligner in aligners}

    unique_svsize = []
    unique_svsize1 = {aligner: 0 for aligner in aligners}

    uniques_svtypes = {aligner: {} for aligner in aligners}

    uniques_count_byaligners = {aligner: 0 for aligner in aligners}

    unique_byregions = {'HiFi': {caller: [0 for i in range(len(regions))] for caller in CALLERS},
                        'ONT': {caller: [0 for i in range(len(regions))] for caller in CALLERS}}


    unique_svsize_plat = {'HiFi': 0, 'ONT': 0}

    for dataset_idx, dataset in enumerate(datasets):
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        for j, caller in enumerate(CALLERS):

            unique_info_out = f'{workdir}/{plat}/{dataset}_aligner_repro/{caller}.aligner-unique.tsv'
            df_unique = pd.read_csv(unique_info_out, sep='\t', header=[0])

            for idx, row in df_unique.iterrows():
                svtype, svlen, aligner, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['CALLER'], row['REGION_TYPE']
                region_idx = regions.index(sv_region)
                uniques_count[caller][dataset_idx] += 1
                uniques_count1[aligner][dataset_idx] += 1
                unique_byregions[plat][caller][region_idx] += 1

                if svlen < 1000 and svlen > 100:
                    unique_svsize1[aligner] += 1
                    unique_svsize_plat[plat] += 1
                    uniques_count_byaligners[aligner] += 1
                    if svtype in uniques_svtypes[aligner]:
                        uniques_svtypes[aligner][svtype] += 1
                    else:
                        uniques_svtypes[aligner][svtype] = 1

                if svtype != 'BND':
                    unique_svsize.append((svlen + 1, f'{plat}-Unique', aligner))



    print(unique_svsize1)

    for key, vals in uniques_count1.items():
        print(f'{key} median: {np.median(vals)}')
    #
    # for plat, region_counts in unique_byregions.items():
    #     print(f'{plat} =====')
    #     for caller, counts in region_counts.items():
    #         print(f'{caller} {counts}')

    # svtypes_list = []
    for aligner, svtypes_dict in uniques_svtypes.items():
        aligner_total = uniques_count_byaligners[aligner]
        print(f'{aligner} ====== ')
        for svtype, count in svtypes_dict.items():
            print(f'{svtype}: {count / aligner_total * 100}')
    #         svtypes_list.append((svtype, count / aligner_total * 100, aligner))
    # df_svtypes = pd.DataFrame(svtypes_list, columns=['svtype', 'pcrt', 'aligner'])
    # fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    # sns.lineplot(data=df_svtypes, x='svtype', y='pcrt', hue='aligner', ax=ax)

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(7, 3))
    df_svsize = pd.DataFrame(unique_svsize, columns=['svlen', 'plat', 'aligner'])
    sns.histplot(data=df_svsize, x='svlen', log_scale=True, hue='aligner', bins=50, ax=axes[0])
    sns.histplot(data=df_svsize, x='svlen', log_scale=True, hue='plat', bins=50, ax=axes[1])

    for ax in axes:
        ax.set_xlabel('')
        # ax.legend(title='')

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{iMACALIGNFIGURE}/aligner_unique_svsize.pdf')

'''

def aligner_unique_loci_by_regions(workdir, platforms, aligners, figdir):

    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']
    xticks = np.arange(len(region_labels))

    fig, axes = plt.subplots(len(CALLERS), len(platforms), sharex='col', sharey='row', figsize=(14, 8))

    for row_idx, caller in enumerate(CALLERS):
        for col_idx, platform in enumerate(platforms):
            unique_loci_rptypes = {aligner: [0 for i in range(len(region_labels))] for aligner in aligners}
            this_ax = axes[row_idx][col_idx]

            if row_idx == 0:
                this_ax.set_title(PLATMAP[platform], fontsize=13)

            if col_idx == 0:
                this_ax.set_ylabel(TOOLMAP[caller])

            for aligner_idx, aligner in enumerate(aligners):
                aligner_unique_out = f'{workdir}/align_repro/{platform}/{caller}.{platform}.{aligner}.unique.tsv'
                df_aligner_unique = pd.read_csv(aligner_unique_out, sep='\t',names=['chrom', 'start', 'end', 'svlen', 'svtype', 'rptype', 'rppcrt', 'aligner', 'platform'])

                for idx, row in df_aligner_unique.iterrows():
                    rptype = row['rptype']
                    region = 'Unique'
                    if rptype == 'VNTR' or rptype == 'STR':
                        region = 'Simple Repeats'
                    elif rptype == 'SegDup':
                        region = 'Segment Dup'
                    elif rptype != 'None':
                        region = 'Repeat Masked'

                    unique_loci_rptypes[aligner][region_labels.index(region)] += 1

            for i, aligner in enumerate(aligners):
                this_aligner_counts = unique_loci_rptypes[aligner]
                if i == 0:
                    this_ax.bar(xticks, this_aligner_counts, label=aligner, color=ALIGNERCOLOR[aligner], edgecolor='w', width=0.4)
                else:
                    bottoms = []
                    for j in range(0, i):
                        bottoms.append(unique_loci_rptypes[aligners[j]])
                    bottom_sum = [sum(x) for x in zip(*bottoms)]
                    this_ax.bar(xticks, this_aligner_counts, bottom=bottom_sum, label=aligner, color=ALIGNERCOLOR[aligner], edgecolor='w', width=0.4)

            this_ax.set_xticks(xticks)
            this_ax.set_xticklabels(region_labels, rotation=90, fontsize=12)
            # this_ax.legend()

    plt.tight_layout()
    plt.show()

    fig.savefig(f'{figdir}/aligner_unique_loci_byregions.pdf')


def aligner_unique_loci_size(workdir, platforms, aligners, figdir):
    sns.set_theme(style="ticks", rc=custom_params, font="Arial", font_scale=1.0)
    fig, axes = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10, 6))

    aligner_legends = [Patch(facecolor=ALIGNERCOLOR[aligner], edgecolor='black', label=aligner) for aligner in aligners]

    hue_color = [ALIGNERCOLOR[aligner] for aligner in aligners]
    row_idx = 0
    col_idx = 0

    for caller_idx, caller in enumerate(CALLERS):

        if caller_idx == 2:
            row_idx = 1
            col_idx = 0

        sv_size_list = []

        for platform_idx, platform in enumerate(platforms):
            for aligner_idx, aligner in enumerate(aligners):
                aligner_unique_out = f'{workdir}/align_repro/{platform}/{caller}.{platform}.{aligner}.unique.tsv'

                df_aligner_unique = pd.read_csv(aligner_unique_out, sep='\t', names=['chrom', 'start', 'end', 'svlen', 'svtype', 'rptype',
                                                       'rppcrt', 'aligner', 'platform'])

                for idx, row in df_aligner_unique.iterrows():
                    svlen = abs(int(row['svlen']))
                    if svlen < 100000:
                        sv_size_list.append((svlen, aligner))

        df_sv_size = pd.DataFrame(sv_size_list, columns=['svlen', 'aligner'])

        this_ax = axes[row_idx][col_idx]

        this_ax.set_xlabel('')
        this_ax.set_title(TOOLMAP[caller], fontsize=13)

        sns.histplot(data=df_sv_size, x='svlen', hue='aligner', hue_order=aligners, palette=hue_color, log_scale=True, ax=this_ax)
        this_ax.legend(handles=aligner_legends, title='')

        col_idx += 1

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{figdir}/aligner_unique_size.pdf')


def aligner_unique_loci_svtypes(workdir, platforms, aligners, figdir):

    svtypes_order = ['DEL', 'INS', 'DUP']

    sns.set_theme(style="ticks", rc=custom_params, font="Arial", font_scale=1.0)

    fig, axes = plt.subplots(3, len(platforms), sharex='col', sharey='row', figsize=(15, 8))

    xticks = np.arange(len(aligners))

    for col_idx, platform in enumerate(platforms):

        svtypes_dict = {svtype: {caller: [0 for i in range(len(aligners))] for caller in CALLERS} for svtype in svtypes_order}

        for caller_idx, caller in enumerate(CALLERS):
            for aligner_idx, aligner in enumerate(aligners):


                aligner_unique_out = f'{workdir}/align_repro/{platform}/{caller}.{platform}.{aligner}.unique.tsv'
                df_aligner_unique = pd.read_csv(aligner_unique_out, sep='\t',names=['chrom', 'start', 'end', 'svlen', 'svtype', 'rptype', 'rppcrt','aligner', 'platform'])

                for idx, row in df_aligner_unique.iterrows():
                    svtype = row['svtype']
                    if svtype in svtypes_dict:
                        svtypes_dict[svtype][caller][aligner_idx] += 1

        for row_idx, svtype in enumerate(svtypes_order):
            this_ax = axes[row_idx][col_idx]
            if row_idx == 0:
                this_ax.set_title(PLATMAP[platform], fontsize=13)

            if col_idx == 0:
                this_ax.set_ylabel(f'Aligner unique {svtype} loci', fontsize=13)

            for caller, counts in svtypes_dict[svtype].items():
                this_ax.plot(xticks, counts, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], markersize=9, lw=2, label=TOOLMAP[caller])

            this_ax.legend()

            if row_idx == 0:
                this_ax.set_ylim([0, 1500])
                this_ax.set_yticks(np.linspace(0, 1500, 6))
                this_ax.set_yticklabels([int(val) for val in np.linspace(0, 1500, 6)], fontsize=12)

            elif row_idx == 1:
                this_ax.set_ylim([0, 2000])
                this_ax.set_yticks(np.linspace(0, 2000, 6))
                this_ax.set_yticklabels([int(val) for val in np.linspace(0, 2000, 6)], fontsize=12)

            else:
                this_ax.set_ylim([0, 500])
                this_ax.set_yticks(np.linspace(0, 500, 6))
                this_ax.set_yticklabels([int(val) for val in np.linspace(0, 500, 6)], fontsize=12)


            this_ax.set_xticks(xticks)
            this_ax.set_xticklabels(aligners, fontsize=12, rotation=90)

            # ax.set_yticks(np.linspace(0, 1, 5))
            # ax.set_yticklabels([f'{int(val*100)}%' for val in np.linspace(0, 1, 5)], fontsize=12)

            this_ax.tick_params(width=2)
            for axis in ['bottom', 'left']:
                this_ax.spines[axis].set_linewidth(2)


    plt.tight_layout()
    plt.show()

    fig.savefig(f'{figdir}/aligner_unique_svtypes.pdf')

def aligner_unique_loci_rptypes(workdir, platforms, aligners, figdir):


    rptypes = ['VNTR', 'STR', 'LTR', 'LINE', 'SINE', 'SegDup', 'None']
    # rptype_legends = [Line2D([0], [0], color=color, label=rptype, lw=2) for rptype, color in REPCOLOR.items()]

    sns.set_theme(style="ticks", rc=custom_params, font="Arial", font_scale=1.0)
    fig, axes = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10, 6))

    row_idx = 0
    col_idx = 0
    for caller_idx, caller in enumerate(CALLERS):

        if caller_idx == 2:
            row_idx = 1
            col_idx = 0

        rptypes = {rep: [[0, 0, 0, 0] for i in platforms] for rep in rptypes}

        total_svs = [[0, 0, 0, 0] for i in platforms]

        for platform_idx, platform in enumerate(platforms):
            for aligner_idx, aligner in enumerate(aligners):

                aligner_unique_out = f'{workdir}/align_repro/{platform}/{caller}.{platform}.{aligner}.unique.tsv'

                df_aligner_unique = pd.read_csv(aligner_unique_out, sep='\t', names=['chrom', 'start', 'end', 'svlen', 'svtype', 'rptype', 'rppcrt', 'aligner', 'platform'])

                for idx, row in df_aligner_unique.iterrows():
                    if row['rptype'] in rptypes:
                        total_svs[platform_idx][aligner_idx] += 1
                        rptypes[row['rptype']][platform_idx][aligner_idx] += 1

        pcrt_list = []
        for j, rptype in enumerate(rptypes):

            for k, ele in enumerate(rptypes[rptype]):
                for n in range(len(aligners)):
                    pcrt_list.append((ele[n] * 100 / total_svs[k][n], aligners[n], rptype))

                # pcrt_list.append((ele[1] * 100 / total_svs[k][1], 'ngmlr', rptype))

        df_pcrt = pd.DataFrame(pcrt_list, columns=['pcrt', 'aligner', 'rptype'])
        this_ax = axes[row_idx][col_idx]

        sns.boxplot(data=df_pcrt, x='rptype', y='pcrt', hue='aligner', hue_order=aligners, palette=[ALIGNERCOLOR[val] for val in aligners], ax=this_ax)

        this_ax.set_title(TOOLMAP[caller], fontsize=13)
        this_ax.legend(title='')

        this_ax.set_xlabel('')

        this_ax.set_ylim(0, 80)
        this_ax.set_yticks(np.linspace(0, 80, 5))
        this_ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 80, 5)], fontsize=12)
        if col_idx == 0:
            this_ax.set_ylabel('% of SVs', fontsize=13)
        else:
            this_ax.set_ylabel('')

        this_ax.tick_params(width=2)
        for axis in ['bottom', 'left']:
            this_ax.spines[axis].set_linewidth(2)

        col_idx += 1

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{figdir}/aligner_unique_rptypes.pdf')

def aligner_concordant_bpshift(workdir, platforms, figdir):

    shift_labels = ['0,10', '10,50', '50,100', '>100']
    sns.set_theme(style="ticks", rc=custom_params, font="Arial", font_scale=1.0)

    fig, axes = plt.subplots(2, len(CALLERS), sharex='col', sharey='row', figsize=(12, 6))
    fig1, axes1 = plt.subplots(2, len(CALLERS), sharex='col', sharey='row', figsize=(12, 6))

    xticks = np.arange(len(platforms))

    svtypes = ['ins', 'del']

    for row_idx, svtype in enumerate(svtypes):
        for col_idx, caller in enumerate(CALLERS):
            minshift_ax = axes[row_idx][col_idx]
            maxshift_ax = axes1[row_idx][col_idx]
            for ax in [minshift_ax, maxshift_ax]:
                if col_idx == 0:
                    ax.set_ylabel(f'Number of {SVTYPEMAP[svtype]} (x1000)', fontsize=13)
                if row_idx == 0:
                    ax.set_title(f'{TOOLMAP[caller]}', fontsize=13)

                ax.set_xticks(xticks)
                ax.set_xticklabels([PLATMAP[plat] for plat in platforms], fontsize=12, rotation=90)
                ax.tick_params(width=2)
                for axis in ['bottom', 'left']:
                    ax.spines[axis].set_linewidth(2)

                ax.set_ylim(0, 10000)
                ax.set_yticks(np.linspace(0, 10000, 6))
                ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 10, 6)], fontsize=12)


            minshift_by_cutoff = {shift: [0 for k in range(len(platforms))] for shift in shift_labels}
            maxshift_by_cutoff = {shift: [0 for k in range(len(platforms))] for shift in shift_labels}

            bp_rptype = {}

            platform_total = []
            acc_total = 0
            inacc_total = 0

            for j, platform in enumerate(platforms):

                overlaps_file = f'{workdir}/align_repro/{platform}/{caller}.aligner-concordant.{svtype}.tsv'
                df_overlaps = pd.read_csv(overlaps_file, sep='\t', header=[0])
                platform_total.append(len(df_overlaps))

                for idx, row in df_overlaps.iterrows():
                    min_bpshift = int(row['MIN_BPDIFF'])
                    if min_bpshift <= 10:
                        minshift_by_cutoff['0,10'][j] += 1
                    elif min_bpshift <= 50:
                        minshift_by_cutoff['10,50'][j] += 1
                    elif min_bpshift <= 100:
                        minshift_by_cutoff['50,100'][j] += 1
                    else:
                        minshift_by_cutoff['>100'][j] += 1

                    max_bpshift = int(row['MAX_BPDIFF'])
                    if max_bpshift <= 10:
                        maxshift_by_cutoff['0,10'][j] += 1
                    elif max_bpshift <= 50:
                        maxshift_by_cutoff['10,50'][j] += 1
                    elif max_bpshift <= 100:
                        maxshift_by_cutoff['50,100'][j] += 1
                    else:
                        maxshift_by_cutoff['>100'][j] += 1

                    rptype = row['RPTYPE']
                    if rptype == 'None':
                        continue

                    if min_bpshift <= 50:
                        acc_total += 1
                        if rptype in bp_rptype:
                            bp_rptype[rptype][0] += 1
                        else:
                            bp_rptype[rptype] = [1, 0]
                    else:
                        inacc_total += 1
                        if rptype in bp_rptype:
                            bp_rptype[rptype][1] += 1
                        else:
                            bp_rptype[rptype] = [0, 1]


            for l, shift in enumerate(shift_labels):
                this_minshift = minshift_by_cutoff[shift]
                if l == 0:
                    minshift_ax.bar(xticks, this_minshift, color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)
                else:
                    bottom = []
                    for m in range(0, l):
                        bottom.append([val for val in minshift_by_cutoff[shift_labels[m]]])
                    bottom_sum = [sum(x) for x in zip(*bottom)]
                    minshift_ax.bar(xticks, this_minshift, bottom=bottom_sum, color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)

                this_maxshift = maxshift_by_cutoff[shift]
                if l == 0:
                    maxshift_ax.bar(xticks, this_maxshift, color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w',
                                    label=shift)
                else:
                    bottom = []
                    for m in range(0, l):
                        bottom.append([val for val in maxshift_by_cutoff[shift_labels[m]]])
                    bottom_sum = [sum(x) for x in zip(*bottom)]
                    maxshift_ax.bar(xticks, this_maxshift, bottom=bottom_sum, color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)

            maxshift_ax.legend(ncol=2)
            minshift_ax.legend(ncol=2)

            acc_pcrt = []
            inacc_pcrt = []
            rptype_labels = []
            for rptype, counts in bp_rptype.items():
                rptype_labels.append(rptype)
                acc_pcrt.append(counts[0] / acc_total)
                inacc_pcrt.append(counts[1] / inacc_total)

            print(acc_pcrt)

    fig.tight_layout()
    fig1.tight_layout()
    plt.show()
    fig.savefig(f'{figdir}/aligner_concordants_minbpshift.pdf')
    fig1.savefig(f'{figdir}/aligner_concordants_maxbpshift.pdf')


def aligner_unique_rptypes(platforms):
    aligners = ['minimap2', 'ngmlr']

    rptype_legends = [Line2D([0], [0], color=color, label=rptype, lw=2) for rptype, color in REPCOLOR.items()]

    fig, axes = plt.subplots(2, len(CALLERS), sharex='col', sharey='row', figsize=(18, 8))
    xticks = np.arange(len(platforms))

    for caller_idx, caller in enumerate(CALLERS):


        rptypes = {'VNTR': [[0, 0] for i in platforms], 'STR': [[0, 0] for i in platforms],
                   'LTR': [[0, 0] for i in platforms], 'LINE': [[0, 0] for i in platforms],
                   'SINE': [[0, 0] for i in platforms]}

        total_svs = [[0, 0] for i in platforms]

        for platform_idx, platform in enumerate(platforms):
            for aligner_idx in range(len(aligners)):
                aligner = aligners[aligner_idx]
                aligner_unique_out = f'{HG002DIR}/survivor/{platform}/{caller}.{platform}.{aligner}.unique.tsv'

                df_aligner_unique = pd.read_csv(aligner_unique_out, sep='\t',
                                                names=['chrom', 'start', 'end', 'svlen', 'svtype', 'rptype',
                                                       'rppcrt', 'aligner', 'platform'])

                for idx, row in df_aligner_unique.iterrows():
                    total_svs[platform_idx][aligner_idx] += 1
                    if row['rptype'] in rptypes:
                        rptypes[row['rptype']][platform_idx][aligner_idx] += 1


        for rptype, count_list in rptypes.items():
            minimap2_rptype_pcrt = []
            ngmlr_rptype_pcrt = []
            for j, counts in enumerate(count_list):
                minimap2_rptype_pcrt.append(counts[0] / total_svs[j][0])
                ngmlr_rptype_pcrt.append(counts[1] / total_svs[j][1])

            axes[0][caller_idx].plot(xticks, minimap2_rptype_pcrt, color=REPCOLOR[rptype], marker='o', markersize=8, lw=2, label=rptype)
            axes[1][caller_idx].plot(xticks, ngmlr_rptype_pcrt, color=REPCOLOR[rptype], marker='o', markersize=8,lw=2, label=rptype)

    for i in range(2):
        for j in range(len(CALLERS)):
            axes[i][j].set_title(TOOLMAP[CALLERS[j]], fontsize=13)
            axes[i][j].set_xticks(xticks)
            axes[i][j].set_xticklabels([PLATMAP[plat] for plat in platforms], rotation=90, fontsize=12)

            axes[i][j].set_ylim(0, 0.3, 4)
            axes[i][j].set_yticks(np.linspace(0, 0.3, 4))
            axes[i][j].set_yticklabels([f'{int(val*100)}%' for val in np.linspace(0, 0.3, 4)], fontsize=12)

            axes[i][j].legend(handles=rptype_legends, ncol=2)

            if i == 0 and j == 0:
                axes[i][j].set_ylabel('% of minimap2 calls', fontsize=13)
            if i == 1 and j == 0:
                axes[i][j].set_ylabel('% of ngmlr calls', fontsize=13)

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{FIGURE}/aligner_unique_rptypes.pdf')
'''