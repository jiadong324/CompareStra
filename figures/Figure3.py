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
plt.rcParams["xtick.minor.width"] = 1
plt.rcParams["xtick.labelsize"] = 11

plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.minor.width"] = 1
plt.rcParams["ytick.labelsize"] = 11

plt.rcParams["axes.linewidth"] = 2


def platform_repro_overview(workdir):

    aligner = 'minimap2'

    plat_assembler = {'HiFi':['flye', 'hifiasm'], 'ONT': ['flye', 'shasta']}

    labels = [1, 2, 3]
    xticks = np.arange(len(labels))

    hue_order = ['Assembly-based', 'Read-based']
    # hue_color = [STRACOLORS[ele.split('-')[0]] for ele in hue_order]

    for flag, assemblers in plat_assembler.items():
        fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(4, 3))

        for col_idx, assembler in enumerate(assemblers):
            this_ax = axes[col_idx]
            supp_pcrt = []
            for caller in ASMCALLERS:
                merged_vcf = f'{workdir}/assm_dataset_repro/{flag}/{caller}.{assembler}.jasmine.merged.vcf'
                supp_dict, merged_total = get_survivor_supp(merged_vcf)

                assm_pcrt = []
                for supp in labels:
                    count = supp_dict[supp]
                    supp_pcrt.append((count / merged_total, f'SUPP={supp}', 'Assembly-based'))
                    assm_pcrt.append(count / merged_total)

                this_ax.plot(xticks, assm_pcrt, color=TOOLCOLORS[caller], marker='o', lw=2, markersize=7, label=f'{assembler}+{TOOLMAP[caller]}')
                this_ax.set_title(flag, fontsize=13)
            # df_supp_pcrt = pd.DataFrame(supp_pcrt, columns=['pcrt', 'supp', 'strategy'])

            # bars = sns.barplot(data=df_supp_pcrt, x='supp', y='pcrt', alpha=0.6, edgecolor='black', capsize=0.15, color=STRACOLORS['Assembly'], ax=this_ax)
            #
            # for bar in bars.patches:
            #     x = bar.get_x()
            #     width = bar.get_width()
            #     center = x + width / 2.
            #     bar.set_facecolor('white')
            #     bar.set_x(center - 0.4 / 2.)
            #     bar.set_width(0.4)
            #     bar.set_lw(2)

        for i, ax in enumerate(axes):
            if i == 0:
                ax.set_ylabel(f'% of concordant ({flag})', fontsize=13)
                ax.set_ylim(0, 1)
                ax.set_yticks(np.linspace(0, 1, 5))
                ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)], fontsize=12)
            else:
                ax.set_ylabel('')

            ax.set_xlabel('')
            # ax.legend()

            ax.set_xticks(np.arange(len(labels)))
            ax.set_xticklabels(labels, fontsize=13)

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            # this_ax1.set_ylim(1, 0)
            # this_ax1.set_yticks(np.linspace(1, 0, 5))
            # this_ax1.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(1, 0, 5)], fontsize=12)
            # this_ax1.set_ylabel('Percent of caller', fontsize=13)
        fig.tight_layout()

        fig.savefig(f'{Figure3}/{flag}_assm_based_datasets_repro.pdf')


    fig1, axes = plt.subplots(1, 2, sharey='row', sharex='col', figsize=(5, 3))

    for col_idx, flag in enumerate(['HiFi', 'ONT']):
        supp_pcrt = []
        this_ax = axes[col_idx]

        for caller in CALLERS:
            merged_vcf = f'{workdir}/read_dataset_repro/{flag}/{caller}.{aligner}.jasmine.merged.vcf'
            supp_dict, merged_total = get_survivor_supp(merged_vcf)

            read_pcrt = []
            for supp in labels:
                count = 0
                if supp in supp_dict:
                    count = supp_dict[supp]
                supp_pcrt.append((count / merged_total, f'SUPP={supp}', 'Read-based'))
                read_pcrt.append(count / merged_total)

            this_ax.plot(xticks, read_pcrt, color=TOOLCOLORS[caller], marker='o', lw=2, markersize=7, label=TOOLMAP[caller])
            this_ax.set_title(flag, fontsize=13)
        # df_supp_pcrt = pd.DataFrame(supp_pcrt, columns=['pcrt', 'supp', 'strategy'])
        #
        # bars = sns.barplot(data=df_supp_pcrt, x='supp', y='pcrt', alpha=0.6, edgecolor='black',
        #                    hue_order=hue_order, capsize=0.15, color=STRACOLORS['Read'], ax=this_ax)
        #
        # for bar in bars.patches:
        #     x = bar.get_x()
        #     width = bar.get_width()
        #     center = x + width / 2.
        #     bar.set_facecolor('white')
        #     bar.set_x(center - 0.4 / 2.)
        #     bar.set_width(0.4)
        #     bar.set_lw(2)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel(f'% of concordant', fontsize=13)
            ax.set_ylim(0, 1)
            ax.set_yticks(np.linspace(0, 1, 5))
            ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)], fontsize=12)
            ax.legend()
        else:
            ax.set_ylabel('')
            ax.legend('')

        ax.set_xlabel('')

        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, fontsize=13)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig1.tight_layout()

    plt.show()
    fig1.savefig(f'{Figure3}/read_based_dataset_repro_of_platforms.pdf')


def platform_unique_region(workdir, aligner):
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']

    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']
    plat_assembler = {'hifi': ['hifiasm', 'flye'], 'ont': ['shasta', 'flye']}
    # unique_pcrt = {region: [[0 for i in range(len(callers))], [0 for i in range(len(callers))]] for region in region_labels}

    asm_unique_region = []
    read_unique_region = []
    all_unique_region = []
    read_unique_sv = []
    assm_unique_sv = []

    for k, plat in enumerate(['hifi', 'ont']):
        for asm_caller in ASMCALLERS:
            # caller_index = callers.index(asm_caller)
            for assembler in plat_assembler[plat]:
                matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{PLATMAP[plat]}/{asm_caller}.{assembler}.{PLATMAP[plat]}.unique.tsv', sep='\t',header=[0])
                region_dict = {region: 0 for region in region_labels}
                for idx, row in matched_info.iterrows():
                    svtype, svlen, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                    region_dict[sv_region] += 1
                    if svlen >= 50:
                        assm_unique_sv.append((svtype, svlen, sv_region, asm_caller, assembler, PLATMAP[plat]))

                for region, count in region_dict.items():
                    # unique_pcrt[region][k][caller_index] += count / len(matched_info) * 100

                    this_region_pcrt = count / len(matched_info) * 100
                    asm_unique_region.append((this_region_pcrt, region, assembler, asm_caller, PLATMAP[plat]))
                    all_unique_region.append((this_region_pcrt, region, TOOLMAP[asm_caller], PLATMAP[plat], 'Assembly-based'))

        for caller in CALLERS:
            # caller_index = callers.index(caller)
            matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.{plat}.unique.tsv',sep='\t', header=[0])
            region_dict = {region: 0 for region in region_labels}
            for idx, row in matched_info.iterrows():
                svtype, svlen, sv_region = row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                region_dict[sv_region] += 1
                if svlen >= 50:
                    read_unique_sv.append((svtype, svlen, sv_region, caller, PLATMAP[plat]))

            for region, count in region_dict.items():
                # unique_pcrt[region][k][caller_index] += count / len(matched_info) * 100
                this_region_pcrt = count / len(matched_info) * 100
                read_unique_region.append((this_region_pcrt, region, caller, PLATMAP[plat]))

                all_unique_region.append((this_region_pcrt, region, TOOLMAP[caller], PLATMAP[plat], 'Read-based'))

    # all_legends = [Patch(facecolor=REGIONCOLORS[region], label=region) for region in region_labels]
    # caller_legends = [Line2D([0], [0], marker='o', mfc='white', mec='black', color='w',
    #                          markersize=10, label=TOOLMAP[caller]) for caller in callers]

    # all_legends.extend(caller_legends)
    #
    # fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    # for region, [hifi_pcrt, ont_pcrt] in unique_pcrt.items():
    #     for i, caller in enumerate(callers):
    #
    #         ax.scatter(hifi_pcrt[i], ont_pcrt[i], color=REGIONCOLORS[region], edgecolor='black', marker=TOOLMARKERS[caller], s=100)
    #
    # ax.set_xlabel('HiFi dataset', fontsize=13)
    # ax.set_xlim(0, 100)
    # ax.set_xticks(np.linspace(0, 100, 5))
    # ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
    #
    # ax.set_ylabel('ONT dataset', fontsize=13)
    # ax.set_ylim(0, 100)
    # ax.set_yticks(np.linspace(0, 100, 5))
    # ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
    #
    # ax.plot([0, 100], [0, 100], color='#c96150', ls='--')
    #
    # ax.legend(handles=all_legends)
    #
    # ax1 = plt.axes([.62, .25, .25, .25])
    # for i in range(1, len(region_labels)):
    #     hifi_pcrt, ont_pcrt = unique_pcrt[region_labels[i]]
    #     for j, caller in enumerate(callers):
    #
    #         ax1.scatter(hifi_pcrt[j], ont_pcrt[j], color=REGIONCOLORS[region_labels[i]], edgecolor='black', marker=TOOLMARKERS[caller], s=100)
    #
    # ax1.set_ylim(0, 25)
    # ax1.set_yticks(np.linspace(0, 25, 2))
    # ax1.set_yticklabels(['0%', '25%'], fontsize=12)
    #
    # ax1.set_xlim(0, 25)
    # ax1.set_xticks(np.linspace(0, 25, 2))
    # ax1.set_xticklabels(['0%', '25%'], fontsize=12)
    # ax1.plot([0, 25], [0, 25], color='#c96150', ls='--')

    df_all_unique_region = pd.DataFrame(all_unique_region, columns=['pcrt', 'region', 'caller', 'plat', 'stra'])
    fig, ax = plt.subplots(1, 1, sharex='col', sharey='row', figsize=(4, 3))

    stra_order = ['Assembly-based', 'Read-based']
    stra_colors = [STRACOLORS[ele.split('-')[0]] for ele in stra_order]
    sns.barplot(data=df_all_unique_region[df_all_unique_region['plat'] == 'HiFi'], y='region', x='pcrt', hue='stra',
                hue_order=stra_order, palette=stra_colors, capsize=.2, ci='sd', ax=ax)

    # axes[0].set_ylabel('% of SVs (HiFi)', fontsize=13)

    # hue_order = [TOOLMAP[caller] for caller in callers]
    # hue_colors = [TOOLCOLORS[caller] for caller in callers]
    # sns.stripplot(data=df_all_unique_region, x='region', y='pcrt', hue='caller', hue_order=hue_order, palette=hue_colors, size=7, ax=axes[1])

    # sns.barplot(data=df_all_unique_region[df_all_unique_region['plat'] == 'ONT'], x='region', y='pcrt', hue='stra',
    #             hue_order=stra_order, palette=stra_colors, capsize=.2, ci='sd', ax=axes[1])
    # axes[1].set_ylabel('% of SVs (ONT)', fontsize=13)

    # for i, ax in enumerate(axes):
    ax.set_xlabel('% of SVs', fontsize=13)
    ax.set_xlim(0, 100)
    ax.set_xticks(np.linspace(0, 100, 5))
    ax.set_xticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)

    ax.set_yticks(np.arange(len(region_labels)))
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)
    ax.set_ylabel('')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=2)

    ax.legend(title='')

    fig.tight_layout()

    df_assm_unique_svs = pd.DataFrame(assm_unique_sv, columns=['svtype', 'svlen', 'region', 'caller', 'assembler', 'plat'])
    df_read_unique_svs = pd.DataFrame(read_unique_sv, columns=['svtype', 'svlen', 'region', 'caller', 'plat'])

    region_order = ['Unique', 'Segment Dup', 'Repeat Masked', 'Simple Repeats']

    region_color = [REGIONCOLORS[ele] for ele in region_order]

    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(7, 3))
    sns.histplot(data=df_assm_unique_svs[df_assm_unique_svs['svlen'] < 1000], x='svlen', hue='region', bins=40, log_scale=[False, True],
                 hue_order=region_order, palette=region_color, kde=True, ax=axes[0])

    sns.histplot(data=df_read_unique_svs[df_read_unique_svs['svlen'] < 1000], x='svlen', hue='region', bins=40, log_scale=[False, True],
                 hue_order=region_order, palette=region_color, kde=True, ax=axes[1])


    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('# of SVs (x1000)', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_xlabel('')
        # ax.legend()

        ax.set_ylim(1, 20000)
        # ax.set_yticks(np.linspace(0, 20000, 5))
        # ax.set_yticklabels([int(val) for val in np.linspace(0, 20, 5)], fontsize=12)


        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig1.tight_layout()
    plt.show()
    fig.savefig(f'{Figure3}/platform_unique_regions.pdf')
    fig1.savefig(f'{Figure3}/platform_unique_svlen_byregions.pdf')


def platform_repro_bpstd(workdir, aligner):

    shift_labels = ['0', '0,10', '10,50', '>50']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']
    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    acc_svs_pcrt = []
    all_svs_bpstd = []
    barwidth = 0.3
    r1 = np.arange(len(callers))
    r2 = [r + barwidth + 0.05 for r in r1]
    # plat_xticks = [r1, r2]

    legends = [Patch(facecolor=BPSHIFTCOLORS[i], edgecolor='black', label=shift) for i, shift in enumerate(shift_labels)]
    legends.append(Patch(facecolor='white', edgecolor='black', hatch='///', label='ONT'))
    legends.append(Patch(facecolor='white', edgecolor='black', label='HiFi'))

    # fig, leftbp_ax = plt.subplots(1, 1, figsize=(5, 3))
    # fig1, rightbp_ax = plt.subplots(1, 1, figsize=(5, 3))


    plat_assembler = {'hifi':['hifiasm', 'flye'], 'ont':['shasta', 'flye']}

    for k, plat in enumerate(['hifi', 'ont']):

        # hatch = ''
        # if plat == 'ont':
        #     hatch = '///'

        start_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
        end_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
        caller_total = [1 for k in range(len(callers))]

        region_label_dict = {region: [] for region in region_labels}

        for asm_caller in ASMCALLERS:
            caller_index = callers.index(asm_caller)
            for assembler in plat_assembler[plat]:

                matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{PLATMAP[plat]}/{asm_caller}.{assembler}.{PLATMAP[plat]}-concordant.info.tsv', sep='\t',header=[0])
                for idx, row in matched_info.iterrows():
                    start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                    region_label_dict[sv_region].append((math.log10(start_std + 1), asm_caller, 'start_std'))
                    region_label_dict[sv_region].append((math.log10(end_std + 1), asm_caller, 'end_std'))

                    if int(row['SUPP']) == 3:
                        caller_total[caller_index] += 1
                        if start_std == 0:
                            start_std_dict['0'][caller_index] += 1
                        elif start_std <= 10:
                            start_std_dict['0,10'][caller_index] += 1
                        elif start_std <= 50:
                            start_std_dict['10,50'][caller_index] += 1
                        else:
                            start_std_dict['>50'][caller_index] += 1

                        if end_std == 0:
                            end_std_dict['0'][caller_index] += 1
                        elif end_std <= 10:
                            end_std_dict['0,10'][caller_index] += 1
                        elif end_std <= 50:
                            end_std_dict['10,50'][caller_index] += 1
                        else:
                            end_std_dict['>50'][caller_index] += 1

        for caller in CALLERS:
            caller_index = callers.index(caller)
            matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{plat}/{caller}.{aligner}.{plat}-concordant.info.tsv',sep='\t', header=[0])
            for idx, row in matched_info.iterrows():
                start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                region_label_dict[sv_region].append((math.log10(start_std + 1), caller, 'start_std'))
                region_label_dict[sv_region].append((math.log10(end_std + 1), caller, 'end_std'))

                if int(row['SUPP']) == 3:
                    caller_total[caller_index] += 1
                    if start_std == 0:
                        start_std_dict['0'][caller_index] += 1
                    elif start_std <= 10:
                        start_std_dict['0,10'][caller_index] += 1
                    elif start_std <= 50:
                        start_std_dict['10,50'][caller_index] += 1
                    else:
                        start_std_dict['>50'][caller_index] += 1

                    if end_std == 0:
                        end_std_dict['0'][caller_index] += 1
                    elif end_std <= 10:
                        end_std_dict['0,10'][caller_index] += 1
                    elif end_std <= 50:
                        end_std_dict['10,50'][caller_index] += 1
                    else:
                        end_std_dict['>50'][caller_index] += 1


        for l, shift in enumerate(shift_labels):

            start_shift = [i / j * 100 for i, j in zip(start_std_dict[shift], caller_total)]
            # end_shift = [i / j * 100 for i, j in zip(end_std_dict[shift], caller_total)]

            for i, pcrt in enumerate(start_shift):
                if i < len(CALLERS):
                    all_svs_bpstd.append((pcrt, PLATMAP[plat], callers[i], shift, 'Read-based'))
                else:
                    all_svs_bpstd.append((pcrt, PLATMAP[plat], callers[i], shift, 'Assembly-based'))


                if shift == '0,10' or shift == '0':
                    if i < len(CALLERS):
                        acc_svs_pcrt.append((pcrt, PLATMAP[plat], TOOLMAP[callers[i]], shift, 'Read-based'))
                        # acc_svs_pcrt.append((end_shift[i], PLATMAP[plat], callers[i], shift, 'Read-based'))
                    else:
                        acc_svs_pcrt.append((pcrt, PLATMAP[plat], TOOLMAP[callers[i]], shift, 'Assembly-based'))
                        # acc_svs_pcrt.append((end_shift[i], PLATMAP[plat], callers[i], shift, 'Assembly-based'))

            # if l == 0:
            #     leftbp_ax.bar(plat_xticks[k], start_shift, color=BPSHIFTCOLORS[l], width=barwidth, edgecolor='black', label=shift, hatch=hatch)
            #     rightbp_ax.bar(plat_xticks[k], end_shift, color=BPSHIFTCOLORS[l], width=barwidth, edgecolor='black', label=shift, hatch=hatch)
            # else:
            #     start_bottom = []
            #     end_bottom = []
            #     for m in range(0, l):
            #         start_bottom.append([i / j * 100 for i, j in zip(start_std_dict[shift_labels[m]], caller_total)])
            #         end_bottom.append([i / j * 100 for i, j in zip(end_std_dict[shift_labels[m]], caller_total)])
            #
            #     leftbp_ax.bar(plat_xticks[k], start_shift, bottom=[sum(x) for x in zip(*start_bottom)], color=BPSHIFTCOLORS[l],
            #                   width=barwidth, edgecolor='black', label=shift, hatch=hatch)
            #     rightbp_ax.bar(plat_xticks[k], end_shift, bottom=[sum(x) for x in zip(*end_bottom)], color=BPSHIFTCOLORS[l],
            #                    width=barwidth, edgecolor='black', label=shift, hatch=hatch)

    # for ax in [leftbp_ax, rightbp_ax]:
    #     ax.set_ylim(0, 100)
    #     ax.set_yticks(np.linspace(0, 100, 5))
    #     ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)])
    #     ax.set_xticks([r + (barwidth + 0.05) / 2 for r in r1])
    #     ax.set_xticklabels([TOOLMAP[caller] for caller in callers], rotation=90)
    #     ax.legend(handles=legends)
    #     ax.spines['left'].set_visible(False)
    #     ax.spines['right'].set_visible(False)
    #     ax.grid(axis='y', ls='--', color='grey', lw=2)
    #
    # fig.tight_layout()
    # fig1.tight_layout()
    #
    # fig.savefig(f'{Figure3}/platform_concordants_left_bpstd_of_callers.pdf')
    # fig1.savefig(f'{Figure3}/platform_concordants_right_bpstd_of_callers.pdf')

    df_all_svs = pd.DataFrame(all_svs_bpstd, columns=['pcrt', 'plat', 'caller', 'shift', 'stra'])
    fig2, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(6, 3))
    plat_order = ['HiFi', 'ONT']
    # plat_color = [PLATCOLORS[plat] for plat in plat_order]
    # plat_legends = [Patch(facecolor=PLATCOLORS[plat], label=plat) for plat in plat_order]

    stra_order = ['Assembly-based', 'Read-based']
    stra_colors = [STRACOLORS[ele.split('-')[0]] for ele in stra_order]

    # shift_order = ['0', '0,10']
    # shift_color = [SHIFTCOLORDICT[shift] for shift in shift_order]

    sns.barplot(data=df_all_svs[df_all_svs['plat'] == 'HiFi'], x='shift', y='pcrt',
                hue='stra', hue_order=stra_order, palette=stra_colors, capsize=.2, ax=axes[0])

    axes[0].set_title('HiFi', fontsize=13)

    sns.barplot(data=df_all_svs[df_all_svs['plat'] == 'ONT'], x='shift', y='pcrt',
                hue='stra', hue_order=stra_order, palette=stra_colors, capsize=.2, ax=axes[1])

    axes[1].set_title('ONT', fontsize=13)

    for i, ax in enumerate(axes):
        # ax.set_ylabel('Percent of SVs', fontsize=13)
        if i == 0:
            ax.set_ylabel('Percent of SVs', fontsize=13)
        else:
            ax.set_ylabel('')

        ax.set_xlabel('')
        ax.set_xticklabels(['0', '0~10', '10~50', '>50'], fontsize=12)

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
        ax.legend(title='')
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig2.tight_layout()

    df_acc_svs = pd.DataFrame(acc_svs_pcrt, columns=['pcrt', 'plat', 'caller', 'shift', 'stra'])
    fig3, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(4, 3))

    caller_order = [TOOLMAP[caller] for caller in callers]
    caller_color = [TOOLCOLORS[caller] for caller in callers]

    sns.stripplot(data=df_acc_svs[df_acc_svs['plat'] == 'HiFi'], x='shift', y='pcrt',
                hue='caller', hue_order=caller_order, palette=caller_color, size=7, ax=axes[0])

    axes[0].set_title('HiFi', fontsize=13)

    sns.stripplot(data=df_acc_svs[df_acc_svs['plat'] == 'ONT'], x='shift', y='pcrt',
                hue='caller', hue_order=caller_order, palette=caller_color, size=7, ax=axes[1])

    axes[1].set_title('ONT', fontsize=13)

    for i, ax in enumerate(axes):
        if i == 0:
            ax.set_ylabel('% of SVs', fontsize=13)
            ax.legend('',frameon=False)
        else:
            ax.set_ylabel('')
            ax.legend(loc='lower right', bbox_to_anchor=(1.04, 0))

        ax.set_xlabel('')
        ax.set_xticklabels(['0', '0~10'], fontsize=12)

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)

        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig3.tight_layout()

    plt.show()
    fig2.savefig(f'{Figure3}/platform_concordants_bpstd.pdf')
    fig3.savefig(f'{Figure3}/platform_concordants_bpstd_of_callers.pdf')