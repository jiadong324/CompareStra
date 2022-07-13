#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/4/6

'''

import math
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import numpy as np
import seaborn as sns
import vcf
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

def plot_bpstd_stacked(workdir, aligner):
    shift_labels = ['0', '0,10', '10,50', '>50']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'pav', 'svimasm']

    bpstd_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
    caller_total = [1 for k in range(len(callers))]


    for asm_caller in ASMCALLERS:
        caller_index = callers.index(asm_caller)
        matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
        for idx, row in matched_info.iterrows():
            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

            if int(row['SUPP']) == 6:
                caller_total[caller_index] += 1
                if start_std == 0:
                    bpstd_dict['0'][caller_index] += 1
                elif start_std <= 10:
                    bpstd_dict['0,10'][caller_index] += 1
                elif start_std <= 50:
                    bpstd_dict['10,50'][caller_index] += 1
                else:
                    bpstd_dict['>50'][caller_index] += 1

                # if end_std <= 10:
                #     bpstd_dict['0,10'][caller_index] += 1
                # elif end_std <= 50:
                #     bpstd_dict['10,50'][caller_index] += 1
                # elif end_std <= 100:
                #     bpstd_dict['50,100'][caller_index] += 1
                # else:
                #     bpstd_dict['>100'][caller_index] += 1

    for caller in CALLERS:
        caller_index = callers.index(caller)
        matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
        for idx, row in matched_info.iterrows():
            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])),row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

            if int(row['SUPP']) == 6:
                caller_total[caller_index] += 1
                if start_std == 0:
                    bpstd_dict['0'][caller_index] += 1
                elif start_std <= 10:
                    bpstd_dict['0,10'][caller_index] += 1
                elif start_std <= 50:
                    bpstd_dict['10,50'][caller_index] += 1
                else:
                    bpstd_dict['>50'][caller_index] += 1

                # if end_std <= 10:
                #     bpstd_dict['0,10'][caller_index] += 1
                # elif end_std <= 50:
                #     bpstd_dict['10,50'][caller_index] += 1
                # elif end_std <= 100:
                #     bpstd_dict['50,100'][caller_index] += 1
                # else:
                #     bpstd_dict['>100'][caller_index] += 1

    fig, ax = plt.subplots(1, 1, figsize=(5, 3))

    # leftbp_ax.set_title('Left breakpoint', fontsize=13)
    ax.set_ylabel('Percent of SVs', fontsize=13)
    xticks = np.arange(len(callers))
    shift_list = []
    for l, shift in enumerate(shift_labels):

        start_shift = [i / j * 100 for i, j in zip(bpstd_dict[shift], caller_total)]

        # for j, count in enumerate(bpstd_dict[shift]):
        #     pcrt = count / caller_total[j] * 100
        #     stra = 'Read-based'
        #     if callers[j] in ASMCALLERS:
        #         stra = 'Assembly-based'
        #     shift_list.append((pcrt, shift, callers[j], stra))


        if l == 0:
            ax.bar(xticks, start_shift, color=SHIFTCOLORDICT[shift_labels[l]], width=0.5, edgecolor='w', label=shift)
        else:
            start_bottom = []
            for m in range(0, l):
                start_bottom.append([i / j * 100 for i, j in zip(bpstd_dict[shift_labels[m]], caller_total)])
            ax.bar(xticks, start_shift, bottom=[sum(x) for x in zip(*start_bottom)], color=SHIFTCOLORDICT[shift_labels[l]],width=0.5, edgecolor='w', label=shift)

    # legends = [Line2D([0], [0], lw=2, color=TOOLCOLORS[caller], label=TOOLMAP[caller]) for caller in callers]

    # df_shift = pd.DataFrame(shift_list, columns=['pcrt', 'shift','caller', 'stra'])
    # hue_colors = [TOOLCOLORS[caller] for caller in callers]
    # sns.lineplot(data=df_shift, x='shift', y='pcrt', hue='caller', hue_order=callers, palette=hue_colors, marker='o', lw=2, ax=ax)
    #
    ax.set_ylim(0, 100)
    ax.set_yticks(np.linspace(0, 100, 5))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
    ax.set_xticks(xticks)
    ax.set_xticklabels([TOOLMAP[caller] for caller in callers], fontsize=13, rotation=90)
    ax.set_xlabel('')
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', ls='--', color='grey', lw=2)


    legend = [Patch(facecolor=SHIFTCOLORDICT['0'], label='0'), Patch(facecolor=SHIFTCOLORDICT['0,10'], label='0~10'),
              Patch(facecolor=SHIFTCOLORDICT['10,50'], label='10~50'), Patch(facecolor=SHIFTCOLORDICT['>50'], label='>50')]

    ax.legend(handles=legend)
    fig.tight_layout()

    plt.show()
    fig.savefig(f'{Figure2}/datasets_concordants_bpstd_of_callers.pdf')

def plot_bpstd(workdir, aligner):
    shift_labels = ['0', '0,10', '10,50', '>50']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']

    bpstd_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
    caller_total = [1 for k in range(len(callers))]


    for asm_caller in ASMCALLERS:
        caller_index = callers.index(asm_caller)
        for assembler in ASSEMBLER:
            matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{assembler}.dastsets-concordant.info.tsv', sep='\t', header=[0])
            for idx, row in matched_info.iterrows():
                start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                if int(row['SUPP']) == 6:
                    caller_total[caller_index] += 1
                    if start_std == 0:
                        bpstd_dict['0'][caller_index] += 1
                    elif start_std <= 10:
                        bpstd_dict['0,10'][caller_index] += 1
                    elif start_std <= 50:
                        bpstd_dict['10,50'][caller_index] += 1
                    else:
                        bpstd_dict['>50'][caller_index] += 1

                # if end_std <= 10:
                #     bpstd_dict['0,10'][caller_index] += 1
                # elif end_std <= 50:
                #     bpstd_dict['10,50'][caller_index] += 1
                # elif end_std <= 100:
                #     bpstd_dict['50,100'][caller_index] += 1
                # else:
                #     bpstd_dict['>100'][caller_index] += 1

    for caller in CALLERS:
        caller_index = callers.index(caller)
        matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
        for idx, row in matched_info.iterrows():
            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])),row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

            if int(row['SUPP']) == 6:
                caller_total[caller_index] += 1
                if start_std == 0:
                    bpstd_dict['0'][caller_index] += 1
                elif start_std <= 10:
                    bpstd_dict['0,10'][caller_index] += 1
                elif start_std <= 50:
                    bpstd_dict['10,50'][caller_index] += 1
                else:
                    bpstd_dict['>50'][caller_index] += 1

                # if end_std <= 10:
                #     bpstd_dict['0,10'][caller_index] += 1
                # elif end_std <= 50:
                #     bpstd_dict['10,50'][caller_index] += 1
                # elif end_std <= 100:
                #     bpstd_dict['50,100'][caller_index] += 1
                # else:
                #     bpstd_dict['>100'][caller_index] += 1


    shift_pcrt = []
    for l, shift in enumerate(shift_labels):
        for caller_idx, caller in enumerate(callers):
            shift_pcrt.append((bpstd_dict[shift][caller_idx] / caller_total[caller_idx] * 100, TOOLMAP[caller], shift))

    df_shift_pcrt = pd.DataFrame(shift_pcrt, columns=['pcrt', 'caller', 'shift'])
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))
    sns.stripplot(data=df_shift_pcrt, x='shift', y='pcrt', hue='caller', hue_order=[TOOLMAP[caller] for caller in callers],
                  palette=[TOOLCOLORS[caller] for caller in callers], size=7, ax=ax)

    ax.legend()
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.grid(axis='y', ls='--', color='grey', lw=2)
    ax.set_ylabel('% of concordant SVs', fontsize=13)

    ax.set_ylim(0, 60)
    ax.set_yticks(np.linspace(0, 60, 4))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 60, 4)], fontsize=12)

    ax.set_xlabel('')
    ax.set_xticks(np.arange(len(shift_labels)))
    ax.set_xticklabels(['0', '0~10', '10~50', '>50'], fontsize=13)

    plt.tight_layout()
    fig.savefig(f'{Figure2}/datasets_concordant_bpstd.pdf')
    plt.show()

def plot_bpstd_radar(workdir, aligner):
    shift_labels = ['0', '0,10', '10,50', '>50']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']

    bpstd_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
    caller_total = [1 for k in range(len(callers))]


    for asm_caller in ASMCALLERS:
        caller_index = callers.index(asm_caller)
        for assembler in ASSEMBLER:
            matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{assembler}.dastsets-concordant.info.tsv', sep='\t', header=[0])
            for idx, row in matched_info.iterrows():
                start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                if int(row['SUPP']) == 6:
                    caller_total[caller_index] += 1
                    if start_std == 0:
                        bpstd_dict['0'][caller_index] += 1
                    elif start_std <= 10:
                        bpstd_dict['0,10'][caller_index] += 1
                    elif start_std <= 50:
                        bpstd_dict['10,50'][caller_index] += 1
                    else:
                        bpstd_dict['>50'][caller_index] += 1

                # if end_std <= 10:
                #     bpstd_dict['0,10'][caller_index] += 1
                # elif end_std <= 50:
                #     bpstd_dict['10,50'][caller_index] += 1
                # elif end_std <= 100:
                #     bpstd_dict['50,100'][caller_index] += 1
                # else:
                #     bpstd_dict['>100'][caller_index] += 1

    for caller in CALLERS:
        caller_index = callers.index(caller)
        matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
        for idx, row in matched_info.iterrows():
            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])),row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

            if int(row['SUPP']) == 6:
                caller_total[caller_index] += 1
                if start_std == 0:
                    bpstd_dict['0'][caller_index] += 1
                elif start_std <= 10:
                    bpstd_dict['0,10'][caller_index] += 1
                elif start_std <= 50:
                    bpstd_dict['10,50'][caller_index] += 1
                else:
                    bpstd_dict['>50'][caller_index] += 1

                # if end_std <= 10:
                #     bpstd_dict['0,10'][caller_index] += 1
                # elif end_std <= 50:
                #     bpstd_dict['10,50'][caller_index] += 1
                # elif end_std <= 100:
                #     bpstd_dict['50,100'][caller_index] += 1
                # else:
                #     bpstd_dict['>100'][caller_index] += 1

    radar_data_dict = {'group': callers}

    for l, shift in enumerate(shift_labels):
        radar_data_dict[shift] = [i / j * 100 for i, j in zip(bpstd_dict[shift], caller_total)]


    df_radar = pd.DataFrame(radar_data_dict)
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
    ax.set_theta_zero_location("NW")

    # Draw one axe per variable + add labels
    plt.xticks(angles[:-1], categories, fontsize=13)

    # Draw ylabels
    ax.set_rlabel_position(270)
    plt.yticks([20, 40, 60], ["20%", "40%", "60%"], color="grey", size=12)
    plt.ylim(0, 60)

    for i, caller in enumerate(callers):
        values = df_radar.loc[i].drop('group').values.flatten().tolist()
        values += values[:1]
        ax.plot(angles, values, linewidth=2, linestyle='solid', label=caller, color=TOOLCOLORS[caller])
        ax.fill(angles, values, TOOLCOLORS[caller], alpha=0.1)

    # plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
    plt.tight_layout()
    fig.savefig(f'{Figure2}/concordant_bpstd_radar.pdf')
    plt.show()


def plot_bpstd_byregions(workdir, aligner):
    # region_label_dict = {'Simple Repeats':[], 'Repeat Masked':[], 'Segment Dup':[], 'Unique':[]}
    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']
    shift_labels = ['0', '0,10', '10,50', '>50']

    # start_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
    # end_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}

    start_region_label = {'Simple Repeats': {shift: [0 for k in range(len(callers))] for shift in shift_labels},
                          'Repeat Masked': {shift: [0 for k in range(len(callers))] for shift in shift_labels},
                          'Segment Dup': {shift: [0 for k in range(len(callers))] for shift in shift_labels},
                          'Unique': {shift: [0 for k in range(len(callers))] for shift in shift_labels}}

    # end_region_label = {'Simple Repeats': {shift: [0 for k in range(len(callers))] for shift in shift_labels},
    #                     'Repeat Masked': {shift: [0 for k in range(len(callers))] for shift in shift_labels},
    #                     'Segment Dup': {shift: [0 for k in range(len(callers))] for shift in shift_labels},
    #                     'Unique': {shift: [0 for k in range(len(callers))] for shift in shift_labels}}

    caller_region_total = {'Simple Repeats': [1 for k in range(len(callers))],
                           'Repeat Masked': [1 for k in range(len(callers))],
                           'Segment Dup': [1 for k in range(len(callers))],
                           'Unique': [1 for k in range(len(callers))]}


    for asm_caller in ASMCALLERS:
        caller_index = callers.index(asm_caller)
        for assembler in ASSEMBLER:
            matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{assembler}.dastsets-concordant.info.tsv', sep='\t',header=[0])

            for idx, row in matched_info.iterrows():
                start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(
                    float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                caller_region_total[sv_region][caller_index] += 1

                if start_std == 0:
                    start_region_label[sv_region]['0'][caller_index] += 1
                elif start_std <= 10:
                    start_region_label[sv_region]['0,10'][caller_index] += 1
                elif start_std <= 50:
                    start_region_label[sv_region]['10,50'][caller_index] += 1
                else:
                    start_region_label[sv_region]['>50'][caller_index] += 1

                # if end_std <= 10:
                #     end_region_label[sv_region]['0,10'][caller_index] += 1
                # elif end_std <= 50:
                #     end_region_label[sv_region]['10,50'][caller_index] += 1
                # elif end_std <= 100:
                #     end_region_label[sv_region]['50,100'][caller_index] += 1
                # else:
                #     end_region_label[sv_region]['>100'][caller_index] += 1

    for caller in CALLERS:
        caller_index = callers.index(caller)
        matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv',sep='\t', header=[0])
        for idx, row in matched_info.iterrows():
            start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), \
                                                           row['TYPE_MATCH'], abs(int(row['SVLEN'])), row[
                                                               'REGION_TYPE']

            caller_region_total[sv_region][caller_index] += 1

            if start_std == 0:
                start_region_label[sv_region]['0'][caller_index] += 1
            elif start_std <= 10:
                start_region_label[sv_region]['0,10'][caller_index] += 1
            elif start_std <= 50:
                start_region_label[sv_region]['10,50'][caller_index] += 1
            else:
                start_region_label[sv_region]['>50'][caller_index] += 1

            # if end_std <= 10:
            #     end_region_label[sv_region]['0,10'][caller_index] += 1
            # elif end_std <= 50:
            #     end_region_label[sv_region]['10,50'][caller_index] += 1
            # elif end_std <= 100:
            #     end_region_label[sv_region]['50,100'][caller_index] += 1
            # else:
            #     end_region_label[sv_region]['>100'][caller_index] += 1

    xticks = np.arange(len(callers))
    fig, axes = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(11, 3))
    legend = [Patch(facecolor=SHIFTCOLORDICT['0'], label='0'),
              Patch(facecolor=SHIFTCOLORDICT['0,10'], label='0~10'),
              Patch(facecolor=SHIFTCOLORDICT['10,50'], label='10~50'),
              Patch(facecolor=SHIFTCOLORDICT['>50'], label='>50')]
    for col_idx, region_label in enumerate(region_labels):
        this_ax = axes[col_idx]
        this_region_dict = start_region_label[region_label]
        this_ax.set_title(region_label)

        this_ax.set_ylim(0, 100)
        this_ax.set_yticks(np.linspace(0, 100, 5))
        this_ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)

        for l, shift in enumerate(shift_labels):

            this_shift = [i / j * 100 for i, j in zip(this_region_dict[shift], caller_region_total[region_label])]

            if l == 0:
                this_ax.bar(xticks, this_shift, color=SHIFTCOLORDICT[shift_labels[l]], width=0.5, edgecolor='w', label=shift)
            else:
                bottom = []
                for m in range(0, l):
                    bottom.append([i / j * 100 for i, j in
                                   zip(this_region_dict[shift_labels[m]], caller_region_total[region_label])])

                bottom_sum = [sum(x) for x in zip(*bottom)]
                this_ax.bar(xticks, this_shift, bottom=bottom_sum, color=SHIFTCOLORDICT[shift_labels[l]], width=0.5, edgecolor='w',
                            label=shift)

        this_ax.set_xticks(xticks)
        this_ax.set_xticklabels([TOOLMAP[caller] for caller in callers], rotation=90, fontsize=13)

        this_ax.spines['left'].set_visible(False)
        this_ax.spines['right'].set_visible(False)
        this_ax.grid(axis='y', ls='--', color='grey', lw=2)

        if col_idx == 0:
            # this_ax.legend()
            this_ax.legend(handles=legend)
            this_ax.set_ylabel('Percent of SVs', fontsize=13)

    fig.tight_layout()
    plt.show()
    fig.savefig(f'{Figure2}/datasets_concordant_bpstd_byregions.pdf')


def plot_left_right_bpstd(workdir, aligners):

    shift_labels = ['0,10', '10,50', '50,100', '>100']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'pav', 'svimasm']

    start_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
    end_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
    caller_total = [1 for k in range(len(callers))]

    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']
    region_label_dict = {region: [] for region in region_labels}

    for aligner in aligners:
        if aligner == 'minimap2':
            for asm_caller in ASMCALLERS:
                caller_index = callers.index(asm_caller)
                matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
                for idx, row in matched_info.iterrows():
                    start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                    region_label_dict[sv_region].append((math.log10(start_std + 1), asm_caller, 'start_std'))
                    region_label_dict[sv_region].append((math.log10(end_std + 1), asm_caller, 'end_std'))

                    if int(row['SUPP']) == 6:
                        caller_total[caller_index] += 1
                        if start_std <= 10:
                            start_std_dict['0,10'][caller_index] += 1
                        elif start_std <= 50:
                            start_std_dict['10,50'][caller_index] += 1
                        elif start_std <= 100:
                            start_std_dict['50,100'][caller_index] += 1
                        else:
                            start_std_dict['>100'][caller_index] += 1

                        if end_std <= 10:
                            end_std_dict['0,10'][caller_index] += 1
                        elif end_std <= 50:
                            end_std_dict['10,50'][caller_index] += 1
                        elif end_std <= 100:
                            end_std_dict['50,100'][caller_index] += 1
                        else:
                            end_std_dict['>100'][caller_index] += 1


        for caller in CALLERS:
            caller_index = callers.index(caller)
            matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
            for idx, row in matched_info.iterrows():
                start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                region_label_dict[sv_region].append((math.log10(start_std + 1), caller, 'start_std'))
                region_label_dict[sv_region].append((math.log10(end_std + 1), caller, 'end_std'))

                if int(row['SUPP']) == 6:
                    caller_total[caller_index] += 1
                    if start_std <= 10:
                        start_std_dict['0,10'][caller_index] += 1
                    elif start_std <= 50:
                        start_std_dict['10,50'][caller_index] += 1
                    elif start_std <= 100:
                        start_std_dict['50,100'][caller_index] += 1
                    else:
                        start_std_dict['>100'][caller_index] += 1

                    if end_std <= 10:
                        end_std_dict['0,10'][caller_index] += 1
                    elif end_std <= 50:
                        end_std_dict['10,50'][caller_index] += 1
                    elif end_std <= 100:
                        end_std_dict['50,100'][caller_index] += 1
                    else:
                        end_std_dict['>100'][caller_index] += 1

    xticks = np.arange(len(callers))
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))

    leftbp_ax = axes[0]
    # leftbp_ax.set_title('Left breakpoint', fontsize=13)
    leftbp_ax.set_ylabel('Left breakpoint std.', fontsize=13)
    rightbp_ax = axes[1]
    # rightbp_ax.set_title('Right breakpoint', fontsize=13)
    rightbp_ax.set_ylabel('Right breakpoint std.', fontsize=13)

    for l, shift in enumerate(shift_labels):

        start_shift = [i / j * 100 for i, j in zip(start_std_dict[shift], caller_total)]
        end_shift = [i / j * 100 for i, j in zip(end_std_dict[shift], caller_total)]

        if l == 0:
            leftbp_ax.bar(xticks, start_shift, color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)
            rightbp_ax.bar(xticks, end_shift, color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)
        else:
            start_bottom = []
            end_bottom = []
            for m in range(0, l):
                start_bottom.append([i / j * 100 for i, j in zip(start_std_dict[shift_labels[m]], caller_total)])
                end_bottom.append([i / j * 100 for i, j in zip(end_std_dict[shift_labels[m]], caller_total)])

            leftbp_ax.bar(xticks, start_shift, bottom=[sum(x) for x in zip(*start_bottom)], color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)
            rightbp_ax.bar(xticks, end_shift, bottom=[sum(x) for x in zip(*end_bottom)], color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)

    for ax in axes:
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)], fontsize=12)
        ax.set_xticks(xticks)
        ax.set_xticklabels([TOOLMAP[caller] for caller in callers], rotation=90, fontsize=13)
        ax.legend()
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', ls='--', color='grey', lw=2)


    fig.tight_layout()

    plt.show()
    fig.savefig(f'{Figure2}/datasets_concordants_left_right_bpstd_of_callers.pdf')


def plot_inacc_svs(workdir, aligner, overlap_pcrt):
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']

    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']
    inacc_svlen = {'Simple Repeats': [], 'Repeat Masked': [], 'Segment Dup': [], 'Unique': []}

    inacc_stras = {'Simple Repeats': [{}, {}], 'Repeat Masked': [{}, {}], 'Segment Dup': [{}, {}], 'Unique': [{}, {}]}

    inacc_svtypes = {'Simple Repeats': {}, 'Repeat Masked': {}, 'Segment Dup': {}, 'Unique': {}}
    inacc_regions = [[0 for i in range(len(callers))] for j in range(len(region_labels))]

    asm_concordants_num = 0
    read_concordants_num = 0


    for asm_caller in ASMCALLERS:
        caller_idx = callers.index(asm_caller)
        for assembler in ASSEMBLER:
            matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{assembler}.dastsets-concordant.info.tsv', sep='\t',header=[0])

            for idx, row in matched_info.iterrows():
                start_std, end_std, svtype, svlen, sv_region, rppcrt = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], \
                                                               abs(int(row['SVLEN'])), row['REGION_TYPE'], float(row['PCRT'])

                if int(row['SUPP']) == 6:
                    asm_concordants_num += 1
                    region_idx = region_labels.index(sv_region)
                    if start_std > 50 and end_std > 50:
                        if rppcrt >= overlap_pcrt:
                            inacc_svlen[sv_region].append((svlen + 1, svtype, asm_caller, 'Assembly-based'))
                            inacc_regions[region_idx][caller_idx] += 1
                        else:
                            inacc_svlen['Unique'].append((svlen + 1, svtype, asm_caller, 'Assembly-based'))
                            inacc_regions[3][caller_idx] += 1

                        if svtype in inacc_stras[sv_region][0]:
                            inacc_stras[sv_region][0][svtype] += 1
                        else:
                            inacc_stras[sv_region][0][svtype] = 1

                        if svtype in inacc_svtypes[sv_region]:
                            inacc_svtypes[sv_region][svtype] += 1
                        else:
                            inacc_svtypes[sv_region][svtype] = 1

    for caller in CALLERS:
        caller_idx = callers.index(caller)
        matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])

        for idx, row in matched_info.iterrows():
            start_std, end_std, svtype, svlen, sv_region, rppcrt = abs(float(row['START_STD'])), abs(float(row['END_STD'])), \
                                                           row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE'], float(row['PCRT'])
            region_idx = region_labels.index(sv_region)
            if int(row['SUPP']) == 6:
                read_concordants_num += 1
                if start_std > 50 and end_std > 50:
                    if rppcrt >= overlap_pcrt:
                        inacc_svlen[sv_region].append((svlen + 1, svtype, caller, 'Read-based'))
                        inacc_regions[region_idx][caller_idx] += 1
                    else:
                        inacc_svlen['Unique'].append((svlen + 1, svtype, caller, 'Read-based'))
                        inacc_regions[3][caller_idx] += 1

                    if svtype in inacc_stras[sv_region][1]:
                        inacc_stras[sv_region][1][svtype] += 1
                    else:
                        inacc_stras[sv_region][1][svtype] = 1

                    if svtype in inacc_svtypes[sv_region]:
                        inacc_svtypes[sv_region][svtype] += 1
                    else:
                        inacc_svtypes[sv_region][svtype] = 1

    for region, inacc_count in inacc_stras.items():
        print(f'{region} \n\tAssembly: {inacc_count[0]} Read: {inacc_count[1]}')
    print(asm_concordants_num, read_concordants_num)


    fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    xticks = np.arange(len(callers))
    for i, region in enumerate(region_labels):
        this_data = inacc_regions[i]
        # if i == 0:
        #     ax.bar(xticks, this_data, label=region, color=REGIONCOLORS[region])
        # else:
        #     bottoms = []
        #     for j in range(0,i):
        #         bottoms.append(inacc_regions[j])
        #     bottom_sum = [sum(x) for x in zip(*bottoms)]
        #     ax.bar(xticks, this_data, label=region, color=REGIONCOLORS[region], bottom=bottom_sum)
        ax.plot(xticks, this_data, color=REGIONCOLORS[region], label=region, marker=REGIONMARKERS[region], lw=2, ms=7)

    ax.legend(ncol=2)
    ax.set_xticks(xticks)
    ax.set_xticklabels([TOOLMAP[caller] for caller in callers], fontsize=13, rotation=90)
    ax.set_ylim(0, 1500)
    ax.set_yticks(np.linspace(0, 1500, 4))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 1500, 4)], fontsize=12)
    ax.set_ylabel('Number of SVs', fontsize=13)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.grid(axis='y', ls='--', color='grey', lw=2)

    fig.tight_layout()

    inacc_types_list = []
    for region, svtype_dict in inacc_svtypes.items():
        for svtype, count in svtype_dict.items():
            inacc_types_list.append((svtype, count, region))
    df_inacc_types = pd.DataFrame(inacc_types_list, columns=['svtype', 'count', 'region'])

    fig1, ax1 = plt.subplots(1, 1, figsize=(5, 3))
    hue_order = ['DEL', 'INS', 'DUP', 'INV', 'BND', 'DUP:TANDEM']
    hue_color = [SVTYPECOLORS[svtype] for svtype in hue_order]
    sns.barplot(data=df_inacc_types, y='count', x='region', hue='svtype', hue_order=hue_order, palette=hue_color, ax=ax1)

    ax1.set_ylabel('# of SV', fontsize=13)
    ax1.set_yscale('log')
    ax1.set_xticks(np.arange(len(region_labels)))
    ax1.set_xticklabels(region_labels, fontsize=13, rotation=90)
    ax1.set_xlabel('')
    ax1.legend(title='', loc='lower right')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    fig1.tight_layout()


    fig2, ax2 = plt.subplots(1, 1, figsize=(5, 4))
    hue_order = ['Assembly-based', 'Read-based']
    hue_color = [STRACOLORS[ele.split('-')[0]] for ele in hue_order]
    legends = [Patch(facecolor=STRACOLORS['Assembly'], label='Assembly-based'), Patch(facecolor=STRACOLORS['Read'], label='Read-based')]

    df_unique_regions = pd.DataFrame(inacc_svlen['Simple Repeats'], columns=['svlen', 'svtype', 'caller', 'stra'])
    sns.histplot(data=df_unique_regions, x='svlen', hue='stra', hue_order=hue_order,
                 palette=hue_color, kde=True, log_scale=True, ax=ax2)

    # ax2.set_title('Simple Repeat regions')
    ax2.set_xlabel('SV size', fontsize=13)
    ax2.set_ylabel('Number of SV', fontsize=13)
    ax2.set_xlim(10, 100000)

    ax2.set_ylim(0, 300)
    ax2.set_yticks(np.linspace(0, 300, 4))
    ax2.set_yticklabels([int(val) for val in np.linspace(0, 300, 4)], fontsize=12)
    ax2.legend(handles=legends)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    fig2.tight_layout()

    fig3, ax3 = plt.subplots(1, 1, figsize=(5, 3))
    caller_colors = [TOOLCOLORS[caller] for caller in callers]
    caller_legend = [Line2D([0], [0], lw=2, label=TOOLMAP[caller], color=TOOLCOLORS[caller]) for caller in callers]
    sns.kdeplot(data=df_unique_regions, x='svlen', hue='caller', hue_order=callers, palette=caller_colors, log_scale=True, ax=ax3)

    # ax3.set_title('Unique regions')
    ax3.set_xlabel('')
    ax3.set_ylabel('Density', fontsize=13)
    ax3.set_xlim(10, 100000)
    ax3.set_ylim(0, 0.4)
    ax3.set_yticks(np.linspace(0, 0.4, 3))
    ax3.set_yticklabels([0, 0.2, 0.4], fontsize=12)
    ax3.legend(handles=caller_legend)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    fig3.tight_layout()


    plt.show()
    if overlap_pcrt == 50:
        fig.savefig(f'{Figure2}/datasets_concordant_bpinacc_highly_reps.pdf')
    else:
        fig.savefig(f'{Figure2}/datasets_concordant_bpinacc_reps.pdf')

    fig1.savefig(f'{Figure2}/datasets_concordant_bpinacc_svtypes.pdf')
    fig2.savefig(f'{Figure2}/datasets_concordant_bpinacc_svsize.pdf')
    fig3.savefig(f'{Figure2}/datasets_concordant_bpinacc_svsize_kde.pdf')

def plot_ins_bpstd(workdir, aligners):

    shift_labels = ['0,10', '10,50', '50,100', '>100']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'pav', 'svimasm']

    start_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
    end_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
    caller_total = [1 for k in range(len(callers))]

    for aligner in aligners:
        if aligner == 'minimap2':
            for asm_caller in ASMCALLERS:
                caller_index = callers.index(asm_caller)
                matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
                for idx, row in matched_info.iterrows():
                    start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                    if svtype != 'INS':
                        continue

                    caller_total[caller_index] += 1

                    if start_std <= 10:
                        start_std_dict['0,10'][caller_index] += 1
                    elif start_std <= 50:
                        start_std_dict['10,50'][caller_index] += 1
                    elif start_std <= 100:
                        start_std_dict['50,100'][caller_index] += 1
                    else:
                        start_std_dict['>100'][caller_index] += 1

                    if end_std <= 10:
                        end_std_dict['0,10'][caller_index] += 1
                    elif end_std <= 50:
                        end_std_dict['10,50'][caller_index] += 1
                    elif end_std <= 100:
                        end_std_dict['50,100'][caller_index] += 1
                    else:
                        end_std_dict['>100'][caller_index] += 1

        for caller in CALLERS:
            caller_index = callers.index(caller)
            matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
            for idx, row in matched_info.iterrows():
                start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                if svtype != 'INS':
                    continue

                caller_total[caller_index] += 1
                if start_std <= 10:
                    start_std_dict['0,10'][caller_index] += 1
                elif start_std <= 50:
                    start_std_dict['10,50'][caller_index] += 1
                elif start_std <= 100:
                    start_std_dict['50,100'][caller_index] += 1
                else:
                    start_std_dict['>100'][caller_index] += 1

                if end_std <= 10:
                    end_std_dict['0,10'][caller_index] += 1
                elif end_std <= 50:
                    end_std_dict['10,50'][caller_index] += 1
                elif end_std <= 100:
                    end_std_dict['50,100'][caller_index] += 1
                else:
                    end_std_dict['>100'][caller_index] += 1

    xticks = np.arange(len(callers))
    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 4))

    leftbp_ax = axes[0]
    leftbp_ax.set_title('Left breakpoint', fontsize=13)
    rightbp_ax = axes[1]
    rightbp_ax.set_title('Right breakpoint', fontsize=13)

    for l, shift in enumerate(shift_labels):

        start_shift = [i / j * 100 for i, j in zip(start_std_dict[shift], caller_total)]
        end_shift = [i / j * 100 for i, j in zip(end_std_dict[shift], caller_total)]

        if l == 0:
            leftbp_ax.bar(xticks, start_shift, color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)
            rightbp_ax.bar(xticks, end_shift, color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)
        else:
            start_bottom = []
            end_bottom = []
            for m in range(0, l):
                start_bottom.append([i / j * 100 for i, j in zip(start_std_dict[shift_labels[m]], caller_total)])
                end_bottom.append([i / j * 100 for i, j in zip(end_std_dict[shift_labels[m]], caller_total)])

            leftbp_ax.bar(xticks, start_shift, bottom=[sum(x) for x in zip(*start_bottom)], color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)
            rightbp_ax.bar(xticks, end_shift, bottom=[sum(x) for x in zip(*end_bottom)], color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)

    for ax in axes:
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)])
        ax.set_xticks(xticks)
        ax.set_xticklabels([TOOLMAP[caller] for caller in callers], rotation=90)
        ax.legend()

    fig.tight_layout()
    plt.show()
    fig.savefig(f'{Figure2}/dataset_concordants_ins_bpstd.pdf')


def plot_del_bpstd(workdir, aligners):

    shift_labels = ['0,10', '10,50', '50,100', '>100']
    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'pav', 'svimasm']

    start_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
    end_std_dict = {shift: [0 for k in range(len(callers))] for shift in shift_labels}
    caller_total = [1 for k in range(len(callers))]


    for aligner in aligners:
        if aligner == 'minimap2':
            for asm_caller in ASMCALLERS:
                caller_index = callers.index(asm_caller)
                matched_info = pd.read_csv(f'{workdir}/assm_dataset_repro/{asm_caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
                for idx, row in matched_info.iterrows():
                    start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']

                    if svtype != 'DEL':
                        continue

                    caller_total[caller_index] += 1

                    if start_std <= 10:
                        start_std_dict['0,10'][caller_index] += 1
                    elif start_std <= 50:
                        start_std_dict['10,50'][caller_index] += 1
                    elif start_std <= 100:
                        start_std_dict['50,100'][caller_index] += 1
                    else:
                        start_std_dict['>100'][caller_index] += 1

                    if end_std <= 10:
                        end_std_dict['0,10'][caller_index] += 1
                    elif end_std <= 50:
                        end_std_dict['10,50'][caller_index] += 1
                    elif end_std <= 100:
                        end_std_dict['50,100'][caller_index] += 1
                    else:
                        end_std_dict['>100'][caller_index] += 1

        for caller in CALLERS:
            caller_index = callers.index(caller)
            matched_info = pd.read_csv(f'{workdir}/read_dataset_repro/{caller}.{aligner}.dastsets-concordant.info.tsv', sep='\t', header=[0])
            for idx, row in matched_info.iterrows():
                start_std, end_std, svtype, svlen, sv_region = abs(float(row['START_STD'])), abs(float(row['END_STD'])), row['TYPE_MATCH'], abs(int(row['SVLEN'])), row['REGION_TYPE']
                if svtype != 'DEL':
                    continue

                caller_total[caller_index] += 1
                if start_std <= 10:
                    start_std_dict['0,10'][caller_index] += 1
                elif start_std <= 50:
                    start_std_dict['10,50'][caller_index] += 1
                elif start_std <= 100:
                    start_std_dict['50,100'][caller_index] += 1
                else:
                    start_std_dict['>100'][caller_index] += 1

                if end_std <= 10:
                    end_std_dict['0,10'][caller_index] += 1
                elif end_std <= 50:
                    end_std_dict['10,50'][caller_index] += 1
                elif end_std <= 100:
                    end_std_dict['50,100'][caller_index] += 1
                else:
                    end_std_dict['>100'][caller_index] += 1

    xticks = np.arange(len(callers))
    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 4))

    leftbp_ax = axes[0]
    leftbp_ax.set_title('Left breakpoint', fontsize=13)
    rightbp_ax = axes[1]
    rightbp_ax.set_title('Right breakpoint', fontsize=13)

    for l, shift in enumerate(shift_labels):

        start_shift = [i / j * 100 for i, j in zip(start_std_dict[shift], caller_total)]
        end_shift = [i / j * 100 for i, j in zip(end_std_dict[shift], caller_total)]

        if l == 0:
            leftbp_ax.bar(xticks, start_shift, color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)
            rightbp_ax.bar(xticks, end_shift, color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)
        else:
            start_bottom = []
            end_bottom = []
            for m in range(0, l):
                start_bottom.append([i / j * 100 for i, j in zip(start_std_dict[shift_labels[m]], caller_total)])
                end_bottom.append([i / j * 100 for i, j in zip(end_std_dict[shift_labels[m]], caller_total)])

            leftbp_ax.bar(xticks, start_shift, bottom=[sum(x) for x in zip(*start_bottom)], color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)
            rightbp_ax.bar(xticks, end_shift, bottom=[sum(x) for x in zip(*end_bottom)], color=BPSHIFTCOLORS[l], width=0.5, edgecolor='w', label=shift)

    for ax in axes:
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)])
        ax.set_xticks(xticks)
        ax.set_xticklabels([TOOLMAP[caller] for caller in callers], rotation=90)
        ax.legend()

    fig.tight_layout()
    plt.show()
    fig.savefig(f'{Figure2}/dataset_concordants_del_bpstd.pdf')