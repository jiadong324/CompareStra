#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2021/10/26

'''
import os
import numpy as np
import pandas as pd
import math
import matplotlib.pylab as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from helpers.Constant import *
from helpers.Functions import *


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



def plot_missed_pav_feature(workdir, samples):

    aligners = ['minimap2', 'ngmlr']
    platforms = ['hifi', 'ont']

    plot_rptypes = ['VNTR', 'STR', 'LTR', 'LINE', 'SINE']

    missed_pav_calls = []
    missed_pav_rptypes = []
    missed_pav_svtypes = []
    miseed_pav_mapq = []


    for platform in platforms:
        for aligner_idx in range(len(aligners)):
            aligner = aligners[aligner_idx]

            for caller_idx in range(len(CALLERS)):
                caller = CALLERS[caller_idx]

                total = 0
                rptypes = {'VNTR': 0, 'STR': 0, 'LTR': 0, 'SINE': 0, 'LINE': 0, 'Others': 0}
                # rptypes = {'Simple_repeats': 0, 'Mobile_element': 0, 'Others': 0}
                svtypes = {}
                mapq_tag = {'No_reads': 0, 'Low_mapq': 0, 'High_mapq': 0}

                for sample in samples:
                    annot_pavs = f'{workdir}/{sample}/{caller}/bed/{sample}.{platform}.{aligner}.{caller}.classified.missed.pavs.tsv'
                    df_annot_pavs = pd.read_csv(annot_pavs, sep='\t', names=['chrom', 'start', 'end', 'svtype', 'svlen', 'rptype', 'rppcrt', 'mapq', 'signum', 'nsignum', 'caller', 'platform', 'aligner'])

                    total += len(df_annot_pavs)

                    for idx, row in df_annot_pavs.iterrows():
                        missed_pav_calls.append((row['svlen'], caller, row['svtype'], f'{platform}-{aligner}'))

                        if row['svtype'] in svtypes:
                            svtypes[row['svtype']] += 1
                        else:
                            svtypes[row['svtype']] = 1

                        rptype = row['rptype']

                        if rptype not in plot_rptypes:
                            rptype = 'Others'

                        rptypes[rptype] += 1

                        mapq = int(row['mapq'])
                        if mapq == -1:
                            mapq_tag['No_reads'] += 1
                        elif mapq < 20:
                            mapq_tag['Low_mapq'] += 1
                        else:
                            mapq_tag['High_mapq'] += 1


                for rptype, count in rptypes.items():
                    missed_pav_rptypes.append((rptype, count * 100 / total, caller, f'{platform}-{aligner}'))

                for svtype, count in svtypes.items():
                    missed_pav_svtypes.append((svtype, count * 100 / total, caller, platform, aligner))

                for mapq, count in mapq_tag.items():
                    miseed_pav_mapq.append((mapq, count * 100 / total, caller, f'{platform}-{aligner}'))

    df_missed_pavs = pd.DataFrame(missed_pav_calls, columns=['svlen', 'caller', 'svtype', 'platform-aligner'])
    df_missed_pavs_rptypes = pd.DataFrame(missed_pav_rptypes, columns=['rptype', 'pcrt', 'caller', 'platform-aligner'])
    df_missed_pavs_mapq = pd.DataFrame(miseed_pav_mapq, columns=['mapq', 'pcrt', 'caller', 'platform-aligner'])
    df_missed_pavs_svtypes = pd.DataFrame(missed_pav_svtypes, columns=['svtype', 'pcrt', 'caller', 'platform', 'aligner'])

    kdecolors = [TOOLCOLORS[caller] for caller in CALLERS]
    custom_params = {"axes.spines.right": False, "axes.spines.top": False}

    sns.set_theme(style="ticks", rc=custom_params, font="Arial", font_scale=1.0)

    # fig = plt.figure(figsize=(6, 4))
    # ax = fig.add_subplot(1,1,1)
    #
    # sns.kdeplot(data=df_missed_pavs, x='svlen', hue='caller', log_scale=True, ax=ax, palette=kdecolors)
    # ax.set_ylim(0, 0.3)
    # ax.set_yticks(np.linspace(0, 0.3, 4))
    # ax.set_yticklabels([round(val, 2) for val in np.linspace(0, 0.3, 4)])
    # plt.tight_layout()
    #
    # fig.savefig(f'{workdir}/Figures/Figure6/samples.missed.com.size.pdf')

    fig1 = plt.figure(figsize=(6, 4))
    ax = fig1.add_subplot(1, 1, 1)

    hue_order = ['DEL', 'INS', 'INV']
    svtype_color = [SVTYPECOLORS[svtype] for svtype in hue_order]
    sns.barplot(data=df_missed_pavs_svtypes, x='caller', y='pcrt', hue='svtype', edgecolor='w', hue_order=hue_order, ax=ax, palette=svtype_color)
    ax.set_ylim(0, 100)
    ax.set_yticks(np.linspace(0, 100, 5))
    ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)])

    ax.set_xticks(np.arange(len(CALLERS)))
    ax.set_xticklabels([TOOLMAP[caller] for caller in CALLERS])

    ax.set_ylabel('')
    ax.set_xlabel('')

    # ax1 = ax.twinx()
    # ax1.set_ylim(1e0, 1e4)
    # ax1.set_yscale('log')

    plt.tight_layout()

    # fig1.savefig(f'{workdir}/Figures/Figure4/samples.missed.com.svtypes.pdf')
    fig1.savefig(f'{workdir}/thesis_figures/Figure7/samples.missed.pav.svtypes.pdf')

    fig2 = plt.figure(figsize=(12, 4))

    caller_legends = [Patch(facecolor=TOOLCOLORS[caller], edgecolor='white', label=TOOLMAP[caller]) for caller in CALLERS]

    ax1 = fig2.add_subplot(1,2,1)
    sns.boxplot(data=df_missed_pavs_rptypes, x='rptype', y='pcrt', hue='caller', ax=ax1, palette=kdecolors)
    ax1.set_ylabel('')
    ax1.set_xlabel('')
    ax2 = fig2.add_subplot(1,2,2)
    sns.boxplot(data=df_missed_pavs_mapq, x='mapq', y='pcrt', hue='caller', ax=ax2, palette=kdecolors)

    for ax in [ax1, ax2]:
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 5))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 5)])
        ax.legend(handles=caller_legends)

    plt.tight_layout()
    plt.show()
    # fig2.savefig(f'{workdir}/Figures/Figure6/samples.missed.com.features.pdf')
    fig2.savefig(f'{workdir}/thesis_figures/Figure7/samples.missed.pav.features.pdf')


def plot_sig_tags_in_high_mapq(workdir, samples):

    sns.set_theme(style="ticks", font="Arial", font_scale=1.0)

    legends = [Line2D([0], [0], marker='o', label='Has_sig', color='#9ac0cd'),
               Line2D([0], [0], marker='X', label='No_sig', color='#fdb716'),
               Line2D([0], [0], marker='*', label='Has_nsig', color='#146533')]


    aligners = ['minimap2', 'ngmlr']
    platforms = ['hifi', 'ont']

    fig = plt.figure(figsize=(10, 7))
    xticks = np.arange(len(CALLERS))


    fig_idx = 1
    for platform in platforms:
        caller_sig_pcrt = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
        caller_nsig_pcrt = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]

        for aligner_idx in range(len(aligners)):

            aligner = aligners[aligner_idx]

            lack_sig = []
            lack_sig_pcrt = []

            enough_sig = []
            enough_sig_pcrt = []

            enough_nsig = []
            enough_nsig_pcrt = []

            for caller_idx in range(len(CALLERS)):
                caller = CALLERS[caller_idx]

                sig_tag = {'Lack_sig_reads': 0, 'Enough_sig_reads': {'Lack_nsig_reads': 0, 'Enough_nsig_reads': 0}}

                total_count = 0
                for sample in samples:
                    annot_pavs = f'{workdir}/{sample}/{caller}/bed/{sample}.{platform}.{aligner}.{caller}.classified.missed.pavs.tsv'
                    df_annot_pavs = pd.read_csv(annot_pavs, sep='\t', names=['chrom', 'start', 'end', 'svtype', 'svlen', 'rptype', 'rppcrt', 'mapq', 'signum', 'nsignum', 'caller', 'platform', 'aligner'])

                    for idx, row in df_annot_pavs.iterrows():
                        sig_num = int(row['signum'])
                        nsig_num = int(row['nsignum'])
                        mapq = int(row['mapq'])
                        if mapq >= 20:
                            total_count += 1
                            if sig_num < 2:
                                sig_tag['Lack_sig_reads'] += 1
                            else:
                                if nsig_num < 2:
                                    sig_tag['Enough_sig_reads']['Lack_nsig_reads'] += 1
                                else:
                                    sig_tag['Enough_sig_reads']['Enough_nsig_reads'] += 1

                for tag, count in sig_tag.items():
                    if tag == 'Lack_sig_reads':
                        lack_sig.append(count / 3)
                        lack_sig_pcrt.append(count / total_count)

                    elif tag == 'Enough_sig_reads':
                        enough_nsig.append(count['Enough_nsig_reads'] / 3)
                        enough_nsig_pcrt.append(count['Enough_nsig_reads'] / total_count)


                        caller_nsig_pcrt[caller_idx][aligner_idx] += count['Enough_nsig_reads'] / total_count

                        enough_sig.append((count['Enough_nsig_reads'] + count['Lack_nsig_reads']) / 3)
                        enough_sig_pcrt.append((count['Enough_nsig_reads'] + count['Lack_nsig_reads']) / total_count)

                        caller_sig_pcrt[caller_idx][aligner_idx] += (count['Enough_nsig_reads'] + count['Lack_nsig_reads']) / total_count


            ax = fig.add_subplot(2,2,fig_idx)

            ax.set_title(f'{platform}-{aligner}')
            ax.grid(color='grey', linestyle=(0, (5, 5)), lw=1, axis='y')

            ax.bar(xticks, enough_sig, color='#9ac0cd')
            ax.bar(xticks, lack_sig, bottom=enough_sig,color='#fdb716')
            ax.bar(xticks, enough_nsig, edgecolor='#146533', hatch='//', color='None', lw=2)

            ax.set_ylim(0, 14000)
            ax.set_yticks(np.linspace(0 , 14000, 5))
            if fig_idx == 1 or fig_idx == 3:
                ax.set_yticklabels([0, 35, 70, 105, 140])
            else:
                ax.set_yticklabels([])

            ax.set_xticks(xticks)
            ax.set_xticklabels(CALLERS)

            ax1 = ax.twinx()
            ax1.plot(xticks, enough_sig_pcrt, color='#9ac0cd', marker='o', lw=2, markersize=8)
            ax1.plot(xticks, lack_sig_pcrt, color='#fdb716', marker='X', lw=2, markersize=10)
            ax1.plot(xticks, enough_nsig_pcrt, color='#146533', marker='*', lw=2, markersize=10)

            ax1.set_ylim(0, 1)
            ax1.set_yticks(np.linspace(0, 1, 5))
            if fig_idx == 2 or fig_idx == 4:
                ax1.set_yticklabels([0, 0.25, 0.50, 0.75, 1])
            else:
                ax1.set_yticklabels([])

            if fig_idx == 4:
                ax.legend(handles=legends)


            fig_idx += 1


    plt.tight_layout()
    plt.show()
    # fig.savefig(f'{workdir}/Figures/Figure5/samples.missed_pavs.sigs.pdf')

def scatter_sig_tags_in_high_mapq(workdir, samples):

    sns.set_theme(style="ticks", font="Arial", font_scale=1.0)
    long_read_callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']


    legends = [Line2D([0], [0], marker='o', label='Has_sig', color='#9ac0cd'),
               Line2D([0], [0], marker='X', label='No_sig', color='#fdb716'),
               Line2D([0], [0], marker='*', label='Has_nsig', color='#146533')]


    aligners = ['minimap2', 'ngmlr']
    platforms = ['hifi', 'ont']

    plt.rcParams['hatch.linewidth'] = 2

    fig1 = plt.figure(figsize=(5, 4))
    fig1_ax = fig1.add_subplot(1,1,1)

    fig1_ax.set_ylim(0.6, 1)
    fig1_ax.set_yticks(np.linspace(0.6, 1, 5))
    fig1_ax.set_yticklabels(['60%', '70%', '80%', '90%', '100%'])
    fig1_ax.set_xlim(0.6, 1)
    fig1_ax.set_xticks(np.linspace(0.6, 1, 5))
    fig1_ax.set_xticklabels(['60%', '70%', '80%', '90%', '100%'])
    fig1_ax.plot(fig1_ax.get_xlim(), fig1_ax.get_ylim(), linestyle='--')

    for platform in platforms:
        caller_nsig_pcrt = [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
        for aligner_idx in range(len(aligners)):

            aligner = aligners[aligner_idx]

            lack_sig = []
            lack_sig_pcrt = []
            for caller_idx in range(len(long_read_callers)):
                caller = long_read_callers[caller_idx]

                sig_tag = {'Lack_sig_reads': 0, 'Enough_sig_reads': {'Lack_nsig_reads': 0, 'Enough_nsig_reads': 0}}

                total_count = 0
                for sample in samples:
                    annot_pavs = f'{workdir}/{sample}/{caller}/bed/{sample}.{platform}.{aligner}.{caller}.classified.missed.pavs.tsv'
                    df_annot_pavs = pd.read_csv(annot_pavs, sep='\t', names=['chrom', 'start', 'end', 'svtype', 'svlen', 'rptype', 'rppcrt', 'mapq', 'signum', 'nsignum', 'caller', 'platform', 'aligner'])

                    for idx, row in df_annot_pavs.iterrows():
                        sig_num = int(row['signum'])
                        nsig_num = int(row['nsignum'])
                        mapq = int(row['mapq'])
                        if mapq >= 20:
                            total_count += 1
                            if sig_num < 2:
                                sig_tag['Lack_sig_reads'] += 1
                            else:
                                if nsig_num < 2:
                                    sig_tag['Enough_sig_reads']['Lack_nsig_reads'] += 1
                                else:
                                    sig_tag['Enough_sig_reads']['Enough_nsig_reads'] += 1

                for tag, count in sig_tag.items():
                    if tag == 'Lack_sig_reads':
                        lack_sig.append(count / 3)
                        lack_sig_pcrt.append(count / total_count)

                    elif tag == 'Enough_sig_reads':
                        caller_nsig_pcrt[caller_idx][aligner_idx] += count['Enough_nsig_reads'] / total_count

        for i in range(len(caller_nsig_pcrt)):
            caller = long_read_callers[i]
            fig1_ax.scatter(caller_nsig_pcrt[i][0], caller_nsig_pcrt[i][1], facecolor=PLATCOLORS[platform], marker=TOOLMARKERS[caller], s=220, edgecolor='black')

    plt.tight_layout()
    plt.show()
    fig1.savefig(f'{workdir}/Figures/Figure5/samples.missed_pavs.nsigs.pdf')

'''

def plot_missed_pav_count(workdir, samples):
    long_read_callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']

    aligners = ['ngmlr', 'minimap2']
    platforms = ['hifi', 'ont']

    missed_pav_calls = []

    for platform in platforms:
        for aligner_idx in range(len(aligners)):
            aligner = aligners[aligner_idx]

            for caller_idx in range(len(long_read_callers)):
                caller = long_read_callers[caller_idx]

                for sample in samples:
                    annot_pavs = f'{workdir}/{sample}/{caller}/bed/{sample}.{platform}.{aligner}.{caller}.classified.missed.pavs.tsv'
                    df_annot_pavs = pd.read_csv(annot_pavs, sep='\t',
                                                names=['chrom', 'start', 'end', 'svtype', 'svlen', 'rptype', 'rppcrt',
                                                       'mapq', 'signum', 'nsignum', 'caller', 'platform', 'aligner'])

                    missed_pav_calls.append((len(df_annot_pavs), caller, platform, aligner))

    df_missed_pavs = pd.DataFrame(missed_pav_calls, columns=['count', 'caller', 'platform', 'aligner'])

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1,1,1)
    sns.lineplot(data=df_missed_pavs, x='caller', y='count', hue='aligner', style='platform', lw=2.5, ax=ax, markers=True, markersize=10, palette=['#1b9e77', '#7570b3'])
    plt.tight_layout()
    plt.show()

def plot_nsig_tag_in_high_mapq(workdir, sample):
    long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']
    # long_read_callers = ['sniffles']

    aligners = ['minimap2']
    platforms = ['hifi']

    plot_rptypes = ['VNTR', 'STR', 'LTR', 'LINE', 'SINE']

    tagged = []
    tagged_rptypes = []

    for platform in platforms:
        for aligner_idx in range(len(aligners)):
            aligner = aligners[aligner_idx]

            for caller_idx in range(len(long_read_callers)):
                caller = long_read_callers[caller_idx]

                annot_pavs = f'{workdir}/{sample}/{caller}/bed/{sample}.{platform}.{aligner}.{caller}.classified.missed.pavs.tsv'
                df_annot_pavs = pd.read_csv(annot_pavs, sep='\t',
                                            names=['chrom', 'start', 'end', 'svtype', 'svlen', 'rptype', 'rppcrt',
                                                   'mapq', 'signum', 'nsignum', 'caller', 'platform', 'aligner'])

                tag = {'Lack_nsig_reads': 0, 'Enough_nsig_reads': 0}

                tag_rptypes = {'Lack_nsig_reads': {'VNTR': 0, 'STR': 0, 'LTR': 0, 'LINE': 0, 'SINE': 0, 'Others': 0},
                               'Enough_nsig_reads': {'VNTR': 0, 'STR': 0, 'LTR': 0, 'LINE': 0, 'SINE': 0, 'Others': 0}}

                for idx, row in df_annot_pavs.iterrows():

                    nsig_num = int(row['nsignum'])
                    rptype = row['rptype']
                    if rptype not in plot_rptypes:
                        rptype = 'Others'

                    mapq = int(row['mapq'])

                    if mapq >= 50:
                        if nsig_num < 5:
                            tag['Lack_nsig_reads'] += 1
                            tag_rptypes['Lack_nsig_reads'][rptype] += 1
                        else:
                            tag['Enough_nsig_reads'] += 1
                            tag_rptypes['Enough_nsig_reads'][rptype] += 1

                for tag, count in tag.items():
                    tagged.append((tag, count / len(df_annot_pavs), caller, f'{platform}-{aligner}'))

                for tag, rp_count_dict in tag_rptypes.items():
                    for rptype, count in rp_count_dict.items():
                        tagged_rptypes.append((tag, rptype, count / len(df_annot_pavs), caller, f'{platform}-{aligner}'))
                        
def plot_mapq_tags(workdir, sample):

    long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']
    # long_read_callers = ['sniffles']

    aligners = ['minimap2']
    platforms = ['hifi']

    plot_rptypes = ['VNTR', 'STR', 'LTR', 'LINE', 'SINE']

    tagged = []
    tagged_rptypes = []
    for platform in platforms:
        for aligner_idx in range(len(aligners)):
            aligner = aligners[aligner_idx]

            for caller_idx in range(len(long_read_callers)):
                caller = long_read_callers[caller_idx]

                annot_pavs = f'{workdir}/{sample}/{caller}/bed/{sample}.{platform}.{aligner}.{caller}.classified.missed.pavs.tsv'
                df_annot_pavs = pd.read_csv(annot_pavs, sep='\t', names=['chrom', 'start', 'end', 'svtype', 'svlen', 'rptype', 'rppcrt', 'mapq', 'signum', 'nsignum', 'caller', 'platform', 'aligner'])


                tag = {'No_reads': 0, 'Low_mapq': 0, 'High_mapq': 0, 'Others': 0}

                tag_rptypes = {'No_reads': {'VNTR':0, 'STR':0, 'LTR':0, 'LINE':0, 'SINE':0, 'Others':0},
                            'Low_mapq': {'VNTR':0, 'STR':0, 'LTR':0, 'LINE':0, 'SINE':0, 'Others':0},
                            'High_mapq': {'VNTR':0, 'STR':0, 'LTR':0, 'LINE':0, 'SINE':0, 'Others':0},
                            'Others': {'VNTR':0, 'STR':0, 'LTR':0, 'LINE':0, 'SINE':0, 'Others':0}}

                for idx, row in df_annot_pavs.iterrows():
                    rptype = row['rptype']
                    if rptype not in plot_rptypes:
                        rptype = 'Others'

                    mapq = int(row['mapq'])
                    if mapq == -1:
                        tag['No_reads'] += 1
                        tag_rptypes['No_reads'][rptype] += 1
                    elif mapq < 20:
                        tag['Low_mapq'] += 1
                        tag_rptypes['Low_mapq'][rptype] += 1
                    elif mapq >= 50:
                        tag['High_mapq'] += 1
                        tag_rptypes['High_mapq'][rptype] += 1
                    else:
                        tag['Others'] += 1
                        tag_rptypes['Others'][rptype] += 1

                for tag, count in tag.items():
                    tagged.append((tag, count * 100 / len(df_annot_pavs), caller, f'{platform}-{aligner}'))

                for tag, rp_count_dict in tag_rptypes.items():
                    for rptype, count in rp_count_dict.items():
                        tagged_rptypes.append((tag, rptype, count * 100 / len(df_annot_pavs), caller, f'{platform}-{aligner}'))

    df_tagged = pd.DataFrame(tagged, columns=['mapq_tag', 'pcrt', 'caller', 'platform-aligner'])
    df_tagged_rptypes = pd.DataFrame(tagged_rptypes, columns=['mapq_tag', 'rptype', 'pcrt', 'caller', 'platform-aligner'])


    fig = plt.figure(figsize=(11, 4))
    ax1 = fig.add_subplot(1, 2, 1)
    sns.barplot(data=df_tagged, x='mapq_tag', y='pcrt', hue='caller', ci=None, ax=ax1)

    ax2 = fig.add_subplot(1, 2, 2)
    sns.barplot(data=df_tagged_rptypes, x='rptype', y='pcrt', hue='mapq_tag', ci=None, ax=ax2)

    plt.tight_layout()
    plt.show()

'''