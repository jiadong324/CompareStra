#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/3/16

'''
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

from helpers.Constant import *

def plot_hifi_ont_uniques(workdir, callers, figdir):

    region_labels = ['Simple Repeats', 'Repeat Masked', 'Segment Dup', 'Unique']

    column_names = ['chrom', 'start', 'id', 'size', 'svtype', 'supp', 'region', 'rptype', 'pcrt', 'aligner']

    sv_info_list = []

    ont_sv_region_pcrt = {caller: [0 for i in range(len(region_labels))] for caller in callers}
    hifi_sv_region_pcrt = {caller: [0 for i in range(len(region_labels))] for caller in callers}

    for caller in callers:
        this_caller_sv_info = []
        ont_svtype_dict = {}
        ont_sv_region_dict = {}

        hifi_svtype_dict = {}
        hifi_sv_region_dict = {}

        df_ont_annot = pd.read_csv(f'{workdir}/plat_repro/minimap2/ont_uniques_loci.{caller}.annot.tsv', sep='\t', names=column_names)
        df_hifi_annot = pd.read_csv(f'{workdir}/plat_repro/minimap2/hifi_uniques_loci.{caller}.annot.tsv', sep='\t', names=column_names)

        for idx, row in df_ont_annot.iterrows():
            svsize, region, svtype = abs(int(row['size'])), row['region'], row['svtype']
            if region in ont_sv_region_dict:
                ont_sv_region_dict[region] += 1
            else:
                ont_sv_region_dict[region] = 1

            this_caller_sv_info.append((caller, svsize, svtype, region, 'ONT'))


        for region, num in ont_sv_region_dict.items():
            region_idx = region_labels.index(region)
            # ont_sv_region_pcrt[caller][region_idx] = num * 100 / len(df_ont_annot)
            ont_sv_region_pcrt[caller][region_idx] = num

        for idx, row in df_hifi_annot.iterrows():
            svsize, region, svtype = abs(int(row['size'])), row['region'], row['svtype']
            if region in hifi_sv_region_dict:
                hifi_sv_region_dict[region] += 1
            else:
                hifi_sv_region_dict[region] = 1

            this_caller_sv_info.append((caller, svsize, svtype, region, 'HiFi'))


        for region, num in hifi_sv_region_dict.items():
            region_idx = region_labels.index(region)
            # hifi_sv_region_pcrt[caller][region_idx] = num * 100 / len(df_hifi_annot)
            hifi_sv_region_pcrt[caller][region_idx] = num


        sv_info_list.append(this_caller_sv_info)


    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    xticks = np.arange(len(region_labels))

    hifi_ax = axes[0]
    hifi_ax.set_title('HiFi', fontsize=13)
    for caller, pcrt in hifi_sv_region_pcrt.items():
        hifi_ax.plot(xticks, pcrt, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])
    hifi_ax.legend()

    hifi_ax.set_xticks(xticks)
    hifi_ax.set_xticklabels(region_labels, rotation=90)
    hifi_ax.set_ylim(0, 400)
    hifi_ax.set_yticks(np.linspace(0, 400, 5))
    hifi_ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 4, 5)], fontsize=13)

    hifi_ax.set_ylabel('# of SVs (x100)')

    ont_ax = axes[1]
    ont_ax.set_title('ONT', fontsize=13)
    for caller, pcrt in ont_sv_region_pcrt.items():
        ont_ax.plot(xticks, pcrt, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])
    ont_ax.legend()

    ont_ax.set_xticks(xticks)
    ont_ax.set_xticklabels(region_labels, rotation=90)

    ont_ax.set_ylim(0, 10000)
    ont_ax.set_yticks(np.linspace(0, 10000, 5))
    ont_ax.set_yticklabels([f'{int(val)}' for val in np.linspace(0, 100, 5)], fontsize=13)


    fig1, axes = plt.subplots(1, 2, sharex='col', sharey='row',figsize=(8, 4))

    for i, caller in enumerate(callers):
        ax = axes[i]
        ax.set_title(TOOLMAP[caller], fontsize=13)
        df_data = pd.DataFrame(sv_info_list[i], columns=['caller', 'size', 'svtype', 'region', 'Platform'])
        sns.histplot(data=df_data, x='size', hue='Platform', log_scale=[True, True], bins=50, ax=ax)
        ax.set_xlabel('')
        if i == 0:
            ax.set_ylabel('# of SVs', fontsize=12)
            # ax.set_ylim(0, 3000)
            # ax.set_yticks(np.linspace(0, 3000, 5))
            # ax.set_yticklabels([int(val) for val in np.linspace(0, 3000, 5)])


    fig.tight_layout()
    fig1.tight_layout()
    plt.show()

    fig.savefig(f'{figdir}/hifi_ont_uniques_by_regions.pdf')
    fig1.savefig(f'{figdir}/hifi_ont_unique_svsize.pdf')