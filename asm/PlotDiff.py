#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/3/3

'''
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import numpy as np
import seaborn as sns

from helpers.Reader import *


def platform_concordant_loci_supp(workdir, callers, aligner, figdir):



    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    supp_labels = [1, 2, 3, 4, 5, 6]
    xticks = np.arange(len(supp_labels))

    for caller_idx, caller in enumerate(callers):
        this_caller_pcrt = []

        platform_concordants = pd.read_csv(f'{workdir}/plat_repro/{aligner}/platform_concordant_loci.{caller}.annot.tsv', sep='\t',
                                           names=['chrom', 'start', 'id', 'svlen', 'svtype', 'supp', 'region', 'rptype', 'pcrt', 'aligner'])
        total_merged = len(platform_concordants)
        this_caller_supps = {i: 0 for i in supp_labels}
        for idx, row in platform_concordants.iterrows():
            supp = int(row['supp'])
            this_caller_supps[supp] += 1

        for supp, count in this_caller_supps.items():
            this_caller_pcrt.append(count * 100 / total_merged)

        ax.plot(xticks, this_caller_pcrt, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, markersize=9, label=TOOLMAP[caller])
    plt.tight_layout()
    plt.show()




def platform_concordants_insdel_bysupp(workdir, callers, aligner, figdir):

    svtypes = ['ins', 'del']
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))


    for col_idx, svtype in enumerate(svtypes):
        for caller in callers:

            hifi_ax = axes[0][col_idx]
            all_dataset_ax = axes[1][col_idx]
            hifi_ax.set_title(SVTYPEMAP[svtype], fontsize=13)

            for ax in [hifi_ax, all_dataset_ax]:
                if col_idx == 0:
                    ax.set_ylabel('% of SVs', fontsize=13)
                ax.set_ylim(0, 1)
                ax.set_yticks(np.linspace(0, 1, 5))
                ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)])

            hifi_merged_vcf = f'{workdir}/plat_repro/{aligner}/{caller}.{aligner}.{svtype}.hifi.merged.vcf'
            hifi_supp_dict, hifi_merged_total = get_survivor_supp(hifi_merged_vcf)

            hifi_data = []
            hifi_labels = [1, 2, 3, 4]
            for supp in hifi_labels:
                count = hifi_supp_dict[supp]
                hifi_data.append(count / hifi_merged_total)

            hifi_ax.plot(np.arange(len(hifi_labels)), hifi_data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, label=TOOLMAP[caller])
            hifi_ax.set_xticks(np.arange(len(hifi_labels)))
            hifi_ax.set_xticklabels([f'SUPP={ele}' for ele in hifi_labels], fontsize=12, rotation=90)
            hifi_ax.legend()

            merged_vcf = f'{workdir}/plat_repro/{aligner}/{caller}.{aligner}.{svtype}.merged.vcf'
            supp_dict, merged_total = get_survivor_supp(merged_vcf)

            data = []
            labels = [1, 2, 3, 4, 5, 6]
            for supp in labels:
                count = supp_dict[supp]
                data.append(count / merged_total)

            all_dataset_ax.plot(np.arange(len(labels)), data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller], lw=2, label=TOOLMAP[caller])
            all_dataset_ax.set_xticks(np.arange(len(labels)))
            all_dataset_ax.set_xticklabels([f'SUPP={ele}' for ele in labels], fontsize=12, rotation=90)
            all_dataset_ax.legend()

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{figdir}/platform_concordants_insdel_pcrt.pdf')

def plot_platform_unique_size(workdir, callers, datasets, aligner, figdir):

    svtypes = ['ins', 'del']

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(8, 4))

    for col_idx, svtype in enumerate(svtypes):

        datasets_svsize = []


        for caller in callers:
            merged_vcf = pysam.VariantFile(f'{workdir}/plat_repro/{aligner}/{caller}.{aligner}.{svtype}.merged.vcf')

            for rec in merged_vcf.fetch():
                supp, supp_vec = int(rec.info['SUPP']), rec.info['SUPP_VEC']

                if supp == 1:
                    svlen = abs(int(rec.info['SVLEN']))
                    datasets_svsize.append([svlen, caller])


        this_ax = axes[col_idx]
        df_svsize = pd.DataFrame(datasets_svsize, columns=['svlen', 'caller'])
        sns.histplot(df_svsize, x='svlen', hue='caller', log_scale=True, ax=this_ax)

    plt.tight_layout()
    plt.show()
