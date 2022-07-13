#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/4/1

'''
import math

import matplotlib.pyplot as plt
import pandas as pd
import pysam
import numpy as np
import seaborn as sns
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


def plot_platrepro_suppvec(workdir, asm_callers, aligner, flag, figdir):


    fig, this_ax = plt.subplots(1, 1, figsize=(6, 4))
    labels = [1, 2, 3, 4, 5, 6]

    this_ax.set_title(aligner, fontsize=13)


    this_ax.set_ylabel('% of SVs', fontsize=13)

    this_ax.set_ylim(0, 1)
    this_ax.set_yticks(np.linspace(0, 1, 5))
    this_ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)])

    if aligner == 'minimap2':
        for caller in asm_callers:
            merged_vcf = f'{workdir}/assm_dataset_repro/{caller}.{aligner}.jasmine.merged.vcf'
            if flag == 'hifi':
                merged_vcf = f'{workdir}/assm_dataset_repro/hifi/{caller}.{aligner}.jasmine.merged.vcf'
            supp_dict, merged_total = get_survivor_supp(merged_vcf)

            data = []
            for supp in labels:
                count = supp_dict[supp]
                data.append(count / merged_total)

            this_ax.plot(np.arange(len(labels)), data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller],
                         lw=2, label=TOOLMAP[caller])
            this_ax.set_xticks(np.arange(len(labels)))
            this_ax.set_xticklabels([f'SUPP={ele}' for ele in labels], fontsize=12, rotation=90)
            this_ax.legend()

    for caller in CALLERS:
        merged_vcf = f'{workdir}/read_dataset_repro/{caller}.{aligner}.jasmine.merged.vcf'
        if flag == 'hifi':
            merged_vcf = f'{workdir}/read_dataset_repro/hifi/{caller}.{aligner}.jasmine.merged.vcf'

        supp_dict, merged_total = get_survivor_supp(merged_vcf)

        data = []
        for supp in labels:
            count = 0
            if supp in supp_dict:
                count = supp_dict[supp]
            data.append(count / merged_total)

        this_ax.plot(np.arange(len(labels)), data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller],
                     lw=2, label=TOOLMAP[caller])
        this_ax.set_xticks(np.arange(len(labels)))
        this_ax.set_xticklabels([f'SUPP={ele}' for ele in labels], fontsize=12, rotation=90)
        this_ax.legend()

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{figdir}/suppvec_pcrt_{flag}_datasets.pdf')

def plot_plat_region_repro_suppvec(workdir, asm_callers, aligners, region):

    fig, axes = plt.subplots(1, len(aligners), sharex='col', sharey='row', figsize=(12, 4))
    labels = [1, 2, 3, 4, 5, 6, 7]


    for col_idx, aligner in enumerate(aligners):
        this_ax = axes[col_idx]
        this_ax.set_title(aligner, fontsize=13)

        if col_idx == 0:
            this_ax.set_ylabel('% of SVs', fontsize=13)

        this_ax.set_ylim(0, 1)
        this_ax.set_yticks(np.linspace(0, 1, 5))
        this_ax.set_yticklabels([f'{int(val * 100)}%' for val in np.linspace(0, 1, 5)])

        if aligner == 'minimap2':
            for caller in asm_callers:
                merged_vcf = f'{workdir}/assm_dataset_repro/{region}_regions/{caller}.{aligner}.jasmine.merged.vcf'

                supp_dict, merged_total = get_survivor_supp(merged_vcf)

                data = []
                for supp in labels:
                    count = supp_dict[supp]
                    data.append(count / merged_total)

                this_ax.plot(np.arange(len(labels)), data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller],
                             lw=2, label=TOOLMAP[caller])
                this_ax.set_xticks(np.arange(len(labels)))
                this_ax.set_xticklabels([f'SUPP={ele}' for ele in labels], fontsize=12, rotation=90)
                this_ax.legend()

        for caller in CALLERS:
            merged_vcf = f'{workdir}/align_dataset_repro/{region}_regions/{caller}.{aligner}.jasmine.merged.vcf'

            supp_dict, merged_total = get_survivor_supp(merged_vcf)

            data = []
            for supp in labels:
                count = 0
                if supp in supp_dict:
                    count = supp_dict[supp]
                data.append(count / merged_total)

            this_ax.plot(np.arange(len(labels)), data, color=TOOLCOLORS[caller], marker=TOOLMARKERS[caller],
                         lw=2, label=TOOLMAP[caller])
            this_ax.set_xticks(np.arange(len(labels)))
            this_ax.set_xticklabels([f'SUPP={ele}' for ele in labels], fontsize=12, rotation=90)
            this_ax.legend()

    plt.tight_layout()
    plt.show()
    fig.savefig(f'{iMACCOMFIGURE}/suppvec_pcrt_{region}_datasets.pdf')
