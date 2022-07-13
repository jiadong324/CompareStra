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

from helpers.Annot import *
from helpers.Functions import *
from com.CheckFDR import *
from helpers.Reader import *


def plot_diff_ins_types(workdir, datasets, aligners, figdir):

    df_disc_ins = pd.read_csv(f'{workdir}/asm_align_disc_ins.tsv', sep='\t', header=0)

    callers_disc_ins = {caller: [[0 for i in range(len(datasets))] for i in range(len(aligners))] for caller in CALLERS}

    for idx, row in df_disc_ins.iterrows():
        dataset, aln_caller, aligner, asm_caller, svtype, pcrt = row['dataset'], row['aln_caller'], row['aligner'], row['asm_caller'], row['svtype'], float(row['pcrt'])
        if svtype == 'DUP':
            aligner_idx = aligners.index(aligner)
            dataset_idx = datasets.index(dataset)
            callers_disc_ins[aln_caller][aligner_idx][dataset_idx] = pcrt

    fig, axes = plt.subplots(1, 4, sharex='col', sharey='row', figsize=(12, 4))

    xticks = np.arange(len(datasets))
    for i, caller in enumerate(CALLERS):
        ax = axes[i]
        for j, aligner in enumerate(aligners):
            ax.plot(xticks, callers_disc_ins[caller][j], label=aligner, color=ALIGNERCOLOR[aligner], marker='X', lw=2, markersize=9)

        ax.set_title(TOOLMAP[caller], fontsize=13)
        ax.set_xticks(xticks)
        ax.set_xticklabels([PLATMAP[plat] for plat in datasets], fontsize=13, rotation=90)
        ax.legend()

        if i == 0:
            ax.set_ylabel('% of DUP', fontsize=13)
        ax.set_ylim(0, 40)
        ax.set_yticks(np.linspace(0, 40, 5))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 40, 5)], fontsize=12)

    plt.tight_layout()
    # plt.show()
    fig.savefig(f'{figdir}/pcrt_disc_ins_type.pdf')
