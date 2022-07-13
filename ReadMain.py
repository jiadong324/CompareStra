#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/1/25
'''


from helpers.Constant import *
from read.BasicRead import *
from read.AlignerRepro import *
from read.DataRepro import *
from figures.Figure6 import *

def main():

    max_size = 100000

    aligners = ['minimap2', 'ngmlr', 'lra', 'winnowmap']
    hifi_datasets = ['hifi_10kb', 'hifi_11kb', 'hifi_15kb', 'hifi_18kb']
    ont_datasets = ['minion_27X', 'promethion_27X']
    datasets = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb', 'ont_9kb', 'ont_19kb', 'ont_30kb']

    '''
        Thesis figure: Aligner concordant calls by SV types (Figure 6)
    '''
    # compare_between_aligner(iMACDIR, aligners, datasets)

    # print_aligner_unionset_stats(iMACDIR, datasets, aligners)

    # overview_aligner_repro(iMACDIR, datasets)
    # aligner_concordant_bpstd(iMACDIR, datasets, aligners)

    # aligner_unique_features(iMACDIR, datasets, aligners)
    # aligner_unique_svtypes(iMACDIR, datasets, aligners)

    aligner_unique_subset_svtypes(iMACDIR, datasets, aligners)

    # plot_truvari_results()

    '''
        Paper figure
    '''
    ## Supp Fig: Breakpoint std. of dataset concordant calls of aligners
    # plot_bpstd(iMACDIR, aligners)

    ## Supp Fig: Read-based dataset concordant rate of aligners
    # overview_dataset_repro(iMACDIR, aligners)

    ## Supp Fig: Read-based inaccurate breakpoint dataset concordant calls of aligners
    # plot_inacc_svs(iMACDIR, aligners, 0)



if __name__ == '__main__':
    main()