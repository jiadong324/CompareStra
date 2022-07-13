#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/3/28

'''

from asm.AssmReader import process_assm_calls, annotate_assm_highconf_svs
from read.AlignReader import process_read_calls, annotate_read_highconf_svs
from com.BasicCom import *
from datasets.GetRepro import *
from com.ComStra import *
from com.CheckFDR import *


def create_bins():

    window_size = 10000000

    chrom_size_dict = {}
    chrom_binned_trees = {}
    chrom_binned_list = []
    most_bins = 0
    with open(CHROMSIZE, 'r') as f:
        for line in f:
            entries = line.strip().split('\t')
            chrom, size = entries[0], int(entries[1])
            if chrom in VALID_CHROMS:
                chrom_size_dict[chrom] = size

                chrom_binned_trees[chrom] = IntervalTree()
                bins = int(size / window_size)
                start = 0
                for ith in range(bins):
                    end = start + window_size
                    chrom_binned_trees[chrom][start: end] = (ith, start)
                    start = end

                if start < size:
                    chrom_binned_trees[chrom][start: size] = (bins, start)

                if chrom == 'chr1':
                    most_bins = bins + 1
                chrom_binned_list.append([0 for i in range(most_bins)])

    return chrom_binned_trees, chrom_binned_list

def main():

    max_size = 100000
    datasets = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb', 'ont_9kb', 'ont_19kb', 'ont_30kb']
    aligners = ['minimap2', 'ngmlr', 'lra', 'winnowmap']
    hifi_datasets = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb']
    ont_datasets = ['ont_9kb', 'ont_19kb', 'ont_30kb']

    ont_assemblers = ['shasta', 'flye']
    hifi_assemblers = ['hifiasm', 'flye']

    '''
        Process raw calls of each caller
    '''

    # process_read_calls(iMACDIR, datasets, aligners, max_size)

    # assm_binned_tree, assm_binned_list = create_bins()

    # process_assm_calls(iMACDIR, max_size)

    # annotate_read_highconf_svs(iMACDIR, datasets, aligners)
    # annotate_assm_highconf_svs(iMACDIR, ont_datasets, ont_assemblers)
    # annotate_assm_highconf_svs(iMACDIR, hifi_datasets, hifi_assemblers)


    '''
        Examine caller reproducibiltiy among datasets
    '''
    ## Reproducibility at whole genome scale
    # merge_datasets_read_calls(iMACDIR, aligners, datasets)
    # merge_read_calls_by_platform(iMACDIR, aligners, {'hifi': hifi_datasets, 'ont': ont_datasets})

    # merge_datasets_assm_calls(iMACDIR, hifi_assemblers, ont_assemblers, datasets)

    # merge_assm_calls_by_platform(iMACDIR, hifi_assemblers, hifi_datasets)
    # merge_assm_calls_by_platform(iMACDIR, ont_assemblers, ont_datasets)

    ## Reproducibility at high confident regions and CMRG
    # merge_datasets_read_at_regions(iMACDIR, aligners, datasets, {'hifi': hifi_datasets, 'ont': ont_datasets})
    # merge_datasets_assm_at_regions(iMACDIR, ['minimap2'], datasets, {'hifi': hifi_datasets, 'ont': ont_datasets})


    '''
        Compare strategy for each dataset
    '''
    ## Strategy comparison at whole genome scale

    # compare_stra(iMACDIR, ont_datasets, aligners)
    # compare_stra_at_regions(iMACDIR, datasets, aligners)

    # merge_assm_uniques(iMACDIR, ['hifi_18kb', 'ont_30kb'])
    annotate_uniques(iMACDIR, 'hifi_18kb', f'{iMACBAM}/hifi_18kb/HG002.hifi.minimap2.sorted.bam', 500)
    annotate_uniques(iMACDIR, 'ont_30kb', f'{iMACBAM}/ont_30kb/HG002.ont.minimap2.sorted.bam', 500)

    # merge_uniques_at_regions(iMACDIR, ['hifi_18kb'], ['minimap2'])

if __name__ == '__main__':
    main()