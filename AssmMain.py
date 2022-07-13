#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/2/22

'''

from asm.PlatDiff import *
from datasets.GetRepro import *

def main():

    max_size = 100000
    callers = ['pav', 'svimasm']
    datasets = ['hifi_10kb', 'hifi_11kb', 'hifi_15kb', 'hifi_15kb_20kb_27X', 'hifi_18kb_27X', 'minion_27X', 'promethion_27X']
    hifi_datasets = ['hifi_10kb', 'hifi_11kb', 'hifi_15kb', 'hifi_15kb_20kb_27X']

    aligners = ['minimap2']



    '''
        Plot basic statistics
    '''

    # plot_sv_num(iMACDIR, hifi_datasets)
    # plot_sv_by_regions(iMACDIR, datasets, iMACASMFIGURE)
    # plot_insdel(iMACDIR, callers, datasets)

    '''
        Concordant and discordant INS/DEL of each caller among datasets
    '''

    # platform_concordant_insdel(iMACDIR, callers, datasets, 'minimap2')
    # hifi_platform_concordant_insdel(iMACDIR, callers, hifi_datasets, 'minimap2')

    # platform_concordant_insdel_bpshift(iMACDIR, callers, 'minimap2', iMACASMFIGURE)


    '''
        Platform concordant and unique loci
    '''

    # platform_concordant_loci(iMACDIR, callers, datasets, 'minimap2')
    # plot_hifi_ont_uniques(iMACDIR, callers, iMACASMFIGURE)

    # platform_concordant_loci_supp(iMACDIR, callers, 'minimap2', iMACASMFIGURE)


    '''
        Analyze concordant INS/DEL among datasets
    '''

    # platform_concordants_insdel_bysupp(iMACDIR, callers, 'minimap2', iMACASMFIGURE)

    '''
        Analyze dataset unique INS/DEL        
    '''

    # plot_platform_unique_size(iMACDIR, callers, datasets, 'minimap2', iMACASMFIGURE)


    '''
        Find caller concordant and discordant calls
    '''

    # caller_concordant_insdel(iMACDIR, callers, datasets, 'minimap2', iMACSIMPLE, iMACRMSK, iMACSD)


if __name__ == '__main__':
    main()