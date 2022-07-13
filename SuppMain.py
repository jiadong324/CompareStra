#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/6/21

'''

from figures.SuppFigures import *

def main():


    max_size = 100000
    assm_aligners = ['minimap2']
    asm_methods = ['pav', 'svimasm']

    aligners = ['minimap2', 'ngmlr', 'lra', 'winnowmap']

    high_cov_datasets = ['hifi_15kb_20kb', 'minion', 'promethion']
    hifi_datasets = ['hifi_10kb', 'hifi_11kb', 'hifi_15kb', 'hifi_18kb']
    datasets = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb', 'ont_9kb', 'ont_19kb', 'ont_30kb']

    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'pav', 'svimasm']


    svnum_supps(iMACDIR, datasets, aligners)


if __name__ == '__main__':
    main()