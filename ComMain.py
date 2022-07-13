#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/2/17
'''

from datasets.ReproSupp import *
from datasets.BpRepro import *
from com.BasicCom import *

from figures.Figure1 import *
from figures.Figure2 import *
from figures.Figure3 import *
from figures.Figure4 import *
from figures.Figure5 import *

def main():


    max_size = 100000
    assm_aligners = ['minimap2']
    asm_methods = ['pav', 'svimasm']

    aligners = ['minimap2', 'ngmlr', 'lra', 'winnowmap']

    high_cov_datasets = ['hifi_15kb_20kb', 'minion', 'promethion']
    hifi_datasets = ['hifi_10kb', 'hifi_11kb', 'hifi_15kb', 'hifi_18kb']
    datasets = ['hifi_10kb', 'hifi_15kb', 'hifi_18kb', 'ont_9kb', 'ont_19kb', 'ont_30kb']

    callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision', 'pav', 'svimasm']

    '''
        Overview of detected SVs among datasets (Figure 1)
    '''
    ## Figure 1
    # overview_svnum(iMACDIR, datasets)

    # sv_regions_pieplot(iMACDIR)
    # svsize_distribution(iMACDIR, datasets, aligners)
    # read_based_svtypes(iMACDIR, datasets, aligners)

    ## Some help functions for Figure 1
    print_sv_num_stats(iMACDIR, datasets, aligners)
    # print_regions_stats()

    # sv_num_stats_of_regions(iMACDIR, datasets, aligners)

    ## Supplementary figures for Figure 1
    # size_distribution_at_region(iMACDIR, datasets)
    # insdel_pcrt_by_aligner(iMACDIR, datasets, aligners, callers, f'{Figure1}/insdel_pcrt_by_aligner.pdf')
    # sv_num_by_aligner(iMACDIR, datasets, aligners, callers, f'{Figure1}/wgs_sv_num.pdf')
    # plot_highconf_cmrg_sv_num(iMACDIR, datasets, aligners, callers)
    # insdel_pcrt_by_regions(iMACDIR, datasets, aligners, callers)


    '''
        Reproducibility of strategies among datasets (Figure 2)
    '''

    ## Figure 2

    # print_platform_unique_stats(iMACDIR, 'minimap2')

    # overview_dataset_repro(iMACDIR, 'minimap2')

    # plot_bpstd(iMACDIR, 'minimap2')
    # plot_bpstd_radar(iMACDIR, 'minimap2')
    # plot_bpstd_byregions(iMACDIR, 'minimap2')
    # plot_inacc_svs(iMACDIR, 'minimap2', 0)


    # overview_datasets_unique(iMACDIR, callers, 'minimap2')
    # datasets_unique_svsize(iMACDIR, 'minimap2')
    # read_dataset_unique_svtype(iMACDIR, 'minimap2')
    # assm_dataset_unique_svtype(iMACDIR)


    ## Supplementary figures for Figure 2

    # plot_left_right_bpstd(iMACDIR, ['minimap2'])
    # plot_del_bpstd(iMACDIR, ['minimap2'])
    # plot_ins_bpstd(iMACDIR, ['minimap2'])

    '''
        Reproducibility of strategies among datasets produced by one platform (Figure 3)
    '''

    # platform_repro_overview(iMACDIR)
    # platform_repro_bpstd(iMACDIR, 'minimap2')
    # platform_unique_region(iMACDIR, 'minimap2')


    '''
        Comparing between strategies (Figure 4)
    '''

    # print_stra_compare_stats(iMACDIR, aligners, datasets)

    # stra_concordant_overview(iMACDIR, aligners, datasets)
    # highconf_stra_compare_overview(iMACDIR, aligners, datasets)

    # stra_concordant_bp_overview(iMACDIR, aligners, datasets)
    # highconf_stra_concordant_bp_overview(iMACDIR, aligners, datasets)

    # stra_concordant_insdel_bpstd(iMACDIR, aligners, datasets)
    # highconf_stra_concordant_insdel_bpstd(iMACDIR, aligners, datasets)

    # inacc_concordant_svsize(iMACDIR, aligners, datasets)

    # stra_concordant_features(iMACDIR, aligners, datasets)
    # stra_concordant_bp_scatter(iMACDIR, aligners, datasets)

    # wgs_highconf_acc_inacc_overview(iMACDIR, aligners, datasets)
    # wgs_highconf_acc_insdel(iMACDIR, aligners, datasets)
    # wgs_highconf_acc_concordant_factors(iMACDIR, aligners, datasets)


    # plot_stra_concordant_features_at_regions(iMACDIR, aligners, datasets)


    '''
        Strategy uniques (Figure 5)
    '''

    # stra_unique_overview(iMACDIR, aligners, datasets, callers)
    # plot_merged_stra_uniques(['hifi_18kb', 'ont_30kb'])

    # overview_unique_annot(iMACDIR, ['hifi_18kb', 'ont_30kb'])

    # plot_assm_annot_results(iMACDIR, ['hifi_18kb', 'ont_30kb'])
    # plot_read_annot_results(iMACDIR, ['hifi_18kb', 'ont_30kb'])

    '''
        Best practice (Figure xxx)
    '''
    # overview_dataset_repro_highconf(iMACDIR, 'minimap2')
    # dataset_repro_at_regions(iMACDIR, 'minimap2')

    '''
        Analyze and plot strategy concordant calls
    '''

    # plot_diff_ins_types(iMACDIR, low_cov_datasets, aligners, iMACCOMFIGURE)
    # plot_insdel_breakpoints_shift(iMACDIR, asm_methods, aligners, hifi_datasets, iMACCOMFIGURE)

    '''
        Analyze and plot strategy unique calls
    '''

    # plot_regioned_unique_loci(iMACDIR, low_cov_datasets, asm_methods, iMACCOMFIGURE)
    # plot_regioned_unique_ins(iMACDIR, low_cov_datasets, asm_methods, iMACCOMFIGURE)
    # plot_regioned_unique_del(iMACDIR, low_cov_datasets, asm_methods, iMACCOMFIGURE)

    # plot_asm_unique_insdel_count(iMACDIR, low_cov_datasets, asm_methods, iMACCOMFIGURE)
    # plot_asm_unique_insdel_size(iMACDIR, low_cov_datasets, aligners, iMACCOMFIGURE)
    # plot_asm_unique_insdel_repeats(iMACDIR, low_cov_datasets, aligners, iMACCOMFIGURE)

    # plot_align_uniques_insdel_size(iMACDIR, low_cov_datasets, aligners, iMACCOMFIGURE)
    # plot_align_unique_insdel_repeats(iMACDIR, low_cov_datasets, aligners, iMACCOMFIGURE)


    '''
        Comparing merged callers to assembly calls
    '''

    # caller_accumulate_merge(iMACDIR, all_datasets)
    # accumlate_merged_to_asm(iMACDIR, all_datasets, 'giab_ctg')
    # recall_precision_of_acc_merged(iMACDIR, low_cov_datasets, high_cov_datasets, iMACCOMFIGURE)

    # get_merged_callers(iMACDIR, all_datasets)
    # compare_merged_insdel_to_asm(iMACDIR, all_datasets, 'giab_ctg', iMACBAM, iMACSIMPLE, iMACRMSK)
    # recall_precision_of_merged(iMACDIR, low_cov_datasets, high_cov_datasets, iMACCOMFIGURE)



    '''
       This step runs script pav_benchmark.sh, producing:
       1. Ouptut recall and precision measure of each caller under /pav_v112/platform_aligner_s5_benchmark_results.txt
       2. For each caller, generating the missed PAV calls under caller/bed/ directory
       '''

    # bedtools_overlap(samples, platforms, aligners)

    # scatter_precision_recall(iMAC, 'HG00733', platforms, iMACHGSVCFIGURE)
    # lineplot_precision_recall(iMAC, 'HG00733', platforms, iMACHGSVCFIGURE)

    # missed_pav_caller_unique(iMAC, samples)
    # plot_missed_pav_feature(MAC, samples)
    # plot_sig_tags_in_high_mapq(iMAC, samples)


    # print('==== Step4 Evaluating caller breakpoint accuracy ======')

    # print('*** 1. Callers breakpoint accuracy comparing with PAV calls')
    # find_nearest_matched_pavs(iMAC, platforms, samples)
    # samples_pav_bpshift_stackplot(MAC, samples)
    # samples_pav_bpshift_stackplot2(MAC, samples, aligners, platforms)

    # simple_reps = '/Users/jiadonglin/Data/genome/repeat_annot/simplerepeat.bed.gz'
    # rmsk = '/Users/jiadonglin/Data/genome/repeat_annot/rmsk.bed.gz'
    # samples_pav_inaccbp(MAC, samples, simple_reps, rmsk)

    # for sample in samples:
    #     pav_bpshift_of_sample(iMAC, sample)

    # for aligner in alignerts:
    #     pieplot_pav_bpshift_count_aligner(iMAC, ['HG00733'], aligner)

    # print('*** 2. Callers breakpoint accuracy comparing with Illumina calls')
    # find_match_ngs(iMAC, platforms, samples)

    # sample_ngs_bpshift_stackplot(MAC, samples, platforms)
    # sample_ngs_bpshift_stackplot2(MAC, samples, platforms)

    # for sample in samples:
    #     samples_ngs_bpshift_of_platform(iMAC, sample, platforms)

    # print('*** 3. Callers breakpoint accuracy of reccurent SVs among samples')

    # merge_samples(iMAC, f'{iMAC}/sample_merged', samples)
    # reccurent_sv_bpdiff(f'{MAC}/sample_merged')

if __name__ == '__main__':
    main()