#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/2/23

'''
import pandas as pd
import pysam
import numpy as np
import vcf

from helpers.Annot import *
from helpers.Reader import *
from com.CheckFDR import *
from datasets.GetRepro import *



def compare_stra(workdir, datasets, aligners):
    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')
    regioned_svs_counts = []

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}



    for readset in datasets:
        plat = 'HiFi'
        if 'ont' in readset:
            plat = 'ONT'

        for assembler in plat_assemblers[plat]:
            for asm_method in ASMCALLERS:
                asm_calls = f'{workdir}/{plat}/minimap2_{readset}/filtered/HG002.{asm_method}.{assembler}.vcf'

                for caller in CALLERS:
                    for aligner in aligners:

                        print(f'Comparing {asm_method}-{assembler} to {caller}-{aligner} on {readset} ...')

                        caller_calls = f'{workdir}/{plat}/{aligner}_{readset}/filtered/HG002.{caller}.filtered.vcf'
                        compare_outdir = f'{workdir}/{plat}/{aligner}_{readset}/filtered/comstra'

                        # if 'ont' in readset:
                        #     caller_calls = f'{workdir}/ONT/{aligner}_{readset}/filtered/HG002.{caller}.filtered.vcf'
                        #     compare_outdir = f'{workdir}/ONT/{aligner}_{readset}/filtered/comstra'

                        if not os.path.exists(compare_outdir):
                            os.mkdir(compare_outdir)

                        tmp_file = f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.txt'
                        tmp_file_writer = open(tmp_file, 'w')
                        print(asm_calls, file=tmp_file_writer)
                        print(caller_calls, file=tmp_file_writer)

                        tmp_file_writer.close()

                        merged_vcf = f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.{readset}.jasmine.merged.vcf'
                        cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_vcf} max_dist=1000 spec_len=50 spec_reads=1'
                        os.system(cmd)

                        # os.remove(tmp_file)

                        # asm_uniques = open(f'{compare_outdir}/{caller}-{asm_method}.{aligner}-minimap2.{asm_method}.uniques.annot.tsv', 'w')
                        # read_uniques = open(f'{compare_outdir}/{caller}-{asm_method}.{aligner}-minimap2.{caller}.uniques.annot.tsv', 'w')
                        matched_info_out = open(f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.{readset}.concordants.info.tsv', 'w')
                        unique_info_out = open(f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.{readset}.uniques.info.tsv', 'w')

                        svs_by_regions = get_caller_compare_info([asm_method, caller], merged_vcf, matched_info_out, unique_info_out, simple_reps, rmsk, sds)

                        # svs_by_regions = {'Simple Repeats': [0, 0, 0], 'Repeat Masked': [0, 0, 0], 'Segment Dup': [0, 0, 0], 'Unique': [0, 0, 0]}

                        # vcf_reader = vcf.Reader(open(merged_vcf, 'r'))
                        # for rec in vcf_reader:
                        #     supp_vec = rec.INFO['SUPP_VEC']
                        #     chrom, start = rec.CHROM, int(rec.POS)
                        #     end, svlen = int(rec.INFO['END']), int(rec.INFO['SVLEN'])
                        #     svtype = rec.INFO['SVTYPE']
                        #
                        #     if start == end:
                        #         end += svlen
                        #
                        #     region_label, rptype, pcrt = annotate_sv_region(chrom, start, end, simple_reps, rmsk, sds)
                        #
                        #     if supp_vec == '10':
                        #         svs_by_regions[region_label][0] += 1
                        #         pav_id = rec.samples[0].data[7]
                                # print(f'{chrom}\t{start}\t{pav_id}\t{svlen}\t{svtype}\t{region_label}\t{rptype}\t{round(pcrt, 2)}\t{aligner}', file=asm_uniques)

                            # elif supp_vec == '11':
                            #     svs_by_regions[region_label][1] += 1
                            # else:
                            #     # print(f'{chrom}\t{start}\t{rec.ID}\t{svlen}\t{svtype}\t{region_label}\t{rptype}\t{round(pcrt, 2)}\t{aligner}', file=read_uniques)
                            #     svs_by_regions[region_label][2] += 1

                        # asm_uniques.close()
                        # read_uniques.close()

                        for region_label, counts in svs_by_regions.items():
                            regioned_svs_counts.append((caller, asm_method, readset, aligner, assembler, region_label, counts[0], counts[1], counts[2]))

    df_regioned_svs = pd.DataFrame(regioned_svs_counts, columns=['caller', 'asm_method', 'dataset', 'aligner', 'assembler', 'region', 'assm_unique','intersects', 'align_unique'])
    df_regioned_svs.to_csv(f'{workdir}/strategy_compare_byregions.tsv', header=True, sep='\t', index=False)


def compare_stra_at_regions(workdir, datasets, aligners):
    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    region = 'highconf'

    plat_assemblers = {'HiFi': ['hifiasm', 'flye'], 'ONT': ['flye', 'shasta']}
    regioned_svs_counts = []

    for readset in datasets:
        plat = 'HiFi'
        if 'ont' in readset:
            plat = 'ONT'

        for asm_method in ASMCALLERS:
            for assembler in plat_assemblers[plat]:
                asm_calls = f'{workdir}/{plat}/minimap2_{readset}/filtered_in_{region}/HG002.{asm_method}.{assembler}.vcf'
                for caller in CALLERS:
                    for aligner in aligners:

                        print(f'Comparing {asm_method} to {caller} on {aligner}-{readset} ...')

                        caller_calls = f'{workdir}/{plat}/{aligner}_{readset}/filtered_in_{region}/HG002.{caller}.vcf'
                        compare_outdir = f'{workdir}/{plat}/{aligner}_{readset}/filtered_in_{region}/comstra'


                        if not os.path.exists(compare_outdir):
                            os.mkdir(compare_outdir)

                        tmp_file = f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.txt'
                        tmp_file_writer = open(tmp_file, 'w')
                        print(asm_calls, file=tmp_file_writer)
                        print(caller_calls, file=tmp_file_writer)

                        tmp_file_writer.close()

                        merged_vcf = f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.{readset}.jasmine.merged.vcf'
                        cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_vcf} max_dist=1000 spec_len=50 spec_reads=1'
                        os.system(cmd)

                        os.remove(tmp_file)
                        matched_info_out = open(f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.{readset}.concordants.info.tsv', 'w')
                        unique_info_out = open(f'{compare_outdir}/{caller}-{aligner}.{asm_method}-minimap2-{assembler}.{readset}.uniques.info.tsv', 'w')

                        svs_by_regions = get_caller_compare_info([asm_method, caller], merged_vcf, matched_info_out, unique_info_out, simple_reps, rmsk, sds)

                        for region_label, counts in svs_by_regions.items():
                            regioned_svs_counts.append((caller, asm_method, readset, aligner, assembler, region_label, counts[0], counts[1], counts[2]))

    df_regioned_svs = pd.DataFrame(regioned_svs_counts, columns=['caller', 'asm_method', 'dataset', 'aligner', 'assember', 'region', 'assm_unique','intersects', 'align_unique'])
    df_regioned_svs.to_csv(f'{workdir}/strategy_compare_byregions_in_{region}.tsv', header=True, sep='\t', index=False)


def merge_uniques_at_regions(workdir, datasets, aligners):
    regions = ['cmrg', 'highconf']

    for region in regions:

        output_dir = f'{workdir}/{region}_merged_uniques'
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        asm_tmp_file = open(f'{output_dir}/asm_tmp.txt', 'w')
        read_tmp_file = open(f'{output_dir}/read_tmp.txt', 'w')

        for readset in datasets:
            plat = 'HiFi'
            if 'ont' in readset:
                plat = 'ONT'
            for asm_method in ASMCALLERS:
                for caller in CALLERS:
                    for aligner in aligners:
                        # print(f'Comparing {asm_method} to {caller} on {aligner}-{readset} ...')
                        compare_outdir = f'{workdir}/{plat}/{aligner}_{readset}/filtered_in_{region}/comstra'
                        merged_vcf = f'{compare_outdir}/{caller}-{asm_method}.{aligner}-minimap2.{readset}.jasmine.merged.vcf'
                        asm_unique_vcf = f'{output_dir}/{caller}-{asm_method}.{aligner}-minimap2.{readset}.{asm_method}.unique.vcf'
                        read_unique_vcf = f'{output_dir}/{caller}-{asm_method}.{aligner}-minimap2.{readset}.{caller}.unique.vcf'

                        separate_stra_merged_vcfs(merged_vcf, asm_unique_vcf, read_unique_vcf)

                        print(asm_unique_vcf, file=asm_tmp_file)
                        print(read_unique_vcf, file=read_tmp_file)


        asm_tmp_file.close()
        read_tmp_file.close()

        asm_merged_uniques = f'{output_dir}/assembly_based.uniques.jasmine.merged.vcf'
        read_merged_uniques = f'{output_dir}/read_based.uniques.jasmine.merged.vcf'

        cmd1 = f'{JASMINE} file_list={output_dir}/asm_tmp.txt out_file={asm_merged_uniques} max_dist=1000 spec_len=50 spec_reads=1'
        os.system(cmd1)
        cmd2 = f'{JASMINE} file_list={output_dir}/read_tmp.txt out_file={read_merged_uniques} max_dist=1000 spec_len=50 spec_reads=1'
        os.system(cmd2)




def compare_stra_insdel(workdir, datasets, aligners):

    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    svtypes = ['ins', 'del']

    for svtype in svtypes:
        regioned_svs_counts = []
        for readset in datasets:
            for asm_method in ASMCALLERS:
                asm_calls = f'{workdir}/HiFi/minimap2_{readset}/filtered/HG002.{asm_method}.{svtype}.vcf'
                if 'ont' in readset:
                    asm_calls = f'{workdir}/ONT/minimap2_{readset}/filtered/HG002.{asm_method}.{svtype}.vcf'

                for caller in CALLERS:
                    for aligner in aligners:

                        print(f'Comparing {asm_method} and {caller} on {aligner}-{readset} ...')

                        # read_calls = f'{workdir}/HiFi/{aligner}_{readset}/filtered/HG002.{caller}.{svtype}.vcf'
                        compare_outdir = f'{workdir}/HiFi/{aligner}_{readset}/filtered/comstra_insdel'
                        #
                        if 'ont' in readset:
                            read_calls = f'{workdir}/ONT/{aligner}_{readset}/filtered/HG002.{caller}.{svtype}.vcf'
                            compare_outdir = f'{workdir}/ONT/{aligner}_{readset}/filtered/comstra_insdel'
                        #
                        # if not os.path.exists(compare_outdir):
                        #     os.mkdir(compare_outdir)
                        #
                        # tmp_file = f'{compare_outdir}/{caller}-{asm_method}.{readset}.txt'
                        # tmp_file_writer = open(tmp_file, 'w')
                        # print(asm_calls, file=tmp_file_writer)
                        # print(read_calls, file=tmp_file_writer)
                        #
                        # tmp_file_writer.close()

                        merged_vcf = f'{compare_outdir}/{caller}-{asm_method}.{aligner}-minimap2.{readset}.{svtype}.jasmine.merged.vcf'
                        # cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_vcf} max_dist=1000 spec_len=50 spec_reads=1'

                        # os.system(cmd)
                        # os.remove(tmp_file)

                        # rep_annotate_uniques(compare_outdir, SVTYPEMAP[svtype], caller, merged_vcf, simple_reps, rmsk, sds)

                        # asm_uniques = open(f'{compare_outdir}/{caller}-{asm_method}.minimap2-{aligner}.{asm_method}.{svtype}.uniques.annot.tsv', 'w')
                        # align_uniques = open(f'{compare_outdir}/{caller}-{asm_method}.minimap2-{aligner}.{caller}.{svtype}.uniques.annot.tsv', 'w')

                        matched_info_out = open(f'{compare_outdir}/{caller}-{asm_method}.{aligner}-minimap2.{readset}.{svtype}.caller_concordant.info.tsv', 'w')
                        unique_info_out = open(f'{compare_outdir}/{caller}-{asm_method}.{aligner}-minimap2.{readset}.{svtype}.caller_unique.info.tsv','w')


                        svs_by_regions = get_caller_compare_info([asm_method, caller], merged_vcf, matched_info_out, unique_info_out, simple_reps, rmsk, sds)

                        # svs_by_regions = {'Simple Repeats': [0, 0, 0], 'Repeat Masked': [0, 0, 0], 'Segment Dup': [0, 0, 0], 'Unique': [0, 0, 0]}

                        # vcf_reader = vcf.Reader(open(merged_vcf, 'r'))
                        # for rec in vcf_reader:
                        #     supp_vec = rec.INFO['SUPP_VEC']
                        #     chrom, start = rec.CHROM, int(rec.POS)
                        #     end, svlen = int(rec.INFO['END']), int(rec.INFO['SVLEN'])
                        #
                        #     if start == end:
                        #         end += svlen
                        #
                        #     region_label, rptype, pcrt = annotate_sv_region(chrom, start, end, simple_reps, rmsk, sds)
                        #
                        #     if supp_vec == '10':
                        #         svs_by_regions[region_label][0] += 1
                        #         pav_id = rec.samples[0].data[7]
                        #         print(f'{chrom}\t{start}\t{pav_id}\t{svlen}\t{svtype}\t{region_label}\t{rptype}\t{round(pcrt, 2)}\t{aligner}', file=asm_uniques)
                        #
                        #     elif supp_vec == '11':
                        #         svs_by_regions[region_label][1] += 1
                        #
                        #     else:
                        #         print(f'{chrom}\t{start}\t{rec.ID}\t{svlen}\t{svtype}\t{region_label}\t{rptype}\t{round(pcrt, 2)}\t{aligner}', file=align_uniques)
                        #         svs_by_regions[region_label][2] += 1

                        for region_label, counts in svs_by_regions.items():
                            regioned_svs_counts.append((caller, asm_method, readset, aligner, region_label, counts[0], counts[1], counts[2]))

        df_regioned_svs = pd.DataFrame(regioned_svs_counts, columns=['caller', 'asm_method', 'dataset', 'aligner', 'region', 'assm_unique', 'intersects', 'align_unique'])
        df_regioned_svs.to_csv(f'{workdir}/strategy_compare_{svtype}_byregions.tsv', header=True, sep='\t', index=False)


# def annotate_assm_uniques(workdir, sample, datasets, aligner, bam_dir):
#
#     simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
#     rmsk = pysam.Tabixfile(iMACRMSK, 'r')
#     sds = pysam.Tabixfile(iMACSD, 'r')
#
#     long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']
#
#     bam_path = f'{bam_dir}/{sample}/{platform}/{sample}.{platform}.{aligner}.sorted.bam'
#
#     bam_file = pysam.AlignmentFile(bam_path)
#
#     for caller in long_read_callers:
#         missed_pavs = f'{workdir}/{caller}/bed/{sample}.{platform}.{aligner}.{caller}.bedtools.missed.pavs.bed'
#         classified_out = open(f'{workdir}/{sample}/{caller}/bed/{sample}.{platform}.{aligner}.{caller}.classified.missed.pavs.tsv', 'w')
#
#         df_missed_pavs = pd.read_csv(missed_pavs, sep='\t', usecols=[0,1,2,3,4], names=['chrom', 'start', 'end', 'svtype', 'svlen'])
#
#         print(f'Process {caller}, # missed PAV calls: {len(df_missed_pavs)}')
#
#         for idx, row in df_missed_pavs.iterrows():
#             chrom, start, end, svtype, svlen = row['chrom'], int(row['start']), int(row['end']), row['svtype'], abs(int(row['svlen']))
#
#             region, rptype, pcrt = annotate_sv_region(chrom, start, end, 0, simple_reps, rmsk, sds)
#             mean_mapq, num_signature_reads, nearest_sigs = check_sv_signature_in(chrom, start, end, bam_file, 50, 100000)
#
#             print(f'{chrom}\t{start}\t{end}\t{svtype}\t{svlen}\t{rptype}\t{pcrt}\t{mean_mapq}\t{num_signature_reads}\t{len(nearest_sigs)}\t{caller}\t{platform}\t{aligner}', file=classified_out)



'''

def unique_file_to_bed(workdir, read_datasets, aligners, asm_methods):
    svtypes = ['ins', 'del']

    for svtype in svtypes:
        for readset in read_datasets:
            for asm_method in asm_methods:

                for caller in CALLERS:
                    for aligner in aligners:
                        compare_outdir = f'{workdir}/{aligner}_{readset}/comasm_insdel'

                        # asm_uniques_bed = open(f'{compare_outdir}/{asm_method}_uniques.{asm_method}-{caller}.{svtype}.annot.bed', 'w')
                        # align_uniques_bed = open(f'{compare_outdir}/{caller}_uniques.{asm_method}-{caller}.{svtype}.annot.bed', 'w')
                        #
                        # with open(f'{compare_outdir}/{asm_method}_uniques.{asm_method}-{caller}.{svtype}.annot.tsv', 'r') as f:
                        #     for line in f:
                        #         entries = line.strip().split('\t')
                        #         chrom, start, svid, length = entries[0], int(entries[1]), entries[2], abs(int(entries[3]))
                        #         print(f'{chrom}\t{start}\t{start + length}\t{svid}', file=asm_uniques_bed)
                        #
                        # with open(f'{compare_outdir}/{caller}_uniques.{asm_method}-{caller}.{svtype}.annot.tsv', 'r') as f:
                        #     for line in f:
                        #         entries = line.strip().split('\t')
                        #         chrom, start, svid, length = entries[0], int(entries[1]), entries[2], abs(int(entries[3]))
                        #         print(f'{chrom}\t{start}\t{start + length}\t{svid}', file=align_uniques_bed)

                        # with open(f'{compare_outdir}/{caller}-{asm_method}.minimap2-{aligner}.{svtype}.matched.info.tsv', 'r') as f:
                        #     for line in f:
                        #         if '#' in line:
                        #             continue
                        #         entries = line.strip().split('\t')
                        #         chrom, start, svid, length = entries[0], int(entries[1]), entries[2], abs(int(entries[3]))
                        #         print(f'{chrom}\t{start}\t{start + length}\t{svid}', file=concordants_bed)
                        os.remove(f'{compare_outdir}/{caller}_uniques.{asm_method}-{caller}.{svtype}.annot.bed')
                        os.remove(f'{compare_outdir}/{asm_method}_uniques.{asm_method}-{caller}.{svtype}.annot.bed')
'''