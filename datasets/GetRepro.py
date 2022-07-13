#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/3/31

'''
import os
import pysam
import math
import pandas as pd

from helpers.Annot import *
from helpers.Constant import *


def merge_datasets_read_calls(workdir, aligners, datasets):
    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    for caller in CALLERS:
        for aligner in aligners:
            ## Matches of all SVs between platforms
            merged_outdir = f'{workdir}/read_dataset_repro'
            if not os.path.exists(merged_outdir):
                os.mkdir(merged_outdir)

            tmp_file = f'{merged_outdir}/{caller}.{aligner}.txt'
            tmp_file_writer = open(tmp_file, 'w')

            for dataset in datasets:
                vcf_file = f'{workdir}/HiFi/{aligner}_{dataset}/filtered/HG002.{caller}.filtered.vcf'
                if 'ont' in dataset:
                    vcf_file = f'{workdir}/ONT/{aligner}_{dataset}/filtered/HG002.{caller}.filtered.vcf'
                print(f'{vcf_file}', file=tmp_file_writer)

            tmp_file_writer.close()

            merged_out_vcf = f'{merged_outdir}/{caller}.{aligner}.jasmine.merged.vcf'
            print(f'Producing {caller} {aligner} merged calls ...')
            cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'
            os.system(cmd)

            os.remove(tmp_file)

            match_info_out = f'{merged_outdir}/{caller}.{aligner}.dastsets-concordant.info.tsv'
            unique_info_out = f'{merged_outdir}/{caller}.{aligner}.dastsets.unique.tsv'
            get_dataset_compare_info(datasets, merged_out_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds)

            hifi_uniques = f'{merged_outdir}/{caller}.{aligner}.hifi.uniques.info.tsv'
            ont_uniques = f'{merged_outdir}/{caller}.{aligner}.ont.uniques.info.tsv'
            get_platform_uniques(merged_out_vcf, hifi_uniques, ont_uniques, simple_reps, rmsk, sds)

def merge_read_calls_by_platform(workdir, aligners, datasets_dict):
    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')


    for platform, datasets in datasets_dict.items():
        for caller in ['svision']:
            for aligner in aligners:
                ## Matches of all SVs between platforms
                merged_outdir = f'{workdir}/read_dataset_repro/{platform}'
                if not os.path.exists(merged_outdir):
                    os.mkdir(merged_outdir)

                tmp_file = f'{merged_outdir}/{caller}.{aligner}.txt'
                tmp_file_writer = open(tmp_file, 'w')

                for dataset in datasets:

                    vcf_file = f'{workdir}/{PLATMAP[platform]}/{aligner}_{dataset}/filtered/HG002.{caller}.filtered.vcf'
                    print(f'{vcf_file}', file=tmp_file_writer)
                tmp_file_writer.close()

                merged_out_vcf = f'{merged_outdir}/{caller}.{aligner}.jasmine.merged.vcf'
                print(f'Producing {caller} {aligner} merged calls ...')
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'
                os.system(cmd)

                os.remove(tmp_file)

                match_info_out = f'{merged_outdir}/{caller}.{aligner}.{platform}-concordant.info.tsv'
                unique_info_out = f'{merged_outdir}/{caller}.{aligner}.{platform}.unique.tsv'
                get_dataset_compare_info(datasets, merged_out_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds)

def merge_assm_calls_by_platform(workdir, assemblers, datasets):
    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    aligner = 'minimap2'

    asm_callers = ['pav', 'svimasm']
    for dataset in datasets:

        platform = 'HiFi'
        if 'ont' in dataset:
            platform = 'ONT'

        for caller in asm_callers:
            for assembler in assemblers:
                ## Matches of all SVs between platforms
                merged_outdir = f'{workdir}/assm_dataset_repro/{platform}'
                if not os.path.exists(merged_outdir):
                    os.mkdir(merged_outdir)

                tmp_file = f'{merged_outdir}/{caller}.{assembler}.txt'
                tmp_file_writer = open(tmp_file, 'w')

                for dataset in datasets:
                    vcf_file = f'{workdir}/{platform}/{aligner}_{dataset}/filtered/HG002.{caller}.{assembler}.vcf'
                    print(f'{vcf_file}', file=tmp_file_writer)
                tmp_file_writer.close()

                merged_out_vcf = f'{merged_outdir}/{caller}.{assembler}.jasmine.merged.vcf'
                print(f'Producing {caller} {assembler} merged calls ...')
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'
                os.system(cmd)

                os.remove(tmp_file)

                match_info_out = f'{merged_outdir}/{caller}.{assembler}.{platform}-concordant.info.tsv'
                unique_info_out = f'{merged_outdir}/{caller}.{assembler}.{platform}.unique.tsv'
                get_dataset_compare_info(datasets, merged_out_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds)

def merge_datasets_assm_calls(workdir, hifi_assemblers, ont_assemblers, datasets):
    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    aligner = 'minimap2'

    assembler_orders = [[hifi_assemblers[0], ont_assemblers[0]],
                        [hifi_assemblers[0], ont_assemblers[1]],
                        [hifi_assemblers[1], ont_assemblers[1]],
                        [hifi_assemblers[1], ont_assemblers[0]]]

    for caller in ['pav', 'svimasm']:
        for (hifi_assembler, ont_assembler) in assembler_orders:
            print(f'Producing {caller} {hifi_assembler}-{ont_assembler} merged calls ...')
            ## Matches of all SVs between platforms
            merged_outdir = f'{workdir}/assm_dataset_repro'
            if not os.path.exists(merged_outdir):
                os.mkdir(merged_outdir)

            tmp_file = f'{merged_outdir}/{caller}.{aligner}.txt'
            tmp_file_writer = open(tmp_file, 'w')

            for dataset in datasets:
                if 'hifi' in dataset:
                    vcf_file = f'{workdir}/HiFi/{aligner}_{dataset}/filtered/HG002.{caller}.{hifi_assembler}.vcf'
                    print(f'{vcf_file}', file=tmp_file_writer)
                if 'ont' in dataset:
                    vcf_file = f'{workdir}/ONT/{aligner}_{dataset}/filtered/HG002.{caller}.{ont_assembler}.vcf'
                    print(f'{vcf_file}', file=tmp_file_writer)

            tmp_file_writer.close()

            merged_out_vcf = f'{merged_outdir}/{caller}.{hifi_assembler}-{ont_assembler}.jasmine.merged.vcf'

            cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'

            os.system(cmd)

            os.remove(tmp_file)

            match_info_out = f'{merged_outdir}/{caller}.{hifi_assembler}-{ont_assembler}.dastsets-concordant.info.tsv'
            unique_info_out = f'{merged_outdir}/{caller}.{hifi_assembler}-{ont_assembler}.dastsets.unique.tsv'
            get_dataset_compare_info(datasets, merged_out_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds)

            hifi_uniques = f'{merged_outdir}/{caller}.{hifi_assembler}-{ont_assembler}.hifi.uniques.info.tsv'
            ont_uniques = f'{merged_outdir}/{caller}.{hifi_assembler}-{ont_assembler}.ont.uniques.info.tsv'

            get_platform_uniques(merged_out_vcf, hifi_uniques, ont_uniques, simple_reps, rmsk, sds)


def merge_datasets_read_at_regions(workdir, aligners, datasets, platform_datasets):
    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    # regions = ['highconf', 'cmrg']
    # for region in regions:
    #     for caller in CALLERS:
    #         for aligner in aligners:
    #             ## Matches of all SVs between platforms
    #             merged_outdir = f'{workdir}/read_dataset_repro/{region}_regions'
    #             if not os.path.exists(merged_outdir):
    #                 os.mkdir(merged_outdir)
    #
    #             tmp_file = f'{merged_outdir}/{caller}.{aligner}.txt'
    #             tmp_file_writer = open(tmp_file, 'w')
    #
    #             for dataset in datasets:
    #                 vcf_file = f'{workdir}/HiFi/{aligner}_{dataset}/filtered_in_highconf/HG002.{caller}.vcf'
    #                 if 'ont' in dataset:
    #                     vcf_file = f'{workdir}/ONT/{aligner}_{dataset}/filtered_in_highconf/HG002.{caller}.vcf'
    #                 print(f'{vcf_file}', file=tmp_file_writer)
    #             tmp_file_writer.close()
    #
    #             merged_out_vcf = f'{merged_outdir}/{caller}.{aligner}.jasmine.merged.vcf'
    #             print(f'Producing {caller} {aligner} merged calls ...')
    #             cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'
    #             os.system(cmd)
    #
    #             os.remove(tmp_file)
    #
    #             match_info_out = f'{merged_outdir}/{caller}.{aligner}.datasets-concordant.info.tsv'
    #             unique_info_out = f'{merged_outdir}/{caller}.{aligner}.datasets.unique.tsv'
    #             get_dataset_compare_info(datasets, merged_out_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds)


    for platform, datasets in platform_datasets.items():
        for caller in CALLERS:
            for aligner in aligners:
                ## Matches of all SVs between platforms
                merged_outdir = f'{workdir}/read_dataset_repro/{platform}/highconf_regions'
                if not os.path.exists(merged_outdir):
                    os.mkdir(merged_outdir)

                tmp_file = f'{merged_outdir}/{caller}.{aligner}.txt'
                tmp_file_writer = open(tmp_file, 'w')

                for dataset in datasets:

                    vcf_file = f'{workdir}/{PLATMAP[platform]}/{aligner}_{dataset}/filtered_in_highconf/HG002.{caller}.vcf'
                    print(f'{vcf_file}', file=tmp_file_writer)
                tmp_file_writer.close()

                merged_out_vcf = f'{merged_outdir}/{caller}.{aligner}.jasmine.merged.vcf'
                print(f'Producing {caller} {aligner} merged calls ...')
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'
                os.system(cmd)

                os.remove(tmp_file)

                match_info_out = f'{merged_outdir}/{caller}.{aligner}.{platform}-concordant.info.tsv'
                unique_info_out = f'{merged_outdir}/{caller}.{aligner}.{platform}.unique.tsv'
                get_dataset_compare_info(datasets, merged_out_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds)

def merge_datasets_assm_at_regions(workdir, aligners, datasets, platform_datasets):
    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    regions = ['highconf', 'cmrg']

    # for region in regions:
    #     for caller in ['pav', 'svimasm']:
    #         for aligner in aligners:
    #             ## Matches of all SVs between platforms
    #             merged_outdir = f'{workdir}/assm_dataset_repro/{region}_regions'
    #             if not os.path.exists(merged_outdir):
    #                 os.mkdir(merged_outdir)
    #
    #             tmp_file = f'{merged_outdir}/{caller}.{aligner}.txt'
    #             tmp_file_writer = open(tmp_file, 'w')
    #
    #             for dataset in datasets:
    #                 vcf_file = f'{workdir}/HiFi/{aligner}_{dataset}/filtered_in_highconf/HG002.{caller}.vcf'
    #                 if 'ont' in dataset:
    #                     vcf_file = f'{workdir}/ONT/{aligner}_{dataset}/filtered_in_highconf/HG002.{caller}.vcf'
    #                 print(f'{vcf_file}', file=tmp_file_writer)
    #             tmp_file_writer.close()
    #
    #             merged_out_vcf = f'{merged_outdir}/{caller}.{aligner}.jasmine.merged.vcf'
    #             print(f'Producing {caller} {aligner} merged calls ...')
    #             cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'
    #             os.system(cmd)
    #
    #             os.remove(tmp_file)
    #
    #             match_info_out = f'{merged_outdir}/{caller}.{aligner}.datasets-concordant.info.tsv'
    #             unique_info_out = f'{merged_outdir}/{caller}.{aligner}.datasets.unique.tsv'
    #             get_dataset_compare_info(datasets, merged_out_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds)


    for platform, datasets in platform_datasets.items():
        for caller in ['pav', 'svimasm']:
            for aligner in aligners:
                ## Matches of all SVs between platforms
                merged_outdir = f'{workdir}/assm_dataset_repro/{platform}/highconf_regions'
                if not os.path.exists(merged_outdir):
                    os.mkdir(merged_outdir)

                tmp_file = f'{merged_outdir}/{caller}.{aligner}.txt'
                tmp_file_writer = open(tmp_file, 'w')

                for dataset in datasets:

                    vcf_file = f'{workdir}/{PLATMAP[platform]}/{aligner}_{dataset}/filtered_in_highconf/HG002.{caller}.vcf'
                    print(f'{vcf_file}', file=tmp_file_writer)
                tmp_file_writer.close()

                merged_out_vcf = f'{merged_outdir}/{caller}.{aligner}.jasmine.merged.vcf'
                print(f'Producing {caller} {aligner} merged calls ...')
                cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'
                os.system(cmd)

                os.remove(tmp_file)

                match_info_out = f'{merged_outdir}/{caller}.{aligner}.{platform}-concordant.info.tsv'
                unique_info_out = f'{merged_outdir}/{caller}.{aligner}.{platform}.unique.tsv'
                get_dataset_compare_info(datasets, merged_out_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds)

def get_dataset_compare_info(datasets, merged_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds):

    matched_list = []
    unique_list = []
    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')

            info_tokens = entries[7].split(';')
            info_dict = {}

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            merged_type = info_dict['SVTYPE']
            merged_id = entries[2]
            supp_vec = info_dict['SUPP_VEC']
            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            if supp > 1:

                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)


                region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)

                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_datasets = datasets[supp_vec.index('1')]
                region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)
                unique_list.append((merged_id, merged_type, unique_datasets, entries[0], int(entries[1]), end, info_dict['SVLEN'], region_label, rptype, pcrt))


    # match_info_out = f'{workdir}/{sample}/{caller}/survivor/{sample}.{caller}.{aligner}.platform-concordant.{svtype}.tsv'
    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN', 'START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(match_info_out, header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns =['ID_MATCH', 'TYPE_MATCH', 'DATASET', '#CHROM', 'POS', 'END', 'SVLEN', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(unique_info_out, sep='\t', header=True, index=False)


def get_caller_compare_info(callers, merged_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds):
    matched_list = []
    unique_list = []

    svs_by_regions = {'Simple Repeats': [0, 0, 0], 'Repeat Masked': [0, 0, 0], 'Segment Dup': [0, 0, 0], 'Unique': [0, 0, 0]}

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')

            info_tokens = entries[7].split(';')
            info_dict = {}

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            merged_type = info_dict['SVTYPE']
            merged_id = entries[2]
            supp_vec = info_dict['SUPP_VEC']

            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            if supp == len(callers):
                start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
                start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)

                region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)
                matched_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), end, info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp == 1:
                unique_caller = callers[supp_vec.index('1')]
                region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)
                unique_list.append((merged_id, merged_type, unique_caller, entries[0], int(entries[1]),
                                    int(info_dict['END']), info_dict['SVLEN'], region_label, rptype, pcrt))

                if supp_vec == '01':
                    svs_by_regions[region_label][2] += 1
                if supp_vec == '10':
                    svs_by_regions[region_label][0] += 1

    df_matched = pd.DataFrame(matched_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN',
                                       'START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_matched.to_csv(match_info_out, header=True, sep='\t', index=False)

    df_uniques = pd.DataFrame(unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'CALLER', '#CHROM', 'POS', 'END', 'SVLEN',
                                       'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_uniques.to_csv(unique_info_out, sep='\t', header=True, index=False)

    return svs_by_regions

def get_platform_uniques(merged_vcf, hifi_info_out, ont_info_out, simple_reps, rmsk, sds):

    hifi_unique_list = []
    ont_unique_list = []
    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')

            info_tokens = entries[7].split(';')
            info_dict = {}

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            merged_type = info_dict['SVTYPE']
            merged_id = entries[2]
            supp_vec = info_dict['SUPP_VEC']

            start_var, end_var = abs(float(info_dict['STARTVARIANCE'])), abs(float(info_dict['ENDVARIANCE']))
            start_std, end_std = math.sqrt(start_var), math.sqrt(end_var)

            end = int(info_dict['END'])
            if merged_type == 'INS':
                end = int(entries[1]) + int(info_dict['SVLEN'])

            region_label, rptype, pcrt = annotate_sv_region(entries[0], int(entries[1]), end, 0, simple_reps, rmsk, sds)
            # HiFi uniques
            if supp_vec[3] == '0' and supp_vec[4] == '0' and supp_vec[5] == '0':
                hifi_unique_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]), int(info_dict['END']), info_dict['SVLEN'], start_std, end_std, region_label, rptype, pcrt))

            if supp_vec[0] == '0' and supp_vec[1] == '0' and supp_vec[2] == '0':
                ont_unique_list.append((merged_id, merged_type, supp_vec, supp, entries[0], int(entries[1]),int(info_dict['END']), info_dict['SVLEN'], start_std, end_std,region_label, rptype, pcrt))


    df_hifi_unqiue = pd.DataFrame(hifi_unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN', 'START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_hifi_unqiue.to_csv(hifi_info_out, header=True, sep='\t', index=False)

    df_ont_uniques = pd.DataFrame(ont_unique_list, columns=['ID_MATCH', 'TYPE_MATCH', 'SUPP_VEC', 'SUPP', '#CHROM', 'POS', 'END', 'SVLEN', 'START_STD', 'END_STD', 'REGION_TYPE', 'RPTYPE', 'PCRT'])
    df_ont_uniques.to_csv(ont_info_out, sep='\t', header=True, index=False)


'''
def merge_hifi_assm_calls(workdir, aligners, datasets):
    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    for caller in ['pav', 'svimasm']:
        for aligner in aligners:
            ## Matches of all SVs between platforms
            merged_outdir = f'{workdir}/assm_dataset_repro/hifi'
            if not os.path.exists(merged_outdir):
                os.mkdir(merged_outdir)

            tmp_file = f'{merged_outdir}/{caller}.{aligner}.txt'
            tmp_file_writer = open(tmp_file, 'w')

            for dataset in datasets:
                vcf_file = f'{workdir}/{aligner}_{dataset}/filtered/HG002.{caller}.vcf'
                print(f'{vcf_file}', file=tmp_file_writer)
            tmp_file_writer.close()

            merged_out_vcf = f'{merged_outdir}/{caller}.{aligner}.jasmine.merged.vcf'
            print(f'Producing {caller} {aligner} merged calls ...')
            cmd = f'{JASMINE} file_list={tmp_file} out_file={merged_out_vcf} max_dist=1000 spec_len=50 spec_reads=1'

            os.system(cmd)

            os.remove(tmp_file)

            match_info_out = f'{merged_outdir}/{caller}.{aligner}.hifi-concordant.info.tsv'
            unique_info_out = f'{merged_outdir}/{caller}.{aligner}.hifi.unique.tsv'
            get_dataset_info(datasets, merged_out_vcf, match_info_out, unique_info_out, simple_reps, rmsk, sds)
'''