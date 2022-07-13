#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/2/24

'''

import matplotlib.pyplot as plt
import pysam
import numpy as np


from helpers.Annot import *
from helpers.Functions import *
from com.CheckFDR import *

def get_merged_callers(workdir, datasets):
    ## Merge all callsets generated with dataset
    svtypes = ['ins', 'del']

    supps = [1, 2, 3, 4]
    for plat_supp in supps:
        for dataset in datasets:
            merged_outdir = f'{workdir}/merged_callers/{dataset}'
            if not os.path.exists(merged_outdir):
                os.mkdir(merged_outdir)

            supp_dir = f'{merged_outdir}/supp_{plat_supp}'
            if not os.path.exists(supp_dir):
                os.mkdir(supp_dir)

            for svtype in svtypes:
                tmp_file = f'{supp_dir}/{dataset}.{svtype}.txt'
                tmp_file_writer = open(tmp_file, 'w')
                for caller in CALLERS:
                    print(f'{workdir}/minimap2_{dataset}/HG002.{caller}.{svtype}.vcf', file=tmp_file_writer)

                tmp_file_writer.close()
                merged_out_vcf = f'{supp_dir}/{dataset}.{svtype}.{plat_supp}.merged.vcf'
                print(f'Merging {dataset} all callers {svtype} ...')
                cmd = f'{SURVIVOR} merge {tmp_file} 1000 {plat_supp} 1 0 0 50 {merged_out_vcf}'
                os.system(cmd)

                os.remove(tmp_file)


def caller_accumulate_merge(workdir, datasets):
    svtypes = ['ins', 'del']
    supps = [2, 3]

    for dataset in datasets:
        for svtype in svtypes:
            for plat_supp in supps:
                counter = 0
                merged_outdir = f'{workdir}/merged_callers/{dataset}/supp_{plat_supp}'

                vcf_to_merge = []
                merged_callers = []
                for i in range(len(CALLERS)):

                    merged_callers.append(CALLERS[i])
                    vcf_to_merge.append(f'{workdir}/minimap2_{dataset}/HG002.{CALLERS[i]}.{svtype}.vcf')
                    counter += 1

                    for j in range(i, len(CALLERS)):
                        if i != j:
                            merged_callers.append(CALLERS[j])
                            vcf_to_merge.append(f'{workdir}/minimap2_{dataset}/HG002.{CALLERS[j]}.{svtype}.vcf')
                            counter += 1
                            if counter == plat_supp:
                                merged_caller_tag = '-'.join(merged_callers)
                                print(f'{svtype}: {merged_caller_tag}')

                                tmp_file = f'{merged_outdir}/{dataset}.{merged_caller_tag}.txt'
                                tmp_file_writer = open(tmp_file, 'w')

                                for vcf in vcf_to_merge:
                                    print(vcf, file=tmp_file_writer)

                                tmp_file_writer.close()

                                merged_vcf = f'{merged_outdir}/{dataset}.{svtype}.{merged_caller_tag}.merged.vcf'
                                cmd = f'{SURVIVOR} merge {tmp_file} 1000 {plat_supp} 1 0 0 50 {merged_vcf}'
                                os.system(cmd)

                                counter -= 1
                                vcf_to_merge = vcf_to_merge[0: -1]
                                merged_callers = merged_callers[0: -1]

                        if j == len(CALLERS) - 1:
                            counter = 0
                            vcf_to_merge = []
                            merged_callers = []

def compare_merged_insdel_to_asm(workdir, datasets, ctg_dataset, simrep_path, rmsk_path, sd_path):
    simple_reps = pysam.Tabixfile(simrep_path, 'r')
    rmsk = pysam.Tabixfile(rmsk_path, 'r')
    sds = pysam.Tabixfile(sd_path, 'r')

    supps = [1, 2, 3, 4]
    svtypes = ['ins', 'del']
    for dataset in datasets:
        for svtype in svtypes:
            for plat_supp in supps:
                pav_calls = f'{workdir}/{ctg_dataset}/pav/{ctg_dataset}.{svtype}.vcf'

                caller_merged = f'{workdir}/merged_callers/{dataset}/supp_{plat_supp}/{dataset}.{svtype}.{plat_supp}.merged.vcf'
                compare_outdir = f'{workdir}/merged_callers/{dataset}/supp_{plat_supp}/com'

                if not os.path.exists(compare_outdir):
                    os.mkdir(compare_outdir)

                merged_vcf = f'{compare_outdir}/{dataset}.pav.{svtype}.{plat_supp}.merged.vcf'

                if not os.path.exists(merged_vcf):
                    tmp_file = f'{compare_outdir}/{dataset}.txt'
                    tmp_file_writer = open(tmp_file, 'w')
                    print(pav_calls, file=tmp_file_writer)
                    print(caller_merged, file=tmp_file_writer)

                    tmp_file_writer.close()

                    cmd = f'{SURVIVOR} merge {tmp_file} 1000 1 1 0 0 50 {merged_vcf}'
                    os.system(cmd)

                    mat_out = f'{compare_outdir}/{dataset}.pav.{svtype}.{plat_supp}.mat.txt'
                    cmd = f'{SURVIVOR} genComp {merged_vcf} 0 {mat_out}'
                    os.system(cmd)

                    os.remove(tmp_file)

                rep_annotate_uniques(compare_outdir, SVTYPEMAP[svtype], dataset, merged_vcf, simple_reps, rmsk, sds)

def accumlate_merged_to_asm(workdir, datasets, ctg_dataset):

    supps = [2, 3]
    svtypes = ['ins', 'del']
    for dataset in datasets:
        for svtype in svtypes:
            for plat_supp in supps:
                merged_outdir = f'{workdir}/merged_callers/{dataset}/supp_{plat_supp}'
                counter = 0
                merged_callers = []
                for i in range(len(CALLERS)):
                    merged_callers.append(CALLERS[i])
                    counter += 1
                    for j in range(i, len(CALLERS)):
                        if i != j:
                            merged_callers.append(CALLERS[j])
                            counter += 1
                            if counter == plat_supp:
                                merged_caller_tag = '-'.join(merged_callers)
                                print(f'{svtype}: {merged_caller_tag}')
                                caller_merged_vcf = f'{merged_outdir}/{dataset}.{svtype}.{merged_caller_tag}.merged.vcf'

                                pav_calls = f'{workdir}/{ctg_dataset}/pav/{ctg_dataset}.{svtype}.vcf'
                                compare_outdir = f'{merged_outdir}/com/caller_combine'

                                if not os.path.exists(compare_outdir):
                                    os.mkdir(compare_outdir)

                                merged_vcf = f'{compare_outdir}/{dataset}.pav.{svtype}.{merged_caller_tag}.merged.vcf'

                                tmp_file = f'{compare_outdir}/{dataset}.{merged_caller_tag}.txt'
                                tmp_file_writer = open(tmp_file, 'w')
                                print(pav_calls, file=tmp_file_writer)
                                print(caller_merged_vcf, file=tmp_file_writer)

                                tmp_file_writer.close()

                                cmd = f'{SURVIVOR} merge {tmp_file} 1000 1 1 0 0 50 {merged_vcf}'
                                os.system(cmd)

                                mat_out = f'{compare_outdir}/{dataset}.pav.{svtype}.{merged_caller_tag}.mat.txt'
                                cmd = f'{SURVIVOR} genComp {merged_vcf} 0 {mat_out}'
                                os.system(cmd)
                                os.remove(tmp_file)

                                counter -= 1
                                merged_callers = merged_callers[0: -1]

                        if j == len(CALLERS) - 1:
                            counter = 0
                            merged_callers = []