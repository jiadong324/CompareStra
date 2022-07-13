#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2021/9/28

'''

import pandas as pd
import vcf
import gzip

from helpers.Functions import *

def read_survivor_vcf(survivor_vcf, max_size):

    exbnd_list = []
    all_sv_list = []
    for line in open(survivor_vcf, 'r'):
        if '#' in line:
            continue
        entries = line.strip().split('\t')

        chrom, start, id = entries[0], int(entries[1]), entries[2]
        info_tokens = entries[7].split(';')

        if 'chr' not in chrom:
            chrom = f'chr{chrom}'

        svtype = info_tokens[3].split('=')[1]

        if chrom not in VALID_CHROMS:
            continue

        supp_vec = info_tokens[1].split('=')[1]

        if svtype != 'TRA':
            svlen = abs(int(info_tokens[2].split('=')[1]))
            if svlen > max_size:
                continue

            end = int(info_tokens[6].split('=')[1])

            if start > end:
                end = start + svlen

            all_sv_list.append((chrom, start, end, svtype, svlen, id))
            exbnd_list.append((chrom, start, end, svtype, supp_vec, svlen, id))

        all_sv_list.append((chrom, start, start + 1, svtype, supp_vec, -1, id))


    df_exbnd = pd.DataFrame(exbnd_list, columns=['chrom', 'start', 'end', 'svtype', 'supp', 'svlen', 'svid'])
    df_allsvs = pd.DataFrame(all_sv_list, columns=['chrom', 'start', 'end', 'svtype', 'supp', 'svlen', 'svid'])

    return df_exbnd, df_allsvs

def read_survivor_gencomp(file_path):

    values = []

    with open(file_path, 'r') as f:
        for line in f:
            values.append([int(val) for val in line.strip().split('\t')])

    return values


def get_survivor_supp(merged_vcf):

    supp_dict = {}
    merged_total = 0
    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')

            info_tokens = entries[7].split(";")
            info_dict = {}
            merged_total += 1

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = int(info_dict['SUPP'])
            if supp in supp_dict:
                supp_dict[supp] += 1
            else:
                supp_dict[supp] = 1

    return supp_dict, merged_total

def get_survivor_suppvec(merged_vcf):

    suppvec_dict = {}
    merged_total = 0
    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split('\t')
            info_tokens = entries[7].split(";")
            info_dict = {}
            merged_total += 1

            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            supp = info_dict['SUPP_VEC']
            if supp in suppvec_dict:
                suppvec_dict[supp] += 1
            else:
                suppvec_dict[supp] = 1

    return suppvec_dict, merged_total

def separate_stra_merged_vcfs(merged_vcf, asm_unique_out, read_unique_out):

    asm_unique_writer = open(asm_unique_out, 'w')
    read_unique_writer = open(read_unique_out, 'w')
    read_counter = 0
    assm_counter = 0

    with open(merged_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                print(line.strip(), file=asm_unique_writer)
                print(line.strip(), file=read_unique_writer)
                continue

            entries = line.strip().split('\t')
            info_tokens = entries[7].split(';')
            info_dict = {}
            for token in info_tokens:
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            if info_dict['SUPP_VEC'] == '10':
                print(line.strip(), file=asm_unique_writer)
                assm_counter += 1

            elif info_dict['SUPP_VEC'] == '01':
                print(line.strip(), file=read_unique_writer)
                read_counter += 1

    print('Read-unique:', read_counter)
    print('Assembly-unique:', read_counter)

    asm_unique_writer.close()
    read_unique_writer.close()


'''

def read_vcf(input_vcf, caller, minsr, max_size, exclude_dict, ref_file):
    re_tag_caller = ['sniffles', 'cutesv']
    sv_exbnd_list = []
    all_sv_list = []
    filtered_sv_list = []

    with open(input_vcf, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            entries = line.strip().split("\t")
            chrom, sv_id = entries[0], entries[2]

            if 'chr' not in chrom:
                chrom = f'chr{chrom}'

            start = int(entries[1])

            if chrom not in VALID_CHROMS:
                continue

            info_tokens = entries[7].split(";")
            info_dict = {}

            for token in info_tokens:
                if "=" not in token:
                    continue
                info_dict[token.split("=")[0]] = token.split("=")[1].replace(">","")

            if minsr != -1:

                if caller in re_tag_caller and int(info_dict['RE']) < minsr:
                    continue

                if caller == 'svim' and int(info_dict['SUPPORT']) < minsr:
                    continue

            try:
                sv_type = info_dict["SVTYPE"]

                if sv_type not in ['BND', 'TRA']:
                    end = int(info_dict["END"])
                    sv_len = end - start
                    if "SVLEN" in info_dict:
                        sv_len = abs(int(info_dict["SVLEN"]))

                    if sv_len < 50 or sv_len >= max_size:
                        filtered_sv_list.append((chrom, start, end, sv_len))
                        continue

                    if "INS" in sv_type:
                        end = start + sv_len

                    if exclude_dict[chrom].overlap(start, end):
                        continue

                    if contains_gaps(chrom, start, end, ref_file):
                        continue

                    sv_exbnd_list.append((chrom, start, end, sv_type, sv_len, sv_id))

                all_sv_list.append((chrom, start, start + 1, sv_type, -1, sv_id))
            except:
                print(line)

    df_svs_exbnd = pd.DataFrame(sv_exbnd_list, columns=['chrom', 'start', 'end', 'svtype', 'svlen', 'svid'])
    df_allsvs = pd.DataFrame(all_sv_list, columns=['chrom', 'start', 'end', 'svtype', 'svlen', 'svid'])

    return df_svs_exbnd, df_allsvs
'''

