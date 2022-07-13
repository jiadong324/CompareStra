#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2021/10/17

'''
import sys


def parse_sniffles(input_vcf, output_vcf):

    out_writer = open(output_vcf, 'w')

    for line in open(input_vcf, 'r'):
        if '#' in line:
            print(line.strip(), file=out_writer)
            continue

        entries = line.strip().split('\t')
        chrom, start = entries[0], int(entries[1])
        info_tokens = entries[7].split(';')[1:]
        info_dict = {}

        # if 'STRANDBIAS' in line :
        #     continue

        for token in info_tokens:
            info_dict[token.split('=')[0]] = token.split('=')[1]

        if 'END' in info_dict:
            if int(info_dict['END']) < start:
                start = int(info_dict['END'])
                new_entry = [chrom, str(start)]
                new_entry.extend(entries[2:])
                new_entry_out = '\t'.join(new_entry)
                print(new_entry_out, file=out_writer)

                continue

        print(line.strip(), file=out_writer)

def parse_svision(input_vcf, output_vcf, minsr):
    out_writer = open(output_vcf, 'w')

    for line in open(input_vcf, 'r'):
        if '#' in line:
            print(line.strip(), file=out_writer)
            continue

        entries = line.strip().split('\t')
        info_tokens = entries[7].split(';')
        info_dict = {}

        for token in info_tokens:
            info_dict[token.split('=')[0]] = token.split('=')[1]

        if int(info_dict['SUPPORT']) < minsr:
            print(line)
            continue

        if entries[6] == 'Covered':
            new_entries = entries[0: 6]
            new_entries.append('PASS')
            new_entries.extend(entries[7: ])
            new_entry_out = '\t'.join(new_entries)

            print(new_entry_out, file=out_writer)


def parse_nanovar(input_vcf, output_vcf):
    out_writer = open(output_vcf, 'w')

    for line in open(input_vcf, 'r'):
        if '#' in line:
            print(line.strip(), file=out_writer)
            continue

        try:
            entries = line.strip().split('\t')
            info_tokens = entries[7].split(';')[1:]
            info_dict = {}

            for token in info_tokens:
                info_dict[token.split('=')[0]] = token.split('=')[1]

            length = abs(int(info_dict['SVLEN']))
            print(line.strip(), file=out_writer)

        except:
            print(line.strip())


def main():

    # sniffles_input = '/Users/apple/Evaluation/HG002/minimap2/hifi/HG002.sniffles.vcf'
    # sniffles_output = '/Users/apple/Evaluation/HG002/minimap2/hifi/HG002.sniffles.tmp.vcf'
    #
    # # parse_sniffles(sniffles_input, sniffles_output)
    #
    # nanovar_input = '/Users/apple/Evaluation/HG002/minimap2/hifi/HG002.hifi.minimap2.sorted.nanovar.pass.vcf'
    # nanovar_output = '/Users/apple/Evaluation/HG002/minimap2/hifi/HG002.nanovar.vcf'
    #
    # # parse_nanovar(nanovar_input, nanovar_output)

    tool = sys.argv[1]
    input_vcf = sys.argv[2]
    output_vcf = sys.argv[3]
    minsr = sys.argv[4]

    if tool == 'sniffles':
        parse_sniffles(input_vcf, output_vcf)

    elif tool == 'nanovar':
        parse_nanovar(input_vcf, output_vcf)

    elif tool == 'svision':
        parse_svision(input_vcf, output_vcf, int(minsr))

if __name__ == '__main__':
    main()