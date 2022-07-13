#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2021/10/23
'''


## iMac work path
iMAC = '/Users/apple/Evaluation'
iMACDIR = f'{iMAC}/HG002/grch37'
iMACBAM = f'/Volumes/mac/Data/HG002'
# iMACREF = '/Users/apple/Data/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'
# iMACEXCLUDE = '/Users/apple/Data/genome/grch38.exclude_regions_cen.bed'

iMACREF = '/Users/apple/Data/genome/hs37d5.fa'
iMACEXCLUDE = '/Users/apple/Data/genome/repeat_annot/grch37/grch37.exclude_regions_cen.bed'
iMACCENTRO = '/Users/apple/Data/genome/repeat_annot/grch37/centromeres.txt'
iMACSIMPLE='/Users/apple/Data/genome/repeat_annot/grch37/simplerepeat.bed.gz'
iMACRMSK='/Users/apple/Data/genome/repeat_annot/grch37/rmsk.bed.gz'
iMACSD = '/Users/apple/Data/genome/repeat_annot/grch37/seg_dup.bed.gz'
CHROMSIZE = '/Users/apple/Data/genome/hs37d5.fa.fai'


READFIGURE = f'{iMAC}/hg002_figures/read'
iMACCOMFIGURE = f'{iMAC}/hg002_figures/com'
ASMFIGURE = f'{iMAC}/hg002_figures/asm'

Figure1 = f'{iMAC}/hg002_figures/Figure1'
Figure2 = f'{iMAC}/hg002_figures/Figure2'
Figure3 = f'{iMAC}/hg002_figures/Figure3'
Figure4 = f'{iMAC}/hg002_figures/Figure4'
Figure5 = f'{iMAC}/hg002_figures/Figure5'
Figure6 = f'{iMAC}/hg002_figures/Figure6'

SUPPFIG = f'{iMAC}/hg002_figures/SuppFigs'

## MacBook work path
# MBP = '/Users/jiadonglin/Evaluation'
# MBPHG002DIR = f'{MBP}/HG002/grch37'
# MBPREF = '/Users/jiadonglin/Data/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'
# MBPEXCLUDE = '/Users/jiadonglin/Data/genome/grch38.exclude_regions_cen.bed'
# MBPSIMPLE='/Users/jiadonglin/Data/genome/repeat_annot/grch37/simplerepeat.bed.gz'
# MBPRMSK='/Users/jiadonglin/Data/genome/repeat_annot/grch37/rmsk.bed.gz'
#
# MBPHG002FIGURE = f'{MBP}/hg002_figures'
# MBPHGSVCFIGURE = f'{MBP}/hgsvc_figures'

## Plot params
TOOLMARKERS = {'pav': 'P', 'pbsv': 'X', 'svim': 'd', 'cutesv': 'o', 'sniffles': 's', 'svision': 'H',
               'svimasm': 'p', 'dipcall': '*'}

TOOLMAP = {'pav': 'PAV', 'pbsv': 'pbsv', 'svim':'SVIM', 'cutesv': 'cuteSV', 'svision': 'SVision', 'sniffles':
    'Sniffles', 'nanovar': 'NanoVar', 'com': 'PAV', 'svimasm': 'SVIM-asm', 'dipcall': 'dipcall'}

TOOLCOLORS = {'svimasm': '#c078aa', 'pbsv': '#d95f02', 'svim': '#7570b3', 'cutesv': '#d51dad', 'svision': '#008066',
              'sniffles': '#386cb0', 'pav': '#c88f1c', 'dipcall': '#e41a1c'}

# REPCOLOR = {'VNTR': '#1b9e77', 'STR': '#d95f02', 'LTR': '#7570b3', 'SINE': '#e7298a', 'LINE': '#386cb0'}
COMCOLOR = {'Assm': '#e29a5e', 'Intersect': '#458bc2', 'Align': '#70a97a'}
ALIGNERCOLOR = {'minimap2': '#ce1181', 'ngmlr': '#1b9e77', 'lra': '#d1664d', 'winnowmap': '#85abd3'}
# ALIGNERLS = {'minimap2': '--', 'ngmlr': '-'}
ALIGNERMK = {'minimap2': 's', 'ngmlr': 'X', 'lra': 'p', 'winnowmap': 'd'}

PLATLS = {'ONT': '--', 'HiFi': '-', 'ontpro': '-.'}
PLATHATCH = {'ONT': '///', 'HiFi':''}
# PLATMARKER = {'ont': 's', 'hifi': 'X', 'ontpro': 'o'}

PLATCOLORS = {'hifi': '#fa850f', 'hifi_10kb': '#fee0c3', 'hifi_11kb': '#fdc287', 'hifi_15kb': '#fc9263',
              'hifi_18kb': '#ba5904', 'minion_27X': '#5691c1',
              'promethion_27X': '#78c679', 'ont': '#4280b3', 'ontpro': '#78c679'
              ,'HiFi': '#cb6627', 'ONT': '#7471ae'}

PLATMAP = {'hifi': 'HiFi', 'ont': 'ONT', 'ontpro':'PromethION', 'hifi_10kb': 'HiFi-10kb', 'hifi_11kb': 'HiFi-11kb', 'hifi_15kb': 'HiFi-15kb',
            'hifi_18kb': 'HiFi-18kb', 'ont_9kb': 'ONT-9kb', 'ont_19kb': 'ONT-19kb', 'ont_30kb': 'ONT-30kb',
           'giab_ctg': 'PAV-GIAB', 'hgsvc_ctg': 'PAV-HGSVC'}

REGIONCOLORS = {'Simple Repeats': '#5691c1', 'Repeat Masked': '#8e7cc3', 'Segment Dup': '#78c679', 'Unique': '#fa850f'}
REGIONMARKERS = {'Simple Repeats': 'o', 'Repeat Masked': 'X', 'Segment Dup': 's', 'Unique': 'd'}
# STRACOLORS = {'Assembly': '#ac5a54', 'Read': '#77a99f'}
STRACOLORS = {'Assembly': '#cf2e33', 'Read': '#389b46'}

SVTYPEMAP = {'ins': 'INS', 'del': 'DEL', 'others': 'Others', 'inv': 'INV', 'dup': 'DUP'}
CALLERS = ['pbsv', 'svim', 'cutesv', 'sniffles', 'svision']
ASMCALLERS = ['pav', 'svimasm']

CMRG = '/Users/apple/Evaluation/HG002/CMRGs/HG002_GRCh37_CMRG_SV_v1.00.bed'
HIGHCONF = '/Users/apple/Evaluation/HG002/HG002_SVs_Tier1_v0.6.pass.bed'

# PLATDICT = {'hifi': 'HiFi', 'ont': 'ONT'}
# SAMPLECOLORS = {'HG00733': '#29a249', 'HG00514': '#1878b6', 'NA19240': '#f57f20'}

SVTYPECOLORS = {'DEL': '#dc9049', 'DUP': '#54b06d', 'INS': '#5090dd', 'INV': '#f2ed4f', 'BND': '#c265e7', 'Others': '#CDB699', 'DUP:TANDEM': '#CDB699'}
# SVTYPECOLORS = {'DEL': '#ff7800', 'DUP': '#00a330', 'INS': '#007bba', 'INV': '#f2ed4f', 'BND': '#c265e7', 'Others': '#CDB699', 'DUP:TANDEM': '#CDB699'}

BPSHIFTCOLORS = ['#2f87d3', '#c07da4', '#dba038', '#878787']
custom_params = {"axes.spines.right": False, "axes.spines.top": False}

SHIFTCOLORDICT = {'0,10':'#c07da4', '10,50': '#dba038', '>50':'#878787', '0': '#2f87d3'}

## Tool path
SURVIVOR= '~/Biotools/SURVIVOR/Debug/SURVIVOR'
JASMINE = '/Users/apple/miniconda3/envs/dnatools/bin/jasmine'


VALID_CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                  "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

AUTOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                  "13", "14", "15", "16", "17", "18", "19", "20", "21", "22"]

ASSEMBLER = ['flye-flye', 'flye-shasta', 'hifiasm-flye', 'hifiasm-shasta']
ASSEMBLERMK = {'flye-flye': 'o', 'flye-shasta': 'D', 'hifiasm-flye': 'X', 'hifiasm-shasta': 's'}

ASSMBLERCOLOR = {'flye': '#005f86', 'hifiasm': '#c06600', 'shasta': '#ff3c4d'}
