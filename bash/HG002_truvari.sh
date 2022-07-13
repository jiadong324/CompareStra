#!/bin/bash

aligner=$1
platform=$2

#export $PATH=/Users/apple/miniconda3/envs/dnatools/bin:$PATH

cd /Users/apple/Somatic/HG002/$aligner/$platform

## Evaluate SVision

#python ~/PycharmProjects/Evaluation/ParseVcf.py svision HG002.svision.s10.vcf tmp.svision.vcf 10
#bcftools sort ./tmp.svision.vcf > ./HG002.svision.vcf
#bgzip -f ./HG002.svision.vcf
#tabix ./HG002.svision.vcf.gz
#
#if [ -d bench_svision ];then
#rm -rf bench_svision
#fi
#
#truvari bench -f ~/Data/genome/hs37d5.fa -b ~/Evaluation/HG002/HG002_SVs_Tier1_v0.6.vcf.gz --includebed ~/Evaluation/HG002/HG002_SVs_Tier1_v0.6.bed -o ./bench_svision --giabreport --passonly -r 1000 -p 0.00 -c HG002.svision.vcf.gz

## Evaluate Sniffles

#python ~/PycharmProjects/Evaluation/ParseVcf.py sniffles HG002.sniffles.vcf tmp.sniffles.vcf 10
#cat <(cat tmp.sniffles.vcf | grep "^#") <(cat tmp.sniffles.vcf | grep -vE "^#" | grep 'DUP\|INS\|DEL' | sed 's/DUP/INS/g' | sort -k1,1 -k2,2g) | bgzip -f -c > HG002.sniffles.vcf.gz
#tabix HG002.sniffles.vcf.gz
#
#if [ -d bench_sniffles ];then
#rm -rf bench_sniffles
#fi
#
#truvari bench -f ~/Data/genome/hs37d5.fa -b ~/Evaluation/HG002/HG002_SVs_Tier1_v0.6.vcf.gz --includebed ~/Evaluation/HG002/HG002_SVs_Tier1_v0.6.bed -o ./bench_sniffles --giabreport --passonly -r 1000 -p 0.00 -c HG002.sniffles.vcf.gz


## Evaluate Nanovar

python ~/PycharmProjects/Somatic/ParseVcf.py nanovar HG002.$platform.$aligner.sorted.nanovar.pass.vcf HG002.nanovar.vcf 10
bgzip -f HG002.nanovar.vcf
tabix HG002.nanovar.vcf.gz

if [ -d bench_nanovar ];then
rm -rf bench_nanovar
fi

truvari bench -f ~/Data/genome/hs37d5.fa -b ~/Somatic/HG002/HG002_SVs_Tier1_v0.6.vcf.gz --includebed ~/Somatic/HG002/HG002_SVs_Tier1_v0.6.bed -o ./bench_nanovar --giabreport --passonly -r 1000 -p 0.00 -c HG002.nanovar.vcf.gz


## Evaluate pbsv

#bgzip -f HG002.pbsv.vcf
#tabix HG002.pbsv.vcf.gz
#
#if [ -d bench_pbsv ];then
#rm -rf bench_pbsv
#fi
#
#truvari bench -f ~/Data/genome/hs37d5.fa -b ~/Evaluation/HG002/HG002_SVs_Tier1_v0.6.vcf.gz --includebed ~/Evaluation/HG002/HG002_SVs_Tier1_v0.6.bed -o ./bench_pbsv --giabreport --passonly -r 1000 -p 0.00 -c HG002.pbsv.vcf.gz

## Evalute SVIM

#bgzip -f HG002.svim.vcf
#tabix HG002.svim.vcf.gz
#
#if [ -d bench_svim ];then
#rm -rf bench_svim
#fi
#
#truvari bench -f ~/Data/genome/hs37d5.fa -b ~/Evaluation/HG002/HG002_SVs_Tier1_v0.6.vcf.gz --includebed ~/Evaluation/HG002/HG002_SVs_Tier1_v0.6.bed -o ./bench_svim --giabreport --passonly -r 1000 -p 0.00 -c HG002.svim.vcf.gz


## Evaluate cutesv

#bgzip -f HG002.cutesv.s10.vcf
#tabix HG002.cutesv.s10.vcf.gz
#
#if [ -d bench_cutesv ];then
#rm -rf bench_cutesv
#fi
#
#truvari bench -f ~/Data/genome/hs37d5.fa -b ~/Evaluation/HG002/HG002_SVs_Tier1_v0.6.vcf.gz --includebed ~/Evaluation/HG002/HG002_SVs_Tier1_v0.6.bed -o ./bench_cutesv --giabreport --passonly -r 1000 -p 0.00 -c HG002.cutesv.s10.vcf.gz


rm tmp.*