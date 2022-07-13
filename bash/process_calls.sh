#!/bin/bash

sample=$1
minsr=$2
platform=$3

cd /Users/apple/CellLines/$sample/svision/
python ~/PycharmProjects/CellLine/Parser.py vcf2bed -v $sample".svision.s5.graph.vcf" -b $sample".svision.s5.exclude.bed" -t svision -s $minsr -x ~/Data/genome/grch38.exclude_regions_cen.bed
#python ~/PycharmProjects/SVision/SVisionFilter.py -v $sample".svision.s5.graph.vcf" -g $sample".graph_exactly_match.txt" -r /Users/apple/Data/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -o ./ -i 0,27 -e ~/Data/genome/grch38.exclude_regions_cen.bed

cd /Users/apple/CellLines/$sample/pbsv/
python ~/PycharmProjects/CellLine/Parser.py vcf2bed -v $sample.$platform".ngmlr.pbsv.s5.vcf" -b $sample".pbsv.s5.exclude.bed" -t pbsv -s $minsr -x ~/Data/genome/grch38.exclude_regions_cen.bed

cd /Users/apple/CellLines/$sample/svim/
python ~/PycharmProjects/CellLine/Parser.py vcf2bed -v $sample.$platform".ngmlr.svim.s5.vcf" -b $sample".svim.s5.exclude.bed" -t svim -s $minsr -x ~/Data/genome/grch38.exclude_regions_cen.bed

cd /Users/apple/CellLines/$sample/sniffles/
python ~/PycharmProjects/CellLine/Parser.py vcf2bed -v $sample.$platform".ngmlr.sniffles.s5.vcf" -b $sample".sniffles.s5.exclude.bed" -t sniffles -s $minsr -x ~/Data/genome/grch38.exclude_regions_cen.bed

cd /Users/apple/CellLines/$sample/cutesv/
python ~/PycharmProjects/CellLine/Parser.py vcf2bed -v $sample.$platform".ngmlr.cutesv.s5.vcf" -b $sample".cutesv.s5.exclude.bed" -t cutesv -s $minsr -x ~/Data/genome/grch38.exclude_regions_cen.bed


