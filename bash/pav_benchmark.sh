#!/bin/bash

sample=$1
aligner=$2
platform=$3

pav=/Users/apple/Evaluation/$sample/pav_v112/pav_"$sample"_CCS_svs.bed
bench=$(awk '{count++} END {print count}' $pav)

comout=/Users/apple/Evaluation/$sample/pav_v112/"$platform"_"$aligner"_s5_benchmark_results.txt

if [ -f "$comout" ];then
rm "$comout"
fi

#declare -a callers=("svision" "pbsv" "svim" "sniffles" "cutesv" "nanovar" )
declare -a callers=("pbsv" "svim" "sniffles" "cutesv" "nanovar" )

for caller in ${callers[@]}; do

dir=/Users/apple/Evaluation/$sample/$caller/bed


calls=$dir/$sample.$platform.$aligner.$caller".s5.exbnd.bed"

echo $calls

# all calls
correct=$(bedtools intersect -c -a $pav -b $calls -f 0.5 -r | awk '$6>0 {count++} END {print count}')
all=$(awk '{count++} END {print count}' $calls)

echo | awk -v a="$correct" -v b="$bench" -v c="$correct" -v d="$all" -v e="$caller" '{printf "%s,%s,%s\n", e, a/b, c/d}' >> $comout

bedtools intersect -c -a $pav -b $calls -f 0.5 -r | awk '$6==0 {print}' > $dir/$sample.$platform.$aligner.$caller".bedtools.missed.pavs.bed"

done
