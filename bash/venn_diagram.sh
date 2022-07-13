#!/bin/bash

sample=$1

export PATH=/Users/apple/miniconda3/envs/dnatools/bin:$PATH
export PATH=/Users/apple/miniconda3/envs/r-env/bin:$PATH


declare -a callers=("svision" "pbsv" "svim" "sniffles" "cutesv" "nanovar" )

#for caller in ${callers[@]}; do
#
#cd /Users/apple/Evaluation/$sample
#
#echo $caller
#intervene venn -i ./pav_v112/"pav_"$sample"_CCS_svs.bed" ./$caller/bed/$sample".hifi.minimap2."$caller".s5.filtered.bed" ./$caller/bed/$sample".hifi.ngmlr."$caller".s5.filtered.bed" ./$caller/bed/$sample".ont.minimap2."$caller".s5.filtered.bed" ./$caller/bed/$sample".ont.ngmlr."$caller".s5.filtered.bed" --names=PAV,HiFi-minimap2,HiFi-ngmlr,ONT-minimap2,ONT-ngmlr -o ./pav_v112/intervene/$caller --bedtools-options f=0.5,r --title $caller"-50%RO"
#
#intervene upset -i ./$caller/bed/$sample".hifi.minimap2."$caller".bedtools.missed.pavs.bed" ./$caller/bed/$sample".hifi.ngmlr."$caller".bedtools.missed.pavs.bed" ./$caller/bed/$sample".ont.minimap2."$caller".bedtools.missed.pavs.bed" ./$caller/bed/$sample".ont.ngmlr."$caller".bedtools.missed.pavs.bed" --names=HiFi-minimap2,HiFi-ngmlr,ONT-minimap2,ONT-ngmlr --bedtools-options f=0.5,r -o /Users/apple/Evaluation/$sample/pav_v112/intervene/$caller --save-overlaps
#
#cd /Users/apple/Evaluation/$sample/pav_v112/intervene/$caller
#Rscript Intervene_upset.R
#
#done

declare -a platforms=("hifi" "ont" )

for plat in ${platforms[@]}; do

cd /Users/apple/Somatic/$sample/pav_v112

echo $plat
intervene upset -i ../svision/bed/$sample"."$plat".minimap2.svision.bedtools.missed.pavs.bed" ../pbsv/bed/$sample"."$plat".minimap2.pbsv.bedtools.missed.pavs.bed" ../cutesv/bed/$sample"."$plat".minimap2.cutesv.bedtools.missed.pavs.bed" ../svim/bed/$sample"."$plat".minimap2.svim.bedtools.missed.pavs.bed" ../sniffles/bed/$sample"."$plat".minimap2.sniffles.bedtools.missed.pavs.bed" ../nanovar/bed/$sample"."$plat".minimap2.nanovar.bedtools.missed.pavs.bed" -o /Users/apple/Somatic/$sample/pav_v112/$plat"_minimap2_missed_pav_upset" --names=SVision,PbSV,CuteSV,SVIM,Sniffles,NanoVar --bedtools-options f=0.5,r --save-overlaps

cd /Users/apple/Somatic/$sample/pav_v112/$plat"_minimap2_missed_pav_upset"
Rscript Intervene_upset.R

cd /Users/apple/Somatic/$sample/pav_v112

intervene upset -i ../svision/bed/$sample"."$plat".ngmlr.svision.bedtools.missed.pavs.bed" ../pbsv/bed/$sample"."$plat".ngmlr.pbsv.bedtools.missed.pavs.bed" ../cutesv/bed/$sample"."$plat".ngmlr.cutesv.bedtools.missed.pavs.bed" ../svim/bed/$sample"."$plat".ngmlr.svim.bedtools.missed.pavs.bed" ../sniffles/bed/$sample"."$plat".ngmlr.sniffles.bedtools.missed.pavs.bed" ../nanovar/bed/$sample"."$plat".ngmlr.nanovar.bedtools.missed.pavs.bed" -o /Users/apple/Somatic/$sample/pav_v112/$plat"_ngmlr_missed_pav_upset" --bedtools-options f=0.5,r --names=SVision,PbSV,CuteSV,SVIM,Sniffles,NanoVar --save-overlaps
cd /Users/apple/Somatic/$sample/pav_v112/$plat"_ngmlr_missed_pav_upset"
Rscript Intervene_upset.R

done