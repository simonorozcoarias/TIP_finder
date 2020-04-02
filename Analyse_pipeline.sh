#####
#Date : 20 décembre 2016
# Analyse automatique des résultats du pipeline 3000genomes
#####

#!/usr/bin/bash

#threshold 2 instead of threashold 5
# for file in *sort.bed;do out=$(echo $file | sed -e "s/\.sort/\_per10kb/");bedtools coverage -counts -b $file -a  ~/Bureau/Database/IRGSP-1.0_10kbpwindows.bed | awk  -F "\t" '{if ($4>=2){print $0}}' > coveragebed\_$out;done



#nom de la famille te analysé
te=HERVK
echo $te

#on analyse seulemement coveragebed files
mkdir final_cov2
mv coveragebed_* final_cov2
cd final_cov2

#recupeartion de toutes les insertions du te donné
for file in *bed;do awk -F "\t" '{print $1"_"$2"_"$3}' $file;done | sort -u > all_insertion_$te.names
#nb tot fenetre insertion
wc -l all_insertion_$te.names

#recupeartion de toutes les insertions du te donné avec seuil 5
for file in coveragebed_*;do awk -F "\t" '{if ($4>=5){print $1"_"$2"_"$3}}' $file;done | sort -u  > all_position_cov5_$te.names
wc -l all_position_cov5_$te.names

#reformatage des sorties
for file in *bed;do n=$(echo $file | sed -e "s/coveragebed_//" | sed -e "s/\-vs-.*_per10kb.bed//"); awk -F "\t" '{print $1"_"$2"_"$3}' $file > $n.txt;done

#suppression des fichiers vides
find ../final_cov2/ -size 0 -exec rm -f {} \;

#script R matrice
R CMD BATCH "--args $te" /home/duran/software/TRACKPOSON/Analyse_pipeline.R

###reformatage de la matrice
sed -e "s/\.txt//g" matrice_final.csv | sed -e "s/FALSE/0/g" | sed -e "s/TRUE/1/g" | sed -e "s/\-/\_/g" | sed -e "s/\./\_/g" | awk 'NR<2{print $0;next}{print $0| "sort -k1,1V"}'  > matrice_final_$te.csv


#script R traditionel + histo
R CMD BATCH "--args $te" /home/duran/software/TRACKPOSON/Analyse_tradi.R
