#!/bin/bash

# to sum all TIPs from species
python contarTIPsXChr.py > TIPscount.csv

# to extract TIPs per chromsome
for seq in `cut -f1 -d"," TIPscount.csv | cut -d"_" -f1 | sort -u`; 
do 
  grep "${seq}_" TIPscount.csv > TIPscount_${seq}.csv; 
done

# to sum number of TIPs for windows 
for file in `ls TIPscount_*.csv`; do 
  echo $file; 
  python contarTIPsXVentana.py $file 1000000; 
done
