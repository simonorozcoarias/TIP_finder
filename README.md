# TIPfinder
A pipeline that aim to find TIPs activity from TE dynamics, using the methodology proposed by TRACKPOSON [1] and huge genomic datasets. TIPFinder works under HPC techniques, parallel programming and has the capacity of scalability over many computational nodes (or servers) and multi-core architectures, which make it especially functional for applications in massive sequencing projects that demands the current (post) genomic era. 

Installation:
We highly recommend to use and install Python packages into a Anaconda enviorenment. To create, execute the following command:

  conda create -n tip_finder
  
So, activate it

  conda activate tip_finder
  
Then install required python packages

  conda install -c anaconda mpi4py
  
  conda install -c anaconda psutil
  
for TIP_finder_utils

  conda install -c anaconda pandas 
  
  conda install -c conda-forge matplotlib
  
  conda install -c anaconda seaborn
  

Usage:

Previos Steps:
- makeblastdb -in reference_genome.fasta -dbtype nucl
- bedtools makewindows -g chr_list.txt -w 10000 > reference_genome_10kbwindows.bed

NOTE: the chr_list.txt must have following structure:
ChName1<TAB>length
ChName2<TAB>length
ChName3<TAB>length  

mpirun -np num_processes -hosts=$SLURM_JOB_NODELIST python3.8 TIP_finder.py -f file_reads.txt -o folder_results -t TE_family_name -b bowtie2_indexed_file -l blast_formated_reference_genome.fasta -w 10_kb_splitted_reference_genome.bed

References:

[1] Carpentier, M. C., Manfroi, E., Wei, F. J., Wu, H. P., Lasserre, E., Llauro, C., ... & Panaud, O. (2019). Retrotranspositional landscape of Asian rice revealed by 3000 genomes. Nature communications, 10(1), 1-12.
