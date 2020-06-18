# TIP_finder
A pipeline that aim to find TIPs activity from TE dynamics, using the methodology proposed by TRACKPOSON [1] accelerating the execution time up to 55 times in huge genomic datasets. TIP_finder applies a parallel strategy to work under HPC techniques efficiently, and has the capacity of scalability over many computational nodes (or servers) and multi-core architectures, which make it especially functional for applications in massive sequencing projects that demands the current (post) genomic era. 

### Prerequisites
TIP_finder used following bioinformatic software: Bowtie2 (v. 2.3.4.1) in order to map the paired reads of genomic data onto indexed consensus sequence of each TEs/HERVs family, Samtools (v. 1.9) to process bowtie2 output and to keep only unmapped reads, bedtools (v. 2.26.0) to split the reference genome into 10kb windows and to count reads in these windows, and NCBI-blastn (v. 2.10.0) or Magic-BLAST to align the unmapped reads.

TIP_finder was developed using Python3 (3.8) and following libraries: sys, time, os, subprocess, argparse, MPI4py. In the other hand, TIP_finder_utils.py requires the additional libraries: math, Pandas, matplotlib, Seaborn and SciPy.

## Installation:
We highly recommend to use and install Python packages into a Anaconda enviorenment. To create, execute the following command:
```
conda create -n tip_finder
```
So, activate it
```
conda activate tip_finder
```
Then install required python packages
```
conda install -c anaconda mpi4py
conda install -c anaconda psutil
conda install -c bioconda magicblast
```
for TIP_finder_utils
```
conda install -c anaconda pandas 
conda install -c conda-forge matplotlib
conda install -c anaconda seaborn
```
## Usage:

### Previos Steps:
- bowtie2-build TYPE_ref_retrotes.fa TYPE_ref_retrotes
- makeblastdb -in reference_genome.fasta -dbtype nucl (if you are using magicblast, The -parse_seqids option is required)
- bedtools makewindows -g chr_list.txt -w 10000 > reference_genome_10kbwindows.bed
- to create the comma-separated read_files.txt, which contains three columns: 1) datasets name, 2) path to the forward-reads file, 3) path to the reverse-reads file. 

### NOTE: 

- the chr_list.txt must have following structure (separated by tabs):
```
Ch1Name<TAB>length
Ch2Name<TAB>length
Ch3Name<TAB>length  
```
- the read_files.txt musth have following structure (separated by commas):
```
dataset1Name,path_to_forward_reads.fastq,path_to_reverse_reads.fastq
dataset2Name,path_to_forward_reads.fastq,path_to_reverse_reads.fastq
dataset3Name,path_to_forward_reads.fastq,path_to_reverse_reads.fastq
dataset4Name,path_to_forward_reads.fastq,path_to_reverse_reads.fastq
```
### Execution
```
mpirun -np num_processes -hosts=server_name python3.8 TIP_finder.py -f file_reads.txt -o folder_results -t TE_family_name -b TYPE_ref_retrotes -l reference_genome.fasta -w reference_genome_10kbwindows.bed
```
Where num_processes are the number of processors available in yout system and server_name is the name of the server where TIP_finder will run.

### NOTE
If you want to run TIP_finder using a SLURM job, to can use following script:
```
#!/bin/bash

#SBATCH --job-name=TIP_finder
#SBATCH -D /path/to/your/working/directory
#SBATCH --output=output_file.out
#SBATCH -e error_file.err
#SBATCH -n number_of_processors
#SBATCH -N 1

# remeber to load the prerequesities such as bowtie2, samtools, etc (for example using module load if it exists in your system).

conda activate tip_finder

echo "Running in $SLURM_JOB_NODELIST"
mpirun -np number_of_processors -hosts=$SLURM_JOB_NODELIST python3.8 TIP_finder.py -f file_reads.txt -o folder_results -t TE_family_name -b TYPE_ref_retrotes -l reference_genome.fasta -w reference_genome_10kbwindows.bed
```
### Help
if you need more information about how to run TIP_finder please execute:
```
python3 TIP_finder.py -h
```

## References:

[1] Carpentier, M. C., Manfroi, E., Wei, F. J., Wu, H. P., Lasserre, E., Llauro, C., ... & Panaud, O. (2019). Retrotranspositional landscape of Asian rice revealed by 3000 genomes. Nature communications, 10(1), 1-12.
