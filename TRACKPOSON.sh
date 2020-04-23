####You need to change the variables
##################################
fq1=$1 #file"1.fq.gz"
fq2=$2 #file"2.fq.gz"
out=$3 ##$3 #out_name
cores=$4
te=GYPSY #$4 #te_name
DIR=/data3/projects/arabica_ltr/TRACK2 #output_DIR
DB=/data3/projects/arabica_ltr/dbs/retroTEs/Gypsy #path_to_TE-bowtie2_index
blast_ref_database=/data3/projects/arabica_ltr/dbs/coffea_arabica_v0.6_06.25.19.fasta #$path_to_blast_ref_database
splitter_script=/data3/projects/arabica_ltr/TRACKPOSON/split_files.py
python_script=/data3/projects/arabica_ltr/TRACKPOSON/parse_blast_file.py #$path_to_find_insertion_point.pl
win="/data3/projects/arabica_ltr/dbs/coffea_arabica_v0.6_06.25.19.fasta.bed" #$path_to_ref_genome_10kbpwindows.bed
########################
#######################

cd $DIR

######mapping reads against TE reference

bowtie2 --time --end-to-end  -k 1 --very-fast -p $cores -x $DB  -1 $fq1 -2 $fq2  | samtools view -bS -@ 8 - > "$out"-vs-"$te".bam # Default Number_Threads -p 6

#######keep only unmap reads with flag unmap/map
#et fasta creation
samtools view "$out"-vs-"$te".bam | awk -F "\t" '{if ( ($1!~/^@/) && (($2==69) || ($2==133) || ($2==165) || ($2==181) || ($2==101) || ($2==117)) ) {print ">"$1"\n"$10}}' > $out-vs-$te.fa

######blast fa against reference genome (IRGSP1.0) for identification insertion point
blastn -db $blast_ref_database -query $out-vs-$te.fa -out $out-vs-$te.fa.bl -outfmt "6 sseqid sstart send qseqid"  -num_threads $cores -evalue 1e-20 #8 -evalue 1e-20
#blastn -db $blast_ref_database -query $out-vs-$te.fa -out $out-vs-$te.fa.bl -num_threads $cores -evalue 1e-20 #8 -evalue 1e-20

######parse blast output to finde TE insertion point ##### filter reads that has more than 1 hit (a match with more than one database sequence)

#python ${splitter_script} $out-vs-$te.fa.bl $cores
lines=`wc -l $out-vs-$te.fa.bl | cut -f1 -d" "`

#bedtools needs a file which format has startPos < endPos
python ${python_script} $out-vs-$te.fa.bl $out-vs-$te.bed $cores $lines

echo "Blast finished"

######sort bed output
sort -k1,1 -k2,2n $out-vs-$te.bed > $out-vs-$te.sort.bed

echo "Sort finished"

######coveragebed by 10kb windows
bedtools coverage -counts -nonamecheck -a $win -b $out-vs-$te.sort.bed | awk -F "\t" '{if ($4>=2){print $0}}' > coveragebed_$out-vs-$te\_per10kb.bed

echo "bedtools finished"
######cleaning temporary files
#rm $out-vs-$te.bam
#rm $out-vs-$te.fa*
#rm $out-vs-$te.bed
