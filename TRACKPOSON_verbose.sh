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

start=`date +%s`
######mapping reads against TE reference
bowtie2 --time --end-to-end  -k 1 --very-fast -p $cores -x $DB  -1 $fq1 -2 $fq2  | samtools view -bS -@ 8 - > "$out"-vs-"$te".bam # Default Number_Threads -p 6
end=`date +%s`
elapsed1=`expr $end - $start`

#######keep only unmap reads with flag unmap/map
start=`date +%s`
#et fasta creation
samtools view "$out"-vs-"$te".bam | awk -F "\t" '{if ( ($1!~/^@/) && (($2==69) || ($2==133) || ($2==165) || ($2==181) || ($2==101) || ($2==117)) ) {print ">"$1"\n"$10}}' > $out-vs-$te.fa
end=`date +%s`
elapsed2=`expr $end - $start`

######blast fa against reference genome (IRGSP1.0) for identification insertion point
start=`date +%s`
blastn -db $blast_ref_database -query $out-vs-$te.fa -out $out-vs-$te.fa.bl -outfmt "6 sseqid sstart send qseqid"  -num_threads $cores -evalue 1e-20 #8 -evalue 1e-20
end=`date +%s`
elapsed3=`expr $end - $start`


######parse blast output to finde TE insertion point ##### filter reads that has more than 1 hit (a match with more than one database sequence)
start=`date +%s`
python ${splitter_script} $out-vs-$te.fa.bl $cores
#bedtools needs a file which format has startPos < endPos
python ${python_script} $out-vs-$te.fa.bl $out-vs-$te.bed $cores
######sort bed output
sort -k1,1 -k2,2n $out-vs-$te.bed > $out-vs-$te.sort.bed
end=`date +%s`
elapsed4=`expr $end - $start`

######coveragebed by 10kb windows
start=`date +%s`
bedtools coverage -counts -nonamecheck -a $win -b $out-vs-$te.sort.bed | awk -F "\t" '{if ($4>=2){print $0}}' > coveragebed_$out-vs-$te\_per10kb.bed
end=`date +%s`
elapsed5=`expr $end - $start`

echo "### execution time per steps;${elapsed1};${elapsed2};${elapsed3};${elapsed4};${elapsed5}"


######cleaning temporary files
rm $out-vs-$te.bam
rm $out-vs-$te.fa*
rm $out-vs-$te.bed
