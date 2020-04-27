"""
 TEIPfinder: Transposable Element Insertion Polymorphism finder
 © Copyright
 Developped by Simon Orozco Arias
 email: simon.orozco.arias@gmail.com
 2020
"""
import sys
import multiprocessing
import time
import os
import subprocess
import argparse

output = multiprocessing.Queue()
output2 = multiprocessing.Queue()

def createDicc(blastfile, id, output, init, end):
	#print("Initialing Thread "+str(id))
	readHits = {}
	fileTh = open(blastfile)
	for i, line in enumerate(fileTh):
		if i >= init and i < end:
			columns = line.split('\t')
			read = columns[3].replace('\n', '')
			chrs = columns[0]
			if read in readHits.keys():
				if chrs not in readHits[read]:
					readHits[read].append(chrs)
			else:
				chrlist = [chrs]
				readHits[read] = chrlist
		elif i >= end:
			break
	#print("Finishing Thread "+str(id))
	fileTh.close()
	output.put(readHits)


def parseBlastOutput(blastfile, DiccReadHits, id, output, init, end):
	#print("Initialing Thread "+str(id))
	partialResults = []
	fileTh = open(blastfile)
	for i, line in enumerate(fileTh):
		if i >= init and i < end:
			columns = line.split('\t')
			read = columns[3].replace('\n', '')
			if len(DiccReadHits[read]) < 2:
				if int(columns[2]) < int(columns[1]):
					partialResults.append(columns[0]+"\t"+columns[2]+"\t"+columns[1]+"\t"+read)
				else:
					partialResults.append(columns[0]+"\t"+columns[1]+"\t"+columns[2]+"\t"+read)
		elif i >= end:
			break
	#print("Finishing Thread "+str(id)+", size of queue: "+str(sys.getsizeof(output)))
	fileTh.close()
	output.put(partialResults)


if __name__ == '__main__':
	print("\n\n##############################################################\n\n")
	print("TEIPfinder: Transposable Element Insertion Polymorphism finder\n\n")
	print("##############################################################\n\n")
	### read parameters
	# usage: python3 TEIPfinder.py -f reads_files.txt -o test -t 12 -bw /data3/projects/arabica_ltr/dbs/retroTEs/Gypsy -bl /data3/projects/arabica_ltr/dbs/coffea_arabica_v0.6_06.25.19.fasta -w /data3/projects/arabica_ltr/dbs/coffea_arabica_v0.6_06.25.19.fasta.bed
	parser = argparse.ArgumentParser()
	parser.add_argument('-f','--reads-file',required=True,dest='readsFilePath',help='file with paths of reads file separated with commas with columns: nameOfSample,forwardReads,reverseReads')
	parser.add_argument('-o','--output-dir',required=True,dest='out',help='path of the output directory')
	parser.add_argument('-t','--threads',dest='threads',help='number of threads')
	parser.add_argument('-b','--db-bowtie',required=True,dest='DB',help='path to TE library bowtie2 index')
	parser.add_argument('-l','--db-blast',required=True,dest='blast_ref_database',help='path to blast reference database')
	parser.add_argument('-w','--windows',required=True,dest='win',help='path to reference genome 10kbp windows in bed format')
	(options,arguments) = parser.parse_args()
	readsFilePath = options.readsFilePath
	out = options.out
	threads = options.threads
	DB = options.DB
	blast_ref_database = options.blast_ref_database
	win = options.win
	if readsFilePath == None:
		print('Missing path of reads file. Exiting')
		sys.exit(0)
	if out == None:
		print('Missing path of output directory. Exiting')
		sys.exit(0)
	if threads == None:
		threads = 1
		print('Missing threads, using only 1')
	if DB == None:
		print('Missing path to TE library. Exiting')
		sys.exit(0)
	if blast_ref_database == None:
		print('Missing path to reference blast genome. Exiting')
		sys.exit(0)
	if win == None:
		print('Missing path to 10kb-splitted reference genome. Exiting')
		sys.exit(0)

	readsFile = open(readsFilePath, 'r').readlines()
	for sample in readsFile:
		te = sample.split(',')[0]
		fq1 = sample.split(',')[1]
		fq2 = sample.split(',')[2].replace('\n', '')

		output = multiprocessing.Queue()
		output2 = multiprocessing.Queue()
		
		### to map reads against TE library
		mappingCommand = "bowtie2 --time --end-to-end  -k 1 --very-fast -p "+str(threads)+" -x "+DB+" -1 "+fq1+" -2 "+fq2+" | samtools view -bS -@ 8 - > "+out+"-vs-"+te+".bam"
		subprocess.run(mappingCommand, shell=True)

		### keep only unmap reads with flag unmap/map
		keepingCommand = 'samtools view '+out+'-vs-'+te+'.bam | awk -F "\\t" \'{if ( ($1!~/^@/) && (($2==69) || ($2==133) || ($2==165) || ($2==181) || ($2==101) || ($2==117)) ) {print ">"$1"\\n"$10}}\' > '+out+'-vs-'+te+'.fa'
		subprocess.run(keepingCommand, shell=True)

		### blast against reference genome for identification insertion point
		blastCommand = 'blastn -db '+blast_ref_database+' -query '+out+'-vs-'+te+'.fa -out '+out+'-vs-'+te+'.fa.bl -outfmt "6 sseqid sstart send qseqid"  -num_threads '+str(threads)+' -evalue 1e-20'
		subprocess.run(blastCommand, shell=True)

		### process blast output to create a bed file
		countLinescommand = 'wc -l '+out+'-vs-'+te+'.fa.bl | cut -f1 -d" "'
		process = subprocess.run(countLinescommand, shell=True, stdout=subprocess.PIPE)
		fileLen = process.stdout
		start = time.time()

		blastfile = out+'-vs-'+te+'.fa.bl'
		outputfile = out+'-vs-'+te+'.bed'
		# execute in multiprocess mode
		processes = []
		lines_per_procs = int(int(fileLen)/int(threads))+1
		remain = int(fileLen) % int(threads)
		# Run processes
		for th in range(int(threads)):
			if th < remain:
				init = th * (lines_per_procs + 1)
				end = init + lines_per_procs + 1
			else:
				init = th * lines_per_procs + remain
				end = init + lines_per_procs
			p = multiprocessing.Process(target=createDicc, args=(blastfile, th, output, init, end))
			p.start()
			processes.append(p)

		print("all Process finished")
		# join all partial results in one
		results = [output.get() for p in processes]
		DiccReadHits = {}
		for dicc in results:
			for key in dicc.keys():
				if key in DiccReadHits.keys():
					for chrs in dicc[key]:
						if chrs not in DiccReadHits[key]:
							DiccReadHits[key].append(chrs)
				else:
					DiccReadHits[key] = dicc[key]
		finish = time.time()
		output = None
		results = None
		print("create dicc with reads done! time="+str(finish - start))

		print("dic size: "+str(sys.getsizeof(DiccReadHits)))

		# searching for reads with maximum 1 hits
		start = time.time()
		processes = []
		for th in range(int(threads)):
			if th < remain:
				init = th * (lines_per_procs + 1)
				end = init + lines_per_procs + 1
			else:
				init = th * lines_per_procs + remain
				end = init + lines_per_procs
			p = multiprocessing.Process(target=parseBlastOutput, args=(blastfile, DiccReadHits, th, output2, init, end))
			p.start()
			processes.append(p)

		# join all partial results in one
		results = [output2.get() for p in processes]
		#print(results)
		openoutputfile = open(outputfile, 'w')
		for partial in results:
			for line in partial:
				openoutputfile.write(line+"\n")
		finish = time.time()
		print("filter reads with one hit done! time="+str(finish - start))
		# closing files
		openoutputfile.close()
		###################################################################################

		### sort bed output
		sortCommand = 'sort -k1,1 -k2,2n '+out+'-vs-'+te+'.bed > '+out+'-vs-'+te+'.sort.bed'
		subprocess.run(sortCommand, shell=True)

		### coveragebed by 10kb windows
		coverageCommand = 'bedtools coverage -counts -nonamecheck -a '+win+' -b '+out+'-vs-'+te+'.sort.bed | awk -F "\\t" \'{if ($4>=2){print $0}}\' > coveragebed_'+out+'-vs-'+te+'_per10kb.bed'
		subprocess.run(coverageCommand, shell=True)

		### remove temporal elements
		os.remove(out+'-vs-'+te+'.bam')
		os.remove(out+'-vs-'+te+'.fa.bl')
		os.remove(out+'-vs-'+te+'.fa')
		os.remove(out+'-vs-'+te+'.bed')
