"""
 TIP_finder: Transposable Element Insertion Polymorphisms finder
 Distributed version
 © Copyright
 Developed by Simon Orozco Arias
 email: simon.orozco.arias@gmail.com
 2020
"""

#!/bin/env python3
import sys
import time
import os
import subprocess
import argparse
from mpi4py import MPI
import pickle

def createDicc(blastfile, id, init, end):
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
	fileTh.close()
	return readHits


def parseBlastOutput(blastfile, DiccReadHits, id, init, end):
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
	fileTh.close()
	return partialResults


def receive_mpi_msg(src=MPI.ANY_SOURCE, t=MPI.ANY_TAG, deserialize=False):
	"""
	Sending MPI messages
	"""
	data_dict = {}
	status = MPI.Status()
	data = float(comm.recv(source=src, tag=t, status=status))
	if deserialize: 
		data = pickle.loads(data)
	data_dict['data'] = data
	data_dict['sender'] = int(status.Get_source())
	return data_dict


def send_mpi_msg(destination, data, serialize=False):
	"""
	Sending MPI messages
	"""
	if serialize: 
		data = pickle.dumps(data)
	comm.send(str(data), dest=destination)

if __name__ == '__main__':

	########################################################
	### Global variables
	global comm, rank, threads, rank_msg

	### Parallel process rank assignment
	comm = MPI.COMM_WORLD
	threads = comm.Get_size()
	rank = comm.Get_rank()
	rank_msg = '[Rank '+str(rank)+' msg]'

	if rank == 0:
		print("\n###################################################################")
		print("#                                                                 #")
		print("# TIP_finder: Transposable Element Insertion Polymorphisms finder #")
		print("#                   Distributed Version                           #")
		print("#                                                                 #")
		print("###################################################################\n")
	
	########################################################
	### read parameters
	# usage: mpirun -np 60 python3 TEIPfinder_distributed.py -f reads_files.txt -o test -t DEL -b /data3/projects/arabica_ltr/dbs/retroTEs/Gypsy -l /data3/projects/arabica_ltr/dbs/coffea_arabica_v0.6_06.25.19.fasta -w /data3/projects/arabica_ltr/dbs/coffea_arabica_v0.6_06.25.19.fasta.bed
	parser = argparse.ArgumentParser()
	parser.add_argument('-f','--reads-file',required=True,dest='readsFilePath',help='file with paths of reads file (in fastq format) separated with commas with columns: nameOfSample,PathForwardReads,PathReverseReads')
	parser.add_argument('-o','--output-dir',required=True,dest='out',help='path of the output directory')
	parser.add_argument('-b','--db-bowtie',required=True,dest='DB',help='path to TE library bowtie2 index')
	parser.add_argument('-l','--db-blast',required=True,dest='blast_ref_database',help='path to blast reference database')
	parser.add_argument('-w','--windows',required=True,dest='win',help='path to reference genome 10kbp windows in bed format')
	parser.add_argument('-t','--te',required=True,dest='te',help='name of the TE family')
	parser.add_argument('-v','--version',action='version', version='%(prog)s v1.0 Distributed (MPI version)')

	options = parser.parse_args()
	readsFilePath = options.readsFilePath
	out = options.out
	DB = options.DB
	blast_ref_database = options.blast_ref_database
	win = options.win
	te = options.te
	if readsFilePath == None:
		if rank == 0:
			print('Missing path of reads file. Exiting')
		sys.exit(0)
	if out == None:
		if rank == 0:
			print('Missing path of output directory. Exiting')
		sys.exit(0)
	if DB == None:
		if rank == 0:
			print('Missing path to TE library. Exiting')
		sys.exit(0)
	if blast_ref_database == None:
		if rank == 0:
			print('Missing path to reference blast genome. Exiting')
		sys.exit(0)
	if win == None:
		if rank == 0:
			print('Missing path to 10kb-splitted reference genome. Exiting')
		sys.exit(0)
	if te == None:
		if rank == 0:
			print('Missing name of the TE family. Exiting')
		sys.exit(0)

	########################################################
	### creating (if does not exist) the output directory
	if rank == 0:
		if os.path.exists(out):
			print("WARNING: output directory "+out+" already exist, be carefully")
		else:
			os.mkdir(out)

	readsFile = open(readsFilePath, 'r').readlines()
	for sample in readsFile:
		indname = sample.split(',')[0]
		fq1 = sample.split(',')[1]
		fq2 = sample.split(',')[2].replace('\n', '')

		extFile1 = os.path.splitext(fq1)[1]
		extFile2 = os.path.splitext(fq2)[1]

		if extFile1 not in [".fastq", ".fq"]:
			if rank == 0:
				print('FATAL ERROR processing '+indname+': forward read file does not have the correct extension (fastq or fq). Exiting')
			sys.exit(0)

		if extFile2 not in [".fastq", ".fq"]:
			if rank == 0:
				print('FATAL ERROR processing '+indname+': reverse read file does not have the correct extension (fastq or fq). Exiting')
			sys.exit(0)
		
		if threads < 2:
			if rank == 0:
				print('FATAL ERROR: Minimum number of process is 2, please check and re-run. Exiting')
			sys.exit(0)

		if rank == 0:
			print("Analyzing "+indname+"...")

		### count lines in forward reads
		countLinescommand = 'wc -l '+fq1+' | cut -f1 -d" "'
		process = subprocess.run(countLinescommand, shell=True, stdout=subprocess.PIPE)
		fileLen = int(process.stdout)

		### count lines in reverse reads
		countLinescommand = 'wc -l '+fq2+' | cut -f1 -d" "'
		process = subprocess.run(countLinescommand, shell=True, stdout=subprocess.PIPE)
		fileLen2 = int(process.stdout)

		if fileLen != fileLen2:
			if rank == 0:
				print("FATAL ERROR processing "+indname+": forward and reverse files have different lengths. Existing")
			sys.exit(0)

		### shared variables
		numReads = int(fileLen)/4
		reads_per_procs = int(int(numReads)/(threads - 1))+1
		remain = int(numReads) % threads

		blastfile = out+"/"+indname+'-vs-'+te+'.fa.bl'
		outputfile = out+"/"+indname+'-vs-'+te+'.bed'

		### master process
		if rank == 0:
			########################################################
			### Launch mapping, keeping unmaped reads and blast processes
			start = time.time()
			for th in range(1, threads):
				send_mpi_msg(th,1,serialize=True)
				# comm.send(1, dest=th)

			finalBlastFile = open(blastfile, 'w')
			for th in range(1, threads):
				data_dict = receive_mpi_msg(deserialize=True)
				data = data_dict['data']
				# data = comm.recv(source=th)
				thFile =  open(out+"/"+indname+'-vs-'+te+'.fa.bl.'+str(th), 'r').readlines()
				for line in thFile:
					finalBlastFile.write(line)
				thFile = None
				os.remove(out+"/"+indname+'-vs-'+te+'.fa.bl.'+str(th))
				os.remove(out+"/"+indname+"-vs-"+te+".bam."+str(th))
				os.remove(out+"/"+indname+'-vs-'+te+'.fa.'+str(th))
			finalBlastFile.close()
			finish = time.time()
			print("Mapping and blast done! time="+str(finish - start))

			########################################################
			### launch dic creation
			start = time.time()
			countLinescommand = 'wc -l '+out+"/"+indname+'-vs-'+te+'.fa.bl | cut -f1 -d" "'
			process = subprocess.run(countLinescommand, shell=True, stdout=subprocess.PIPE)
			fileLen = int(process.stdout)
			DiccReadHits = {}
			for th in range(1, threads):
				send_mpi_msg(th,fileLen,serialize=True)
				# comm.send(fileLen, dest=th)

			for th in range(1, threads):
				data_dict = receive_mpi_msg(deserialize=True)
				partialDic = data_dict['data']
				# partialDic = comm.recv(source=th)
				# join all partial results
				for key in partialDic.keys():
					if key in DiccReadHits.keys():
						for chrs in partialDic[key]:
							if chrs not in DiccReadHits[key]:
								DiccReadHits[key].append(chrs)
					else:
						DiccReadHits[key] = partialDic[key]

			finish = time.time()
			print("create dicc with reads done! time="+str(finish - start))

			########################################################
			### process blast output to create a bed file
			start = time.time()
			for th in range(1, threads):
				send_mpi_msg(th, DiccReadHits, serialize=True)
				# comm.send(DiccReadHits, dest=th)

			openoutputfile = open(outputfile, 'w')
			for th in range(1, threads):
				data_dict = receive_mpi_msg(deserialize=True)
				partialLines = data_dict['data']
				# partialLines = comm.recv(source=th)
				for line in partialLines:
					openoutputfile.write(line+"\n")
			openoutputfile.close()
			finish = time.time()
			print("filter reads with one hit done! time="+str(finish - start))

			########################################################
			### sort bed output
			start = time.time()
			sortCommand = 'sort -k1,1 -k2,2n '+out+"/"+indname+'-vs-'+te+'.bed > '+out+"/"+indname+'-vs-'+te+'.sort.bed'
			subprocess.run(sortCommand, shell=True)

			########################################################
			### coveragebed by 10kb windows
			coverageCommand = 'bedtools coverage -counts -nonamecheck -a '+win+' -b '+out+"/"+indname+'-vs-'+te+'.sort.bed | awk -F "\\t" \'{if ($4>=2){print $0}}\' > '+out+"/"+'coveragebed_'+indname+'-vs-'+te+'_per10kb.bed'
			subprocess.run(coverageCommand, shell=True)

			########################################################
			### remove temporal elements
			os.remove(out+"/"+indname+'-vs-'+te+'.fa.bl')
			os.remove(out+"/"+indname+'-vs-'+te+'.bed')
			finish = time.time()
			print("pos-processing TIPs done! time="+str(finish - start))
		
		### slave processes
		else:
			data_dict = receive_mpi_msg(deserialize=True)
			data = data_dict['data']
			# data = comm.recv(source=0)
			th = rank - 1
			if th < remain:
				init = th * (reads_per_procs + 1)
				end = init + reads_per_procs + 1
			else:
				init = th * reads_per_procs + remain
				end = init + reads_per_procs

			########################################################
			### to split reads in chuncks between initial and final line to process reads (in fastq format)
			lineInit = init * 4
			lineEnd = end * 4
			forward = open(out+"/"+indname+"_1.fastq."+str(rank), "w")
			reverse = open(out+"/"+indname+"_2.fastq."+str(rank), "w")
			readFiles = open(fq1)
			for i, line in enumerate(readFiles):
				if i >= lineInit and i < lineEnd:
					forward.write(line)
			forward.close()
			readFiles = open(fq2)
			for i, line in enumerate(readFiles):
				if i >= lineInit and i < lineEnd:
					reverse.write(line)
			reverse.close()

			########################################################
			### to map reads against TE library
			mappingCommand = "bowtie2 --time --end-to-end -k 1 --very-fast -p 1 -x "+DB+" -1 "+out+"/"+indname+"_1.fastq."+str(rank)+" -2 "+out+"/"+indname+"_2.fastq."+str(rank)+" | samtools view -bS -@ 8 - > "+out+"/"+indname+"-vs-"+te+".bam."+str(rank)
			subprocess.run(mappingCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

			########################################################
			### remove splitted read files
			os.remove(out+"/"+indname+"_1.fastq."+str(rank))
			os.remove(out+"/"+indname+"_2.fastq."+str(rank))
							
			########################################################
			### keep only unmap reads with flag unmap/map
			keepingCommand = 'samtools view '+out+"/"+indname+'-vs-'+te+'.bam.'+str(rank)+' | awk -F "\\t" \'{if ( ($1!~/^@/) && (($2==69) || ($2==133) || ($2==165) || ($2==181) || ($2==101) || ($2==117)) ) {print ">"$1"\\n"$10}}\' > '+out+"/"+indname+'-vs-'+te+'.fa.'+str(rank)
			subprocess.run(keepingCommand, shell=True)

			########################################################
			### blast against reference genome for identification insertion point
			blastCommand = 'blastn -db '+blast_ref_database+' -query '+out+"/"+indname+'-vs-'+te+'.fa.'+str(rank)+' -out '+out+"/"+indname+'-vs-'+te+'.fa.bl.'+str(rank)+' -outfmt "6 sseqid sstart send qseqid"  -num_threads 1 -evalue 1e-20'
			subprocess.run(blastCommand, shell=True)
			send_mpi_msg(0, 1, serialize=True)
			# comm.send(1, dest=0)

			########################################################
			### launch dic creation
			data_dict = receive_mpi_msg(deserialize=True)
			fileLen = data_dict['data']
			# fileLen = comm.recv(source=0)
			lines_per_procs = int(int(fileLen)/int(threads-1))+1
			remain = int(fileLen) % int(threads)
			th = rank - 1
			if th < remain:
				init = th * (lines_per_procs + 1)
				end = init + lines_per_procs + 1
			else:
				init = th * lines_per_procs + remain
				end = init + lines_per_procs
			partial = createDicc(blastfile, rank, init, end)
			send_mpi_msg(0, partial, serialize=True)
			# comm.send(partial, dest=0)
			
			########################################################
			# searching for reads with maximum 1 hits
			data_dict = receive_mpi_msg(deserialize=True)
			diccFromMaster = data_dict['data']
			# diccFromMaster = comm.recv(source=0)
			th = rank - 1
			if th < remain:
				init = th * (lines_per_procs + 1)
				end = init + lines_per_procs + 1
			else:
				init = th * lines_per_procs + remain
				end = init + lines_per_procs
			partial = parseBlastOutput(blastfile, diccFromMaster, rank, init, end)
			send_mpi_msg(0, partial, serialize=True)
			# comm.send(partial, dest=0)

		if rank == 0:
			print("Done!")
