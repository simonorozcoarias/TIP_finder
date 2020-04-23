import sys
import multiprocessing
import time
import os
import gc
from mpi4py import MPI

def createDicc(blastfile, id, init, end):
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
	return readHits


def parseBlastOutput(blastfile, DiccReadHits, id, init, end):
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
	return partialResults

if __name__ == '__main__':
	blastfile = sys.argv[1]
	outputfile = sys.argv[2]
	fileLen = int(sys.argv[3])

	### Global variables
	global comm, rank, size, rank_msg

	### Parallel process rank assignment
	comm = MPI.COMM_WORLD
	threads = comm.Get_size()
	rank = comm.Get_rank()
	rank_msg = '[Rank '+str(rank)+' msg]'
	print(rank_msg)
	lines_per_procs = int(fileLen/(threads - 1))+1
	remain = fileLen % threads

	# master process
	if rank == 0:
		print("starting master process")
		### launch dic creation
		DiccReadHits = {}
		for th in range(1, threads):
			comm.send(1, dest=th)

		for th in range(threads):
			partialDic = comm.recv(source=th)
			# join all partial results
			for key in partialDic.keys():
				if key in DiccReadHits.keys():
					for chrs in partialDic[key]:
						if chrs not in DiccReadHits[key]:
							DiccReadHits[key].append(chrs)
				else:
					DiccReadHits[key] = partialDic[key]
			print("Process "+str(th)+" finished")
		print("all Process finished")
		########################################################
		# searching for reads with maximum 1 hits
		print("starting master process in searching one hit")
		for th in range(1, threads):
			comm.send(DiccReadHits, dest=th)

		openoutputfile = open(outputfile, 'w')
		for th in range(threads):
			partialLines = comm.recv(source=th)
			for line in partialLines:
				openoutputfile.write(line+"\n")
			print("Process "+str(th)+" finished")
		openoutputfile.close()
		print("all Process finished")

	else: # slave processes
		### launch dic creation
		data = comm.recv(source=0)
		print("starting process "+str(rank)+" in dic creation")
		th = rank - 1
		if rank < remain:
			init = th * (lines_per_procs + 1)
			end = init + lines_per_procs + 1
		else:
			init = th * lines_per_procs + remain
			end = init + lines_per_procs
		partial = createDicc(blastfile, rank, init, end)
		comm.send(partial, dest=0)
		########################################################
		# searching for reads with maximum 1 hits
		diccFromMaster = comm.recv(source=0)
		print("starting process "+str(rank)+" in searching one hit")
		th = rank - 1
		if th < remain:
			init = th * (lines_per_procs + 1)
			end = init + lines_per_procs + 1
		else:
			init = th * lines_per_procs + remain
			end = init + lines_per_procs
		partial = parseBlastOutput(blastfile, diccFromMaster, rank, init, end)
		comm.send(partial, dest=0)
