import sys
import multiprocessing
import time
import os
import gc

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
	blastfile = sys.argv[1]
	outputfile = sys.argv[2]
	threads = int(sys.argv[3])
	fileLen = int(sys.argv[4])
	start = time.time()
	# execute in multiprocess mode
	processes = []
	lines_per_procs = int(fileLen/threads)+1
	remain = fileLen % threads
	# Run processes
	for th in range(threads):
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
	for th in range(threads):
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
