import sys
import multiprocessing
import time
import os

def createDicc(blastf, id):
	print("Initialing Thread "+str(id))
	fileTh = open(blastfile+"."+str(id),'r')
	lines = fileTh.readlines()
	readHits = {}
	for line in lines:
		columns = line.split('\t')
		read = columns[3].replace('\n', '')
		chrs = columns[0]
		if read in readHits.keys():
			if chrs not in readHits[read]:
				readHits[read].append(chrs)
		else:
			chrlist = [chrs]
			readHits[read] = chrlist
	print("Finishing Thread "+str(id))
	fileTh.close()
	os.remove(blastfile+"."+str(id))
	return readHits


def parseBlastOutput(lines, DiccReadHits, lines_per_procs, x):
	print("Initialing Thread "+str(id))
	partialResults = []
	for line in lines:
		columns = line.split('\t')
		read = columns[3].replace('\n', '')
		if len(DiccReadHits[read]) < 2:
			if int(columns[2]) < int(columns[1]):
				partialResults.append(columns[0]+"\t"+columns[2]+"\t"+columns[1]+"\t"+read)
			else:
				partialResults.append(columns[0]+"\t"+columns[1]+"\t"+columns[2]+"\t"+read)
	print("Finishing Thread "+str(id))
	return partialResults

if __name__ == '__main__':
	blastfile = sys.argv[1]
	outputfile = sys.argv[2]
	threads = int(sys.argv[3])

	start = time.time()
	# reading input blast file
	blastf = open(blastfile, 'r')
	openfile = blastf.readlines()
	finish = time.time()
	print("Read blast file done! time="+str(finish - start))

	# splitting original blast file in chuncks for each thread
	start = time.time()
	lines_per_procs = int(len(openfile)/threads)+1
	for th in range(threads):
		init = th*lines_per_procs + th
		end = (th*lines_per_procs + th) + lines_per_procs - 1
		fileTh = open(blastfile+"."+str(th),'w')
		for line in openfile[init:end]:
			fileTh.write(line)
		fileTh.close()
	blastf.close()
	openfile = None
	finish = time.time()
	print("Split file in chuncks done! time="+str(finish - start))

	# calculating how many lines will be processed by each thread
	start = time.time()
	# execute in multiprocess mode
	pool = multiprocessing.Pool(processes=threads)
	localresults = [pool.apply_async(createDicc, args=(blastf, x)) for x in range(threads)]
	# join all partial results in one
	results = [p.get() for p in localresults]
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
	print("create dicc with reads done! time="+str(finish - start))

	# searching for reads with maximum 1 hits
	start = time.time()
	pool = multiprocessing.Pool(processes=threads)
	localresults = [pool.apply_async(parseBlastOutput, args=(openfile[(x*lines_per_procs + x): ((x*lines_per_procs + x) + lines_per_procs - 1)], DiccReadHits, lines_per_procs, x)) for x in range(threads)]
	results = [p.get() for p in localresults]
	openoutputfile = open(outputfile, 'w')
	for partial in results:
		for line in partial:
			openoutputfile.write(line+"\n")
	finish = time.time()
	print("filter reads with one hit done! time="+str(finish - start))
	# closing files
	openoutputfile.close()
