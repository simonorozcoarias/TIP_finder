import sys
import multiprocessing
import os

def createDicc(blastfile, id):
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
	fileTh.close()
	return readHits


def parseBlastOutput(blastfile, DiccReadHits, id):
        fileTh = open(blastfile+"."+str(id),'r')
	lines = fileTh.readlines()
	partialResults = []
	for line in lines:
		columns = line.split('\t')
		read = columns[3].replace('\n', '')
		if len(DiccReadHits[read]) < 2:
			if int(columns[2]) < int(columns[1]):
				partialResults.append(columns[0]+"\t"+columns[2]+"\t"+columns[1]+"\t"+read)
			else:
				partialResults.append(columns[0]+"\t"+columns[1]+"\t"+columns[2]+"\t"+read)
	os.remove(blastfile+"."+str(id))
	return partialResults

if __name__ == '__main__':
	blastfile = sys.argv[1]
	outputfile = sys.argv[2]
	threads = int(sys.argv[3])

	# calculating how many lines will be processed by each thread
	# execute in multiprocess mode
	pool = multiprocessing.Pool(processes=threads)
	localresults = [pool.apply_async(createDicc, args=(blastfile, x)) for x in range(threads)]
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

	# searching for reads with maximum 1 hits
	pool = multiprocessing.Pool(processes=threads)
	localresults = [pool.apply_async(parseBlastOutput, args=(blastfile, DiccReadHits, x)) for x in range(threads)]
	results = [p.get() for p in localresults]
	openoutputfile = open(outputfile, 'w')
	for partial in results:
		for line in partial:
			openoutputfile.write(line+"\n")
	# closing files
	openoutputfile.close()
