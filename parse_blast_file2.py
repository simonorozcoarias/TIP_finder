import sys
import threading
import time

class DiccCreationWorker(threading.Thread):
	def __init__(self, idThread, lines):
		threading.Thread.__init__(self)
		self.idThread = idThread
		self.lines = lines
		self.readHits = {}

	def run(self):
		print("Initialing Thread "+str(self.idThread))
		for line in self.lines:
			columns = line.split('\t')
			read = columns[3].replace('\n', '')
			chrs = columns[0]
			if read in self.readHits.keys():
				if chrs not in self.readHits[read]:
					self.readHits[read].append(chrs)
			else:
				chrlist = [chrs]
				self.readHits[read] = chrlist
		print("Finishing Thread "+str(self.idThread))


class BlastParserWorker(threading.Thread):
	def __init__(self, idThread, lines, DiccReadHits):
		threading.Thread.__init__(self)
		self.idThread = idThread
		self.lines = lines
		self.readHits = {}
		self.DiccReadHits = DiccReadHits
		self.partialResults = []

	def run(self):
		print("Initialing Thread "+str(self.idThread))
		for line in self.lines:
			columns = line.split('\t')
			read = columns[3].replace('\n', '')
			if len(self.DiccReadHits[read]) < 2:
				if int(columns[2]) < int(columns[1]):
					self.partialResults.append(columns[0]+"\t"+columns[2]+"\t"+columns[1]+"\t"+read)
				else:
					self.partialResults.append(columns[0]+"\t"+columns[1]+"\t"+columns[2]+"\t"+read)
		print("Finishing Thread "+str(self.idThread))


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

	# execute in multithread mode
	start = time.time()
	lines_per_procs = int(len(openfile)/threads) + 1
	workers = []
	for x in range(threads):
		thread = DiccCreationWorker(x, openfile[(x*lines_per_procs + x): ((x*lines_per_procs + x) + lines_per_procs - 1)]) 
		thread.start()
		workers.append(thread)

	DiccReadHits = {}
	for worker in workers:
		worker.join()
		for key in worker.readHits.keys():
			if key in DiccReadHits.keys():
				for chrs in worker.readHits[key]:
					if chrs not in DiccReadHits[key]:
						DiccReadHits[key].append(chrs)
			else:
				DiccReadHits[key] = worker.readHits[key]

	finish = time.time()
	print("create dicc with reads done! time="+str(finish - start))

	# searching for reads with maximum 1 hits
	start = time.time()
	workers = []
	for x in range(threads):
		thread = BlastParserWorker(x, openfile[(x*lines_per_procs + x): ((x*lines_per_procs + x) + lines_per_procs - 1)], DiccReadHits) 
		thread.start()
		workers.append(thread)

	openoutputfile = open(outputfile, 'w')
	for worker in workers:
		worker.join()
		for line in worker.partialResults:
			openoutputfile.write(line+"\n")

	finish = time.time()
	print("filter reads with one hit done! time="+str(finish - start))
	openoutputfile.close()
