import sys

def contar(filetips, windows, maximum):
	copia = [0 for x in range(0, maximum, windows)]
	gypsy = [0 for x in range(0, maximum, windows)]
	tips = open(filetips, 'r').readlines()
	for line in tips:
		columns = line.split(',')
		end = int(columns[0].split('_')[2])
		interval = int(end/windows)
		copia[interval] += int(columns[1])
		gypsy[interval] += int(columns[2])
	
	outfile = open(filetips+'.window', 'w')
	for i in range(len(copia)):
		interval = (i+1)*windows
		outfile.write(str(interval)+','+str(copia[i])+','+str(gypsy[i])+'\n')


if __name__ == '__main__':
        TIPsCount = sys.argv[1]
        windows = int(sys.argv[2])
        maximum = 61000000
        contar(TIPsCount, windows, maximum)
