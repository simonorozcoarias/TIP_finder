def contar(fileenfermos, filesanos):
	enfermos = open(fileenfermos, 'r')
	sanos = open(filesanos, 'r')

	sumaEnfermos = {}
	linesEnf = enfermos.readlines()
	i = 0
	for line in linesEnf:
        	if i > 0:
                	line = line.replace('\n', '')
                	columns = line.split('\t')
                	columnsInt = [int(x) for x in columns[1:]]
                	sumaEnfermos[columns[0]] = sum(columnsInt)
        	i += 1

	sumaSanos = {}
	linesSan = sanos.readlines()
	i = 0
	for line in linesSan:
        	if i > 0:
                	line = line.replace('\n', '')
                	columns = line.split('\t')
                	columnsInt = [int(x) for x in columns[1:]]
               		sumaSanos[columns[0]] = sum(columnsInt)
        	i += 1

	for key in sumaEnfermos.keys():
        	tipssanos = 0
        	if key in sumaSanos.keys():
                	tipssanos = sumaSanos[key]
        	print key+","+str(sumaEnfermos[key])+","+str(tipssanos)

if __name__ == '__main__':
        enfermos = 'matrice_final_COPIA_onlyChrs.csv'
        sanos = 'matrice_final_GYPSY_onlyChrs.csv'
        contar(enfermos, sanos)
