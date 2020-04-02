import sys

if __name__ == '__main__':
        blastfile = sys.argv[1]
        threads = int(sys.argv[2])

        # reading input blast file
        blastf = open(blastfile, 'r')
        openfile = blastf.readlines()

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

