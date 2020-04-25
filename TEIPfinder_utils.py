import sys
import os
import pandas as pd

if __name__ == '__main__':
	### read parameters
	te = sys.argv[1] # TE family name
	directory = sys.argv[2] # directory with all coveraged files
	outputDir = sys.argv[3] # directory which will contain results

	### change value of these parameters before execute
	map_th = 5
	##############################################################
	
	if os.path.exists(outputDir):
		print("output directory: "+outputDir+" already exist, be carefully")
	else:
		os.mkdir(outputDir)

	##############################################################

	### creating final matrix
	sample_names = []
	ins_samples = {}
	file_number = 0
	for covf in os.listdir(directory):
		if "coveragebed_"+te in covf:
			file_number += 1
			lines = open(directory+"/"+covf, 'r').readlines()
			sample_name = covf.replace("coveragebed_"+te+"-vs-", "").replace("_per10kb.bed", "")
			sample_names.append(sample_name)
			for line in lines:
				columns = line.split('\t')
				if int(columns[3]) >= map_th:
					new_ins = columns[0]+"_"+columns[1]+"_"+columns[2]
					if new_ins in ins_samples.keys():
						ins_samples[new_ins].append(sample_name)
					else:
						ins_samples[new_ins] = [sample_name]
	
	columns_names = "insertion,"+",".join(sample_names)
	final_matrix = pd.DataFrame(0, index=range(len(ins_samples.keys())), columns=range(len(sample_names)))
	final_matrix.columns = sample_names
	final_matrix.index = list(ins_samples.keys())
	print(final_matrix)
	for k in ins_samples.keys():
		print(k)
		for s in ins_samples[k]:
			print(s)
			final_matrix[k][s] = 1
	print(final_matrix)
