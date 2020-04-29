"""
 TEIPfinder: Transposable Element Insertion Polymorphism finder Utilities
 © Copyright
 Developped by Simon Orozco Arias
 email: simon.orozco.arias@gmail.com
 2020
"""

#!/bin/env python3
import sys
import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
from scipy.stats import chi2
import argparse


def createFinalMatrix(te, directory, outputDir, map_th):
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
		if "coveragebed_" in covf and "-vs-"+te in covf:		
			file_number += 1
			lines = open(directory+"/"+covf, 'r').readlines()
			sample_name = covf.replace("coveragebed_", "").replace("-vs-"+te, "").replace("_per10kb.bed", "")
			sample_names.append(sample_name)
			for line in lines:
				columns = line.split('\t')
				if int(columns[3]) >= map_th:
					new_ins = columns[0]+"_"+columns[1]+"_"+columns[2]
					if new_ins in ins_samples.keys():
						ins_samples[new_ins].append(sample_name)
					else:
						ins_samples[new_ins] = [sample_name]
	
	final_matrix = pd.DataFrame(0, index=range(len(ins_samples.keys())), columns=range(len(sample_names)))
	final_matrix.columns = sample_names
	final_matrix.index = list(ins_samples.keys())
	for k in ins_samples.keys():
		for s in ins_samples[k]:
			final_matrix.loc[k, s] = 1
	
	final_matrix.to_csv(outputDir+'/final_matrix_'+te+'.csv')


def histograms(outputDir, matrixPath):
	final_matrix = pd.read_csv(matrixPath, index_col = 0) 
	final_matrix['sum'] = final_matrix.sum(axis=1)

	### plotting frequencies in log2 per each insertion site
	print("plotting frequencies in log base 2 per each insertion site")
	plt.figure(figsize=(9, 8))
	plt.xlabel('Frequencies (log)')
	sns.distplot([math.log(x,2) for x in final_matrix['sum']], color='b', hist_kws={'alpha': 0.4}, rug=True);
	plt.savefig(outputDir+'/insertion_frequencies_log.png', dpi=300)

	### plotting frequencies per each insertion site
	print("plotting frequencies per each insertion site")
	final_matrix = final_matrix.drop(columns="sum")
	final_matrix.loc['samplesSum'] = final_matrix.sum()
	plt.figure(figsize=(9, 8))
	sns.distplot(final_matrix.loc['samplesSum'], color='b', hist_kws={'alpha': 0.4}, rug=True);
	plt.xlabel('Frequencies')
	plt.savefig(outputDir+'/inserts_per_samples.png', dpi=300)

def countPerChrs(matrixPathCase1, matrixPathCase2, outputDir, graph):
	final_matrixcase1 = pd.read_csv(matrixPathCase1, index_col = 0)
	final_matrixcase2 = pd.read_csv(matrixPathCase2, index_col = 0)

	totalIndividuals = len(final_matrixcase1.columns) + len(final_matrixcase2.columns)

	### sum TIPs per case
	final_matrixcase1['sum'] = final_matrixcase1.sum(axis=1)
	final_matrixcase2['sum'] = final_matrixcase2.sum(axis=1)

	insertion_sites = list(final_matrixcase1.index.values)
	insertion_sites2 = list(final_matrixcase2.index.values)
	for ins in insertion_sites2:
		if ins not in insertion_sites:
			insertion_sites.append(ins)

	final_count = pd.DataFrame(0, index=range(len(insertion_sites)), columns=range(2))
	final_count.columns = ["Case 1", "Case 2"]
	final_count.index = insertion_sites
	for ins in insertion_sites:
		contCase1 = 0
		contCase2 = 0
		if ins in final_matrixcase1.index:
			contCase1 = final_matrixcase1.loc[ins, 'sum']
		if ins in final_matrixcase2.index:
			contCase2 = final_matrixcase2.loc[ins, 'sum']
		final_count.loc[ins, "Case 1"] = contCase1
		final_count.loc[ins, "Case 2"] = contCase2

	if graph == True:
		print("plotting frequencies of each case")
		plt.figure(figsize=(9, 8))
		sns.distplot(final_count['Case 1'], hist=False, rug=True, color='b');
		sns.distplot(final_count['Case 2'], hist=False, rug=True, color='r');
		plt.xlabel('Frequencies')
		plt.savefig(outputDir+'/inserts_per_case.png', dpi=300)
	final_count.to_csv(outputDir+'/TIPscount_cases.csv')
	return totalIndividuals - 2


def countPerWindow(filetips, windows):
	final_count = pd.read_csv(filetips, index_col = 0)
	### to find the biggest position of a TIP
	maximum = max([int(x.split("_")[2]) for x in list(final_count.index.values)])
	ymax1 = max(list(final_count["Case 1"])) 
	ymax2 = max(list(final_count["Case 2"])) 
	ymax = max(ymax1,ymax2)
	max([int(x.split("_")[2]) for x in list(final_count.index.values)])
	#print(maximum)
	### to extract uniq id of each chromosome
	chr_ids = list(dict.fromkeys([x.split("_")[0] for x in list(final_count.index.values)])) 

	for ch in chr_ids:
		case1 = [0 for x in range(0, int(maximum), int(windows))]
		case2 = [0 for x in range(0, int(maximum), int(windows))]
		for ins in list(final_count.index.values):
			if ch+"_" in ins:
				end = int(ins.split('_')[2])
				interval = int(end/windows)
				#print("interval "+str(interval)+", total Length: "+str(len(case1))+", insertion: "+ins)
				case1[interval] += int(final_count.loc[ins, "Case 1"])
				case2[interval] += int(final_count.loc[ins, "Case 2"])

		plt.figure(figsize=(9, 8))
		plt.plot([x for x in range(0, maximum, windows)], case1, 'b', [x for x in range(0, maximum, windows)], case2, 'r')
		plt.xlabel('Chromosome length')	
		axes = plt.gca()
		axes.set_ylim([0, ymax])
		plt.savefig(outputDir+'/chr_'+ch+'.png', dpi=300)
		plt.close()
		

def associationTest(filetips, individuals, prob):
	final_count = pd.read_csv(filetips, index_col = 0)
	assocations = pd.DataFrame(0, index=range(len(list(final_count.index))), columns=range(4))
	assocations.columns = ["Case 1", "Case 2", "P-value", "Association"]
	assocations.index = list(final_count.index.values)

	for ins in list(final_count.index.values):
		cases1Pos = int(final_count.loc[ins, "Case 1"])
		cases1Neg = individuals - int(final_count.loc[ins, "Case 1"])
		cases2Pos = int(final_count.loc[ins, "Case 2"])
		cases2Neg = individuals - int(final_count.loc[ins, "Case 2"])
		# contingency table
		table = [[cases1Pos,cases2Pos], [cases1Neg, cases2Neg]]
		stat, p, dof, expected = chi2_contingency(table)
		critical = chi2.ppf(prob, dof)
		# interpret p-value
		alpha = 1.0 - prob
		if p <= alpha:
			measureAsoc = math.sqrt(float(stat)/(float(stat) + individuals))
			assocations.loc[ins, "Case 1"] = cases1Pos
			assocations.loc[ins, "Case 2"] = cases2Pos
			assocations.loc[ins, "P-value"] = p
			assocations.loc[ins, "Association"] = measureAsoc

	### remove TIPs with no statistical assocation
	assocations.drop(assocations[assocations.Association == 0].index, inplace=True)
	print("TIPs with statistical assocation: "+str(len(assocations)))
	assocations.to_csv(outputDir+'/TIPS_with_association.csv')

if __name__ == '__main__':
	print("\n##################################################################################")
	print("#                                                                                #")
	print("# TEIPfinder Utils: Transposable Element Insertion Polymorphism finder utilities #")
	print("#                                                                                #")
	print("##################################################################################\n")

	### read parameters
	parser = argparse.ArgumentParser()
	parser.add_argument('-u','--util',required=True,dest='util',help='Utility to be used, must be: finalMatrix or histograms or peaks or association')
	parser.add_argument('-t','--te',dest='te',help='TE family name (used in finalMatrix util)')
	parser.add_argument('-o','--output-dir',required=True,dest='outputDir',help='Path of the output directory')
	parser.add_argument('-d','--directory',dest='directory',help='Directory which contains coveraged files. (used in finalMatrix util).')
	parser.add_argument('-m','--map-thr',dest='map_th',help='Minimum number of maps. Default 5. (used in finalMatrix util).')
	parser.add_argument('-f','--matrix-path',dest='matrixPath',help='Path to the final matrix. (used in histograms util).')
	parser.add_argument('-1','--case1-matrix',dest='matrixPathCase1',help='Path to the final matrix of case 1. (used in peaks and association utilities).')
	parser.add_argument('-2','--case2-matrix',dest='matrixPathCase2',help='Path to the final matrix of case 2. (used in peaks and association utilities).')
	parser.add_argument('-w','--windows',type=int,dest='windows',help='Window length to graph TIPs. Dafault 1000000. (used in peaks util)')
	parser.add_argument('-n','--confidence-level',type=float,dest='prob',help='confidence level used in association tests. Dafault 0.95. (used in association util)')
	parser.add_argument('-v','--version',action='version', version='%(prog)s v1.0')

	options = parser.parse_args()
	util = options.util
	te = options.te
	outputDir = options.outputDir
	directory = options.directory
	map_th = options.map_th
	matrixPath = options.matrixPath
	matrixPathCase1 = options.matrixPathCase1
	matrixPathCase2 = options.matrixPathCase2
	prob = options.prob
	windows = options.windows

	if util == None:
		print('Missing util parameter (-u or --util). Exiting')
		sys.exit(0)
	if outputDir == None:
		print('Missing output directory (-o or --output-dir). Exiting')
		sys.exit(0)
	
	########################################################
	### execution of utilities
	if util == "finalMatrix":
		if directory == None:
			print('Missing input directory (-d or --directory). Exiting')
			sys.exit(0)
		if te == None:
			print('Missing name of the TE family (-t or --te). Exiting')
			sys.exit(0)
		if map_th == None:
			print('Missing minimum number of maps (-m or --map-thr). Using by default 5')
			map_th = 5
		createFinalMatrix(te, directory, outputDir, int(map_th))
	elif util == "histograms":
		if matrixPath == None:
			print('Missing path to the final matrix (-f or --matrix-path). Exiting')
			sys.exit(0)
		histograms(outputDir, matrixPath)
	elif util == "peaks":
		if matrixPathCase1 == None:
			print('Missing path to the matrix of individuals from case 1 (-1 or --case1-matrix). Exiting')
			sys.exit(0)
		if matrixPathCase2 == None:
			print('Missing path to the matrix of individuals from case 2 (-2 or --case2-matrix). Exiting')
			sys.exit(0)
		if windows == None:
			print('Missing window length to be used in graphs (-w or --windows). Using by default 1000000')
			windows = 1000000
		### to count number of TIPs by cases
		individuals = countPerChrs(matrixPathCase1, matrixPathCase2, outputDir, 1) # graph 1 indicates create graphs, otherwise indicates do not generate graphs.
		print(str(individuals)+" individuals processed")
		### using file created by countPerChrs, count number of TIPs per a given window size
		countPerWindow(outputDir+"/TIPscount_cases.csv", int(windows))
	elif util == "association":
		if matrixPathCase1 == None:
			print('Missing path to the matrix of individuals from case 1 (-1 or --case1-matrix). Exiting')
			sys.exit(0)
		if matrixPathCase2 == None:
			print('Missing path to the matrix of individuals from case 2 (-2 or --case2-matrix). Exiting')
			sys.exit(0)
		if prob == None:
			print('Missing confidence level (-n or --confidence-level). Using by default 0.95')
			prob = 0.95
		### to count number of TIPs by cases
		individuals = countPerChrs(matrixPathCase1, matrixPathCase2, outputDir, 0)
		print(str(individuals)+" individuals processed")
		associationTest(outputDir+"/TIPscount_cases.csv", individuals, float(prob))		
	else:
		print("Util "+util+" did not found")
		print("Util must be: finalMatrix or histograms or peaks or association")
	
	
