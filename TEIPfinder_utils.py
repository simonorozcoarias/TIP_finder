import sys
import os
import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
from scipy.stats import chi2


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
		plt.xlabel('Chromose length')		
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
	### read parameters
	util = sys.argv[1] # util to be used 
	te = sys.argv[2] # TE family name
	outputDir = sys.argv[3] # directory which will contain results
	
	### change value of these parameters before execute
	map_th = 5 # used in final Matrix util
	windows = 1000000 # used in peaks util
	prob = 0.95 # used in association analysis (confidence level)
	##############################################################

	print("##############################################################")
	print("TEIPfinder: Transposable Element Insertion Polymorphism finder")
	print("##############################################################")

	if util == "finalMatrix":
		directory = sys.argv[4] # directory with all coveraged files
		createFinalMatrix(te, directory, outputDir, map_th)
	elif util == "histograms":
		matrixPath = sys.argv[4]
		histograms(outputDir, matrixPath)
	elif util == "peaks":
		matrixPathCase1 = sys.argv[4]
		matrixPathCase2 = sys.argv[5]
		### to count number of TIPs by cases
		individuals = countPerChrs(matrixPathCase1, matrixPathCase2, outputDir, 1) # graph 1 indicates create graphs, otherwise indicates do not generate graphs.
		print(str(individuals)+" individuals processed")
		### using file created by countPerChrs, count number of TIPs per a given window size
		countPerWindow(outputDir+"/TIPscount_cases.csv", windows)
	elif util == "association":
		matrixPathCase1 = sys.argv[4]
		matrixPathCase2 = sys.argv[5]
		### to count number of TIPs by cases
		individuals = countPerChrs(matrixPathCase1, matrixPathCase2, outputDir, 0)
		print(str(individuals)+" individuals processed")
		associationTest(outputDir+"/TIPscount_cases.csv", individuals, prob)		
	else:
		print("Util "+util+" did not found")
		print("Util must be: finalMatrix or histograms or peaks or association")
	
	