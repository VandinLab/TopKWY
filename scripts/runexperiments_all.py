#!/usr/bin/python
import math
import numpy as np
import sys
import os
import re
import time
from random import randint

jp = 10000;
alpha = 0.05;
maxram = 100000;
flag = "final"
runs = 10

K_ = [10, 100, 1000, 10000, 100000, 1000000]

for K in K_:

	datasets = ['covtype','retail','mushroom','bms-web1','bms-web2','cod-rna',
	'breast-cancer','ijcnn1','a9a','svmguide3','T10I4D100K','chess','bms-pos',
	'T40I10D100K','phishing','connect','accidents','pumb-star','susy'];

	# store metrics associated to datasets
	times = dict()
	memory = dict()
	tested = dict()
	runs_count = dict()

	for i in datasets:
		times[i] = np.zeros(runs)
		memory[i] = np.zeros(runs)
		tested[i] = np.zeros(runs)
		runs_count[i] = 0


	def launch_experiment(dat_name , K , jp , alpha , flag , run_id):
		print "Run new experiment "+str(run_id)+" for " + dat_name
		print "K = "+str(K)+" jp = "+str(jp)+" alpha = "+str(alpha)+" flag = "+str(flag)
		# remove topkwy executable from dataset folder
		cmd = "rm datasets/" + dat_name + "/" + topkwy_name
		os.system(cmd)
		# copy the topkminer executable into the dataset folder
		cmd = "cp topkwy datasets/" + dat_name + "/" + topkwy_name
		os.system(cmd)
		# launch the experiment
		out_file_name = 'experiment_'+flag+'_run'+str(run_id)+"_"+str(alpha)+"_"+str(jp)+'_'+str(K)+'.txt';
		cmd = "cd datasets/" + dat_name + "/ && ./"+topkwy_name+" " + dat_name + "_new.spec "+str(maxram)+" "+str(K)+" "+str(jp)+" "+str(alpha)+" > "+out_file_name+" 2>&1 &"
		os.system(cmd)
		#print cmd


	topkwy_name = 'topkwy';

	if len(flag) > 0:
		out_file_name = 'experiment_'+flag+'_run'+str(0)+"_"+str(alpha)+"_"+str(jp)+'_'+str(K)+'.txt';
		topkwy_name = topkwy_name +'_'+ flag
	else:
		#out_file_name = 'experiment_'+str(jp)+'_'+str(K)+'.txt';
		out_file_name = 'experiment_'+flag+'_run'+str(0)+"_"+str(alpha)+"_"+str(jp)+'_'+str(K)+'.txt';

	# launch the first experiments
	for dat_name in datasets:
		launch_experiment(dat_name , K , jp , alpha , flag , 0)
		time.sleep(randint(10, 20))

	running_datasets = list()
	for i in datasets:
		running_datasets.append(i)

	running = 1
	while (running == 1):

		time.sleep( 10 )

		# check if experiment concluded and collect results
		for dat_name in running_datasets[:]:

			out_file_name = 'experiment_'+flag+'_run'+str(runs_count[dat_name])+"_"+str(alpha)+"_"+str(jp)+'_'+str(K)+'.txt';
			outname = "datasets/" + dat_name + "/"+out_file_name
			cmd = "tail -250 "+outname+" > tail_script.txt"
			os.system(cmd)
			#print cmd

			f_in = open('tail_script.txt', 'r')

			ok = 1
			found = 0
			while (ok==1) :
				line=f_in.readline()
				if len(line) > 0:

					searchObj = re.match( r'Total running time time:.*', line, re.M|re.I)
					if searchObj:
						time_txt = searchObj.group()
						time_txt = time_txt[25:]

						time_value = -1.0
						try:
							time_value = float(time_txt)
						except ValueError:
							time_value = -1.0

						print "  Running time: " + str(time_value)

						results = times[dat_name]
						results[runs_count[dat_name]] = time_value
						times[dat_name] = results

						print "  Tested: " + str(tested_value)

						results = tested[dat_name]
						results[runs_count[dat_name]] = tested_value
						tested[dat_name] = results

						running_datasets.remove(dat_name)
						found = 1

					searchObj = re.match( r'tested .*', line, re.M|re.I)
					if searchObj:
						tested_txt = searchObj.group()
						tested_txt = tested_txt[7:]
						tested_value = int(tested_txt)

					searchObj = re.match( r'Number of significant itemsets: .*', line, re.M|re.I)
					if searchObj:
						res_txt = searchObj.group()
						res_txt = res_txt[32:]
						res_value = int(res_txt)
						print dat_name + " finished run "+str(runs_count[dat_name]+1)
						print "  Results: " + str(res_value)

					searchObj = re.match( r'total peak memory usage: .*', line, re.M|re.I)
					if searchObj:
						memory_txt = searchObj.group()
						memory_txt = memory_txt[25:]
						memory_txt = memory_txt[:6]

						memory_value = -1.0
						try:
							memory_value = float(memory_txt)
						except ValueError:
							memory_value = -1.0

						results = memory[dat_name]
						results[runs_count[dat_name]] = memory_value
						memory[dat_name] = results


						print "  Memory usage: " + str(memory_value)
				else:
					ok = 0
			f_in.close()

		# check if all experiments concluded
		# this happens when all datasets performed all the runs
		if len(datasets) == 0:
			running = 0

		# experiments which are not running and not have done all the runs, start a new run
		for i in datasets:
			if i not in running_datasets:
				runs_count[i] = runs_count[i] + 1
				if runs_count[i] < runs:
					time.sleep(randint(10, 20))
					launch_experiment(i , K , jp , alpha , flag , runs_count[i])
					running_datasets.append(i)

					# print results of this run on file
					file_out_log = open(flag+"_"+str(K)+"_"+str(jp)+".txt",'a')
					file_out_log.write(str(K)+";"+str(jp)+";"+str(alpha)+";"+i+";"+str(times[i][runs_count[i]-1])+";"+str(memory[i][runs_count[i]-1])+";"+str(tested[i][runs_count[i]-1])+";"+str(runs_count[i]-1)+";\n")
					file_out_log.close()

				else:
					print " all runs done "
					# all runs concluded for this dataset
					# compute average of metrics and variance and print it to file
					b = times[i]
					times[i] = b[b>0]
					mean_time = np.mean(times[i])
					variance_time = np.var(times[i])

					b = memory[i]
					memory[i] = b[b>0]
					mean_memory = np.mean(memory[i])
					variance_memory = np.var(memory[i])

					all_results_file = open("all_results.csv",'a')
					print i+" done with time mean "+str(mean_time)+" and variance "+str(variance_time)
					print i+" done with memory mean "+str(mean_memory)+" and variance "+str(variance_memory)
					all_results_file.write(str(K)+";"+str(jp)+";"+str(alpha)+";"+i+";"+str(mean_time)+";"+str(variance_time)+";"+str(mean_memory)+";"+str(variance_memory)+";\n")
					datasets.remove(i)
					all_results_file.close()
