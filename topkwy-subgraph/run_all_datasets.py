import math
import numpy as np
import sys
import os
import re
import time


datasets = ["AIDS","BZR","COX2","DHFR","Mutagenicity","PTC_FM","PTC_FR","PTC_MM","PTC_MR","DD","ENZYMES","MUTAG","NCI1","NCI109","Tox21_AHR"]
bounded_1 = ["AIDS","BZR","COX2","DHFR","NCI1","NCI109","Tox21_AHR"]
bounded_2 = ["ENZYMES"]

ks = [10, 100, 1000, 10000, 100000, 1000000]
jp = 10000
runs = 10
small_wait = 5
medium_wait = 180
long_wait = 600
execute_wy_light = 0
execute_unbounded = 0

runnning_processes = 0
max_running = 25

bounds = dict()
for d in datasets:
    bounds[d] = -1
for d in bounded_1:
    bounds[d] = 20
for d in bounded_2:
    bounds[d] = 8

def get_file_suffix(light , j , dataset , k , run_id):

    suffix_ = ""
    if light == 1:
        suffix_ = suffix_ + "light_"

    suffix_ = suffix_ + str(dataset)+"_j10e"+str(int(math.log10(j)))

    if k > 0:
        suffix_ = suffix_+"_k10e"+str(int(math.log10(k)))
    if bounds[dataset] > 0:
        suffix_ = suffix_+"_m"+str(bounds[dataset])

    suffix_ = suffix_+"_r"+str(run_id)+".txt"

    return suffix_

def get_cmd(light , dataset , j , k, run_id):

    if light == 1:
        cmd_ = "./sgmine_light"
    else:
        cmd_ = "./sgmine"

    cmd_ = cmd_+ " -w"
    if bounds[dataset] > 0:
        cmd_ = cmd_+ " -m "+str(bounds[dataset])
    if k > 0:
        cmd_ = cmd_+ " -k "+str(k)

    suff = get_file_suffix(light , j , dataset , k , run_id)

    cmd_ = cmd_+ " -j "+str(j)+" -i datasets/"+str(dataset)+"/"+str(dataset)+"_transactions.txt -c datasets/"+str(dataset)+"/"+str(dataset)+"_labels.txt"
    cmd_ = cmd_+ " -o exp_output_"+str(suff)

    cmd_ = cmd_+ " > stats_"+str(suff)+" 2>&1 &"

    return cmd_


def check_results(light , dataset , j , k, run_id):

    suff = get_file_suffix(light , j , dataset , k , run_id)
    stats_file_path = "stats_"+str(suff)
    cmd = "tail -50 "+stats_file_path+" > temp_stats.txt"
    os.system(cmd)

    #store results
    results = -1
    running_time = 0.0
    wy_time = 0.0
    peak_memory = 0

    f_in = open("temp_stats.txt" , 'r')
    for line in f_in:

        #print line.replace('\n','')

        searchObj = re.match( r'Number of significant subgraphs: .*', line, re.M|re.I)
        if searchObj:
            temp_ = searchObj.group()
            temp_ = temp_[len("Number of significant subgraphs: "):]
            results = int(float(temp_))
            print str(dataset)+" finished run " + str(run_id_) + " "+str(stats_file_path)
            print "  Results : " + str(results)

        searchObj = re.match( r'Runtime for correction .*', line, re.M|re.I)
        if searchObj:
            temp_ = searchObj.group()
            temp_ = temp_[len("Runtime for correction (s): "):]
            running_time = float(temp_)
            print "  Running time : " + str(running_time)

        searchObj = re.match( r'Runtime for WY computations .*', line, re.M|re.I)
        if searchObj:
            temp_ = searchObj.group()
            temp_ = temp_[len("Runtime for WY computations (s): "):]
            wy_time = float(temp_)
            print "  Running time (WY) : " + str(wy_time)

        searchObj = re.match( r'Peak Memory.*', line, re.M|re.I)
        if searchObj:
            temp_ = searchObj.group()
            temp_ = temp_[len("Peak Memory (MB): "):]
            peak_memory = int(float(temp_))
            print "  Memory : " + str(peak_memory)

    if results > 0:
        #print results for this run to file
        f_out = open("all_runs.csv" , 'a')
        f_out.write(str(light)+";"+str(k)+";"+str(j)+";"+str(dataset)+";"+str(run_id)+";"+str(results)+";"+str(running_time)+";"+str(wy_time)+";"+str(peak_memory)+";\n")
        f_out.close()
        return 1
    return 0


# main script

flag = 1
if len(sys.argv) > 1:
    flag = int(sys.argv[1])

to_execute = list()
running_cmds = set()

for run_id_ in range(runs):
    for dataset in datasets:

        if flag == 1:
            # experiment for wy light
            if execute_wy_light == 1:
                cmd = get_cmd(1 , dataset , jp , 0 , run_id_)
                to_execute.append(cmd)
                #os.system(cmd)

            if execute_unbounded == 1:
                # experiment for topkwy with k=inf
                cmd = get_cmd(0 , dataset , jp , 0 , run_id_)
                to_execute.append(cmd)
                #os.system(cmd)

            for k in ks:
                # experiment for topkwy and k < inf
                cmd = get_cmd(0 , dataset , jp , k , run_id_)
                to_execute.append(cmd)
                #os.system(cmd)

    while len(to_execute) > 0 and len(running_cmds) < max_running:
        cmd = to_execute.pop()
        os.system(cmd)
        print cmd
        running_cmds.add(cmd)
        print len(running_cmds)
        time.sleep(np.random.randint(small_wait,2*small_wait))

    while len(running_cmds) > 0:

        # check all datasets before next run
        for dataset in datasets:

            # check result for wy light
            cmd = get_cmd(1 , dataset , jp , 0 , run_id_)
            if cmd in running_cmds:
                if check_results(1 , dataset , jp , 0 , run_id_) == 1:
                    running_cmds.remove(cmd)
                time.sleep(np.random.randint(small_wait,2*small_wait))

            # check result for topkwy with k=inf
            cmd = get_cmd(0 , dataset , jp , 0 , run_id_)
            if cmd in running_cmds:
                if check_results(0 , dataset , jp , 0 , run_id_) == 1:
                    running_cmds.remove(cmd)
                time.sleep(np.random.randint(small_wait,2*small_wait))

            for k in ks:
                # check result for topkwy with k < inf
                cmd = get_cmd(0 , dataset , jp , k , run_id_)
                if cmd in running_cmds:
                    if check_results(0 , dataset , jp , k , run_id_) == 1:
                        running_cmds.remove(cmd)
                    time.sleep(np.random.randint(small_wait,2*small_wait))

        while len(to_execute) > 0 and len(running_cmds) < max_running:
            cmd = to_execute.pop()
            os.system(cmd)
            running_cmds.add(cmd)
            print len(running_cmds)
            print cmd
            time.sleep(np.random.randint(small_wait,2*small_wait))

        time.sleep(np.random.randint(medium_wait,2*medium_wait))
