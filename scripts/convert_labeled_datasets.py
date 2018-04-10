import math
import numpy as np
import sys
import os
import re
import time

# convert these unlabeled dataets to labeled
datasets = ['accidents','connect', 'retail', 'bms-pos', 'bms-web1', 'bms-web2','chess','T10I4D100K','T40I10D100K'];

k = 10


# scan the dataset to retrieve maximum item
for dat_name in datasets[:]:

	datfilename = "datasets/" + dat_name + "/"+ dat_name + ".dat"
	newdatfilename = "datasets/" + dat_name + "/"+ dat_name + "_new.dat"
	newlabelsfilename = "datasets/" + dat_name + "/"+ dat_name + "_new.labels"

	f_in = open(datfilename, 'r')

	ok = 1
	maxitem = 0
	numtrans = 0
	maxtranslength = 0
	while (ok==1) :
		line=f_in.readline()
		if len(line) > 0:

			numtrans = numtrans + 1

			items=line.split( )

			for i in items:
				if int(i)>maxitem:
					maxitem=int(i)
		else:
			ok = 0

	f_in.close()

	#write the spec file
	newspecfilename = "datasets/" + dat_name + "/"+ dat_name + "_new.spec"
	f_out_spec = open(newspecfilename,'w')
	f_out_spec.write(dat_name + "_new.dat" + "\n")
	f_out_spec.write(str(maxitem+1) + "\n")
	f_out_spec.write(str(maxtranslength+1) + "\n")
	f_out_spec.write(str(numtrans+1) + "\n")
	f_out_spec.write(dat_name + "_new.labels" + "\n")
	f_out_spec.close()
	print "spec file done "+newspecfilename


	f_in = open(datfilename, 'r')

	# we know maxitem, create array to count occurences of singletons
	counts = np.zeros([maxitem+1 , 2] , dtype=np.int)

	ok = 1
	while (ok==1) :
		line=f_in.readline()
		if len(line) > 0:

			items=line.split( )

			if maxtranslength < len(items):
				maxtranslength = len(items)

			for i in items:
				counts[int(i) , 0] = counts[int(i) , 0] + 1
				if counts[int(i) , 0] == 1:
					counts[int(i) , 1] = i
		else:
			ok = 0

	f_in.close()

	# sort items by support and print results

	counts = counts[np.argsort(-counts[:, 0])]

	print "num of trans of " + dat_name + " = " + str(numtrans)
	print "max trans length of " + dat_name + " = " + str(maxtranslength)
	print "items frequency of " + dat_name
	for i in counts[0:k,:]:
		print str(i)

	found = 0
	current_item=0
	label_item=-1
	label_item_support=-1
	min_supp_diff=numtrans
	midfreq=int(float(numtrans)/2.0)
	while found==0:
		if abs(counts[current_item,0]-midfreq) < min_supp_diff:
			min_supp_diff = abs(counts[current_item,0]-midfreq)
			label_item = counts[current_item,1]
			label_item_support = counts[current_item,0]
			current_item=current_item+1
		else:
			# i can stop here since all other items will not improve min_supp_diff
			found = 1

	# print selected item for label construction
	print "selected item for " + dat_name + " = " + str(label_item) + " support = " + str(label_item_support)


	#now create the new dataset and the class labels

	print "opening "+datfilename

	f_in = open(datfilename, 'r')
	f_out_trans = open(newdatfilename,'w')
	f_out_labels = open(newlabelsfilename,'w')

	ok = 1
	while (ok==1) :
		line=f_in.readline()
		if len(line) > 0:

			items=line.split( )
			line_to_write=""
			label_to_write=0
			# we use set to avoid duplicates in the same transaction
			items_set = set()
			for i in items:
				items_set.add(i)
			items_list = list(items_set)
			for i in items_list:
				if int(i)!=label_item:
					#write this item to the new dataset
					line_to_write=line_to_write + str(i) + " "
				else:
					#do not write this item, but assign the label to 1
					label_to_write =  1

			#write the transaction and the label
			if len(line_to_write) > 0:
				f_out_trans.write(line_to_write+"\n")
				f_out_labels.write(str(label_to_write)+"\n")
		else:
			ok = 0

	f_in.close()
	f_out_trans.close()
	f_out_labels.close()

	print "created successfully "+newdatfilename+" and "+newlabelsfilename
