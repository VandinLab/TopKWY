import math
import numpy as np
import sys
import os
import re
import time


dataset_name = sys.argv[1]
print dataset_name

adjacency_file = dataset_name+"_A.txt"
transactions_file = dataset_name+"_graph_indicator.txt"
glabels_file = dataset_name+"_graph_labels.txt"
vlabels_file = dataset_name+"_node_labels.txt"
elabels_file = dataset_name+"_edge_labels.txt"

edges = dict()
vertexes = dict()
labels = dict()


f_in = open(transactions_file, 'r')
f_in_2 = open(vlabels_file, 'r')

vertex_id = 1
ok=1
while ok==1:
	line = f_in.readline()
	line.replace('\n','')
	if len(line) > 0:
		#print line + " - " + str(len(line))
		t_id = int(line)
		vlabel = int(f_in_2.readline())
		if t_id in vertexes:
			t = vertexes[t_id]
			t.add( (vertex_id , vlabel) )
			vertexes[t_id] = t
		else:
			new_t = set()
			new_t.add( (vertex_id , vlabel) )
			vertexes[t_id] = new_t

		vertex_id += 1
	else:
		ok = 0

f_in.close()
f_in_2.close()

#print vertexes

f_in = open(glabels_file, 'r')

t_id = 1
ok=1
while ok==1:
	line = f_in.readline()
	line.replace('\n','')
	if len(line) > 0:
		#print line + " - " + str(len(line))
		t_label = int(line)
		if t_id in vertexes:
			labels[t_id] = t_label
		else:
			print "problem with transaction of id "+str(t_id)
			exit()

		t_id += 1
	else:
		ok = 0
f_in.close()

#print labels

if len(labels) != len(vertexes):
	print "problem with size of transactions!"
	exit()


f_in = open(adjacency_file, 'r')
edges_exists = 0
try:
	f_in_2 = open(elabels_file, 'r')
	edges_exists = 1
except IOError:
	edges_exists = 0


ok=1
while ok==1:
	line = f_in.readline()
	line.replace('\n','')
	if len(line) > 0:
		#print line + " - " + str(len(line))
		v_ids_ = line.split(',')
		v_ids_ = ( int(v_ids_[0]) , int(v_ids_[1]) )
		v_ids = ( min(v_ids_[0] , v_ids_[1]) , max(v_ids_[0] , v_ids_[1]) )
		e_label = 0
		if edges_exists == 1:
			e_label = int(f_in_2.readline())
		if v_ids[0] in edges:
			new_edges = edges[v_ids[0]]
			new_edges.add( (v_ids[1] , e_label) )
			edges[v_ids[0]] = new_edges
		else:
			new_edges = set()
			new_edges.add( (v_ids[1] , e_label) )
			edges[v_ids[0]] = new_edges

	else:
		ok = 0
f_in.close()

#print edges


# all info stored, now print the dataset in GASTON format
transactions_out_file = dataset_name+"_transactions.txt"
labels_out_file = dataset_name+"_labels.txt"
f_out = open(transactions_out_file, 'w')
f_out2 = open(labels_out_file, 'w')

for t_id in range(1 , len(vertexes)+1):

	#print transaction id
	#f_out.write("# "+str(t_id)+"\n")
	f_out.write("t # "+str(t_id-1)+"\n")
	#if labels[t_id] < 0 or labels[t_id] > 1:
	#	labels[t_id] = 0
	f_out2.write(str(labels[t_id])+"\n")

	vertex_map = dict()
	v_counter = 0

	transaction = vertexes[t_id]

	# check that all vertexes are connected
	connected_vertexes = set()
	for vertex_info in transaction:
		v_id = vertex_info[0]

		if v_id in edges:
			connected_vertexes.add(v_id)
			edges_of_v = edges[v_id]

			for edge in edges_of_v:
				connected_vertexes.add(edge[0])

	# print all vertexes

	for vertex_info in transaction:
		if vertex_info[0] in connected_vertexes:
			#print vertex_info
			f_out.write("v "+str(v_counter)+" "+str(vertex_info[1])+"\n")
			vertex_map[vertex_info[0]] = v_counter
			v_counter += 1

	# print all edges

	for vertex_info in transaction:
		v_id = vertex_info[0]

		if v_id in edges:
			edges_of_v = edges[v_id]

			for edge in edges_of_v:
				f_out.write("e "+str(vertex_map[v_id])+" "+str(vertex_map[edge[0]])+" "+str(edge[1])+"\n")

f_out.close()
