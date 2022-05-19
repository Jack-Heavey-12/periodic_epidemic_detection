print('Starting Program!')

import networkx as nx
import sys
from gurobipy import *
import math
import numpy as np
import time
import pandas as pd
import random
import pickle


from collections import defaultdict

from matplotlib import pyplot as plt




def rounding(x, n, big_n):
	c = 2 + math.log(big_n, n+1)
	d = c * (math.log(n+1) ** 2)

	limit = np.random.uniform(0,1)

	#just a slightly cleaner way of doing the math, the other way should have still worked, though
	if limit <= np.min([1, x*d]):
		return 1
	return 0



def LinearProgram(graph, k, sources, num_samples=len(graph)):
	#graph is a list of cascade graphs or sampled sub_graphs
	#k is the defined budget of average number of tests per day
	
	### NOTE
	#sources is a list of nodes where sources[i] is the source for graph[i] 
	#Currently encoded to be an individual node, but might want to change that to take a list of nodes later?

	#num_samples can just be defaulted as length of graph, which is th number of samples. Don't really need to provide this, but did previously so *shrug*


	x = {} #Dict for defining variables, defined as in the paper
	y = {} #Dict for defining variables, defined as in the paper
	M = len(graph) #number of sample graphs in list graph
	nodes = graph[0].number_of_nodes()

	src = sources

	#could condense this, but want to make sure that it's the same list so it truly is random
	#could also be faster I suppose at the cost of memory
	lst_of_nodes = list(graph[0].nodes())


	m = Model('MinDelSS_Periodic')


	
	#Defines the V_id sets in a dictionary with first key i second key is d, which will return a list of nodes at that level.
	#The way this is constructed, 
	v_set = {}
	for i in range(M):
		empt_lst = lambda:[]
		v_set[i] = defaultdict(empt_lst)
		#Adds a meta node that connects to the source node. More useful for multiple sources, but ensures our correct distance definition.
		graph[i].add_node('meta')
		graph[i].add_edge('meta', src[i])
		#This will return 1 for source nodes and the distance to the rest of them
		paths = nx.shortest_path_length(graph[i], source='meta')
		for u in graph[i]:
			#allows for periodic testing of a period up to n
			for p in range(nodes):
				try:
					# paths[u] returns distance d. 
					dist = paths[u]
					# Defines d based on the function provided in the paper
					d = min(dist + (p - dist) % p, nodes + 1)
					v_set[i][d].append(u)
					#returns n+1 as the distance if there is an error (happens if there is not a path between the source and the current node).
				except KeyError:
					v_set[i][nodes + 1].append(u)
		#Only the meta node has a distance of 0, which is not a real node.
		del v_set[i][0]
		graph[i].remove_node('meta')


	for i in range(M):
		for d in v_set[i].keys():
			for u in v_set[i][d]:
				if u not in x:
					#adding all of the nodes to the x dictionary - only need to add each node once despite all of the samples.
					x[str(u)] = m.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'x['+str(u)+']')
					#x[str(i)+','+str(d)+','+str(u)] = m.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'x['+str(i)+','+str(d)+','+str(u)']')
			y[str(i)+','+str(d)] = m.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, ub = 1.0, name = 'y['+str(i)+','+str(d)+']')
	m.update()

	for i in range(M):
		for d in v_set[i].keys():
			m.addConstr(quicksum(x[str(u)] for u in v_set[i][d]) >= y[str(i)+','+str(d)], name='C1_sample='+str(i)+'_dist='+str(d))
	m.update()

	'''for i in range(M):
		for d in v_set[i].keys():
			m.addConstr(quicksum(x[str(u)] for u in v_set[i][d]) <= k, name='C2_sample='+str(i)+'_dist='+str(d))
	m.update()'''

	#in this instance, u is a tuple where the first value is the node label, the second value is the period as an integer
	#Thus, we pull out u[1]
	m.addConstr(quicksum(1/u[1] * x[str(u)] for u in list(x.keys())) <= k, name='C2_budget')


	for i in range(M):
		m.addConstr(quicksum(y[str(i)+','+str(d)] for d in v_set[i].keys()) == 1, name='C3_sample='+str(i))

	m.update()

	#Objective Function
	#First fraction is 1/N, need to somehow sum over all d values as well, think that is done in the comprehension right now but maybe need more clarity there
	m.setObjective(1/M * quicksum((y[str(i)+','+str(d)] * d) for i in range(M) for d in v_set[i].keys()), GRB.MINIMIZE)
	m.update()
	m.optimize()

	print('Model Optimized')



	#x is the dictionary where the values of x are stored, producing the rounding here
	n = len(x)
	x_prime_dict = {}
	for j in x.keys():
		x_prime_dict[j] = rounding(x[j].x, n, num_samples)


	#producing the number of infections here by saving the length of the connected component to the selected source
	infection_list = []
	for i in range(M):
		infection_list.append(len(nx.node_connected_component(graph[i], src[i])))



	#Want to calculate the objective value now for the rounded values
	#first, create a smaller set which has all of the nodes where they are in s_r
	s_r = set()
	for i in list(x_prime_dict.keys()):
		if x_prime_dict[i] == 1:
			s_r.add(i)

	degree_set = []

	node_lst = sorted(H.degree, key=lambda x: x[1], reverse=True)
	for i in range(int(sum(x_prime_dict.values()))):
		degree_set.append(node_lst[i][0])

	random_nodes = random.sample(list(H.nodes()), int(sum(x_prime_dict.values())))





	obj_vals = []
	for i in range(M):
		paths = nx.shortest_path_length(graph[i], source=src[i])
		shortest = nodes + 1
		for u in list(s_r):
			try:
				if paths[u]+1 < shortest:
					shortest = paths[u]
			except KeyError:
				continue
		obj_vals.append(shortest)

	print('New obj_vals found')

	obj_vals_degree = []
	for i in range(M):
		paths = nx.shortest_path_length(graph[i], source=src[i])
		shortest = nodes + 1
		for u in degree_set:
			try:
				if paths[u]+1 < shortest:
					shortest = paths[u]
			except KeyError:
				continue
		obj_vals_degree.append(shortest)

	obj_vals_random = []
	for i in range(M):
		paths = nx.shortest_path_length(graph[i], source=src[i])
		shortest = nodes + 1
		for u in random_nodes:
			try:
				if paths[u]+1 < shortest:
					shortest = paths[u]
			except KeyError:
				continue
		obj_vals_random.append(shortest)

	'''s_r_greedy = greedy(graph, int(sum(x_prime_dict.values())), src, v_set)

	obj_vals_greedy = []
	for i in range(M):
		#paths = nx.shortest_path_length(graph[i], source=src[i])
		shortest = nodes + 1
		for u in s_r_greedy:
			try:
				if paths[u]+1 < shortest:
					shortest = paths[u]
			except KeyError:
				continue
		obj_vals_greedy.append(shortest)

	print('Baselines Calculated')'''



	#Return stuff here, so want to calculate everything first so that we can return everything with just one function call


	#return np.mean(obj_vals), np.mean(obj_vals_degree), np.mean(obj_vals_random), np.mean(obj_vals_greedy), x_prime_dict, degree_set, random_nodes, #s_r_greedy
	return np.mean(obj_vals), np.mean(obj_vals_degree), np.mean(obj_vals_random), int(sum(x_prime_dict.values())), degree_set, random_nodes, #s_r_greedy

#Function for sampling edges with probability p. Will return a list of sample graphs.
def sampling(num_samples, graph, p):

	lst = []

	#Generates one graph per sample
	for i in range(num_samples):
		TempGraph = nx.Graph()

		TempGraph.add_nodes_from(graph.nodes())

		#Adds an edge if it is randomly selected with probability p
		for j in graph.edges():
			r = np.random.random()

			if r <= p:
				TempGraph.add_edge(j[0], j[1])

		lst.append(TempGraph)
		print(i)

	lst_of_nodes = list(graph.nodes())
	sources = []
	nodes = len(graph.nodes())

	for i in range(num_samples):
		val = np.random.randint(0, nodes-1)
		sources.append(lst_of_nodes[val])

	return lst, sources



def produce_plot(input_name, output_string):
	df = pd.read_csv(input_name)

	df['Percent'] = df['Rounded X Value'] / df['Value of K']

	fig = plt.plot(df['Value of K'], df['Value of K'], color='r', marker='o')

	plt.plot(df['Value of K'], df['Rounded X Value'], color='b', marker='o')
	plt.legend(['Value of Budget K', 'Size of Selected Set S_r'])
	plt.title('Violation of Budget Constraints')

	plt.savefig(output_string + '.png')
	plt.clf()

	fig2 = plt.plot(df['Value of K'], df['Percent'], 'o')
	plt.title('Violation of K as a Percentage of K')
	plt.savefig(output_string +'_percent.png')
	plt.clf()


	fig3 = plt.plot(df['Rounded X Value'], df['new_objVal'], color='m', marker='o')
	plt.title('Size of Sensor Set vs. Rounded Mean Detection Time')
	plt.savefig(output_string + '_objvK.png')
	plt.clf()


#main function
if __name__ == '__main__':
	#dataset = sys.argv[1] #input the original dataset/graph name
	outputfile = 'mindelss_output' #input the output file name
	
	#G = []
	paths = []
	levels = [] 
	sources = []
	inedges = []
 
	print("Generating Simulations")
	
	#Change num samples here
	#numSamples = 20
	
	#Change probability of transmission here
	#Used to determine the probability of a given edge being sampled
	p = 0.15

	#Change allowable infection rate here
	alpha = 0.3



	#So this graph is only 75 nodes, 1138 edges.
	#df = pd.read_csv('hospital_contacts', sep='\t', header=None)
	#df.columns = ['time', 'e1', 'e2', 'lab_1', 'lab_2']
	#H = nx.from_pandas_edgelist(df, 'e1', 'e2')

	#This graph is 2500 nodes, 9388 edges, but is randomly constructed
	#H = nx.read_adjlist('random_adj_list')

	#theoretical physics collaboration network
	#8638 nodes, 24827 edges
	#H_prime = nx.read_adjlist('arxiv_collab_network.txt')
	#H = H_prime.subgraph(max(nx.connected_components(H_prime))).copy()
	#del H_prime
	
	#UVA Hospital Network, 30 weeks from end. 9949 nodes, 399495 edges
	network = open('personnetwork_post', 'r')
	lines = network.readlines()
	lst = []
	for line in lines:
		lst.append(line.strip())
	network.close()
	H = nx.parse_edgelist(lst)
	del lst

	#UVA hospital Network, 30 weeks from beginning (skipping to at least time 10,000). 10789 nodes, 291881 edges
	#network = open('personnetwork_exp', 'r')
	#lines = network.readlines()
	#lst = []
	#for line in lines:
	#	lst.append(line.strip())
	#network.close()
	#H = nx.parse_edgelist(lst)
	#del lst

	#with open('Carilion_1hour.pkl', 'rb') as pkl:
	#	H = pickle.load(pkl)


	#nodes = len(H.nodes()) #small n

	k_arr = []
	obj_val_lst = []
	obj_val_degree_lst = []
	obj_val_random_lst = []
	obj_val_greedy_lst = []

	#prob_lst = [.001, .005, .01, .03, .05, .1, .15, .2, .25, .3, .35, .4]

	#num_samples = 3000
	#G, sources = sampling(num_samples, H, p)
	#paths = []
	#for i in G:
	#	paths.append(dict(nx.shortest_path_length(i)))

	G = []
	sources = []
	file = open('sources', 'r')
	lines = file.readlines()
	for line in lines:
		sources.append(line.strip())

	for i in range(3000):
		G.append(nx.read_edgelist('graphs/post_cov_p05_'+str(i)))


	#for i in range(25, 275, 25):
	for i in range(25,275,25):
	#for i in prob_lst:

		k = i #VALUE_FOR_K
		#k = np.floor((i/1000) * nodes)
		#p = i
		

		#UPDATE - Deleted LIST_OF_SOURCES as an argument, defined in function randomly.
		#obj_val_rounded, obj_val_degree, obj_val_random, obj_val_greedy, x_prime_dict, degree_set, random_nodes, s_r_greedy = LinearProgram_baselines(H, G, k, sources)
		obj_val_rounded, obj_val_degree, obj_val_random, x_prime_dict_sum, degree_set, random_nodes = LinearProgram_baselines(G, k, sources, 3000)


		#NEED TO EVALUATE THIS WRITE STATEMENT BELOW HERE
		k_arr.append(k)
		obj_val_lst.append(obj_val_rounded)
		obj_val_degree_lst.append(obj_val_degree)
		obj_val_random_lst.append(obj_val_random)
		obj_val_greedy_lst.append(len(random_nodes))





		textfile = open('final_output_speedy.csv', 'w')
		textfile.write('Value of K,Rounded Obj Value,Degree Obj Value,Random Obj Value,Greedy Obj Val\n')
		for val in range(len(k_arr)):
			#textfile.write(str(k_arr[val])+','+str(obj_val_lst[val])+','+str(obj_val_degree_lst[val])+
			#	','+str(obj_val_random_lst[val])+','+str(obj_val_greedy_lst[val])+'\n')
			textfile.write(str(k_arr[val])+','+str(obj_val_lst[val])+','+str(obj_val_degree_lst[val])+
				','+str(obj_val_random_lst[val])+str(obj_val_greedy_lst[val])+'\n')
		textfile.close()
		print(obj_val_greedy_lst)


	print('Done!')

	
