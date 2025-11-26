import sys
import os
import argparse
import time
from math import exp, log, fabs
import random
from treeswift import *
import treeswift
import json
import numpy as np
import cvxpy as cp

def __label_tree__(tree_obj):
	is_labeled = True
	i = 0
	labels = set()
	for node in tree_obj.traverse_preorder():
		if not node.is_root():
			if not node.edge_length:
				node.edge_length = 0
		if node.is_leaf():
			continue
		if not node.label or node.label in labels or node.label[0] != 'I': 
			is_labeled = False
			node.label = 'I' + str(i)
			i += 1        
		labels.add(node.label)
	return is_labeled

def compute_all_pairwise_LCAs(tree):
	lca_dict = {}
	node_set = {}
	for node in tree.traverse_postorder():
		node_set[node.label] = [node.label]
		if not node.is_leaf():
			for i in range(len(node.child_nodes())):
				c1 = node.child_nodes()[i]
				c1_node_set = node_set[c1.label]
				for l in c1_node_set:
					lca_dict[(l, node.label)] = node.label
				node_set[node.label] += c1_node_set

				for j in range(i+1, len(node.child_nodes())):
					c2 = node.child_nodes()[j]
					c2_node_set = node_set[c2.label]
					for l1 in c1_node_set:
						for l2 in c2_node_set:
							lca_dict[(l1, l2)] = node.label
	return lca_dict

def get_number_of_leaves(tree):
	leaf_num_dict = {}
	for node in tree.traverse_postorder():
		if node.is_leaf():
			leaf_num_dict[node.label] = 1
		else:
			num = 0
			for c in node.child_nodes():
				num += leaf_num_dict[c.label]
			leaf_num_dict[node.label] = num
	return leaf_num_dict

def get_leaf_set(tree):
	leaf_dict = {}
	for node in tree.traverse_postorder():
		if node.is_leaf():
			leaf_dict[node.label] = [node.label]
		else:
			leaves = []
			for c in node.child_nodes():
				leaves += leaf_dict[c.label]
			leaf_dict[node.label] = leaves
	return leaf_dict

def compute_A_transpose_A(tree):
	AtA = {}
	lca_dict = compute_all_pairwise_LCAs(tree)
	leaf_num_dict = get_number_of_leaves(tree)
	num_leaves = leaf_num_dict[tree.root.label]
	for n in tree.traverse_preorder():
		for m in tree.traverse_preorder():
			if not n.is_root() and not m.is_root():
				if (m.label, n.label) not in AtA:
					if m == n:
						AtA[(n.label, m.label)] = leaf_num_dict[m.label] * (num_leaves - leaf_num_dict[n.label])
					else:
						if (n.label, m.label) in lca_dict:
							lca = lca_dict[(n.label, m.label)]
						else:
							lca = lca_dict[(m.label, n.label)]
						if lca != n.label and lca != m.label:
							AtA[(n.label, m.label)] = leaf_num_dict[n.label] * leaf_num_dict[m.label]
						elif lca == n.label:
							AtA[(n.label, m.label)] = leaf_num_dict[m.label] * (num_leaves - leaf_num_dict[n.label])
						else:
							AtA[(n.label, m.label)] = leaf_num_dict[n.label] * (num_leaves - leaf_num_dict[m.label])
	return AtA

def compute_A_transpose_d(tree, refTree):
	leaf_set_dict = get_leaf_set(tree)
	all_leaves = leaf_set_dict[tree.root.label]
	Atd = {}
	for n in tree.traverse_preorder():
		# print("query tree edge: ", n.label)
		leaf_count = {}
		length_sum = {}
		total_sum = 0
		existing_pairs = 0
		bottom_leaves = leaf_set_dict[n.label]
		true_pairs = len(bottom_leaves) * (len(all_leaves) - len(bottom_leaves))
		for node in refTree.traverse_postorder():
			if node.is_leaf():
				if node.label in bottom_leaves:
					leaf_count[node.label] = [1, 0]
					length_sum[node.label] = [node.edge_length, 0]
				else:
					leaf_count[node.label] = [0, 1]
					length_sum[node.label] = [0, node.edge_length]
			else:
				count = [0, 0]
				bl_sum = [0, 0]
				for c1 in node.child_nodes():
					c1_count = leaf_count[c1.label]
					c1_sum = length_sum[c1.label]
					count[0] += c1_count[0]
					count[1] += c1_count[1]
					bl_sum[0] += c1_sum[0]
					bl_sum[1] += c1_sum[1]
					for c2 in node.child_nodes():
						if c1 != c2:
							c2_count = leaf_count[c2.label]
							total_sum += c1_sum[0] * leaf_count[c2.label][1] + c1_sum[1] * leaf_count[c2.label][0]
							existing_pairs += c1_count[0] * c2_count[1]
				if not node.is_root():
					bl_sum[0] += node.edge_length * count[0]
					bl_sum[1] += node.edge_length * count[1]
					leaf_count[node.label] = count
					length_sum[node.label] = bl_sum
		if not n.is_root():
			if existing_pairs == 0:
				Atd[n.label] = total_sum
			else:
				Atd[n.label] = total_sum * true_pairs / existing_pairs
	return Atd

def find_closest_leaf(tree, refTree):
	ref_leaf_set = [l.label for l in refTree.traverse_leaves()]
	dist_mat = tree.distance_matrix(leaf_labels=True)

	closest_leaf = {}
	for l in tree.traverse_leaves():
		if l.label not in ref_leaf_set:
			dist_mat[l.label] = {key: dist_mat[l.label][key] for key in dist_mat[l.label] if key in ref_leaf_set}
			closest_leaf[l.label] = min(dist_mat[l.label], key=dist_mat[l.label].get)
	return closest_leaf, ref_leaf_set

def compute_A_transpose_d_combined(tree, refTree):
	closest_leaf, ref_leaf_set = find_closest_leaf(tree, refTree)
	# ref_leaf_set = [l.label for l in refTree.traverse_leaves()]
	leaf_set_dict = get_leaf_set(tree)
	all_leaves = leaf_set_dict[tree.root.label]
	Atd = {}
	for n in tree.traverse_preorder():
		# print("query tree edge: ", n.label)
		leaf_count = {}
		length_sum = {}
		total_sum = 0
		# existing_pairs = 0
		bottom_leaves = leaf_set_dict[n.label]
		rest_leaves = [l for l in all_leaves if l not in bottom_leaves]
		true_pairs = len(bottom_leaves) * (len(all_leaves) - len(bottom_leaves))
		ref_pairs = len([l for l in bottom_leaves if l in ref_leaf_set]) * len([l for l in rest_leaves if l in ref_leaf_set])

		if ref_pairs == 0:
			missing_bottom = {}
			missing_rest = {}
			for l in bottom_leaves:
				if l in closest_leaf:
					cl = closest_leaf[l]
					if cl not in missing_bottom:
						missing_bottom[cl] = 1
					else:
						missing_bottom[cl] += 1

			for l in rest_leaves:
				if l in closest_leaf:
					cl = closest_leaf[l]
					if cl not in missing_rest:
						missing_rest[cl] = 1
					else:
						missing_rest[cl] += 1

			for node in refTree.traverse_postorder():
				if node.is_leaf():
					mb = 0
					mr = 0
					if node.label in missing_bottom:
						mb = missing_bottom[node.label]
					if node.label in missing_rest:
						mr = missing_rest[node.label]
					if node.label in bottom_leaves:
						leaf_count[node.label] = [1 + mb, mr]
						length_sum[node.label] = [node.edge_length * (1 + mb), node.edge_length * mr]
						# existing_pairs += (1 + mb) * mr
					elif node.label in all_leaves:
						leaf_count[node.label] = [mb, 1 + mr]
						length_sum[node.label] = [node.edge_length * mb, node.edge_length * (1 + mr)]
						# existing_pairs += (1 + mr) * mb
					else:
						leaf_count[node.label] = [0, 0]
						length_sum[node.label] = [0, 0]
				else:
					count = [0, 0]
					bl_sum = [0, 0]
					for c1 in node.child_nodes():
						c1_count = leaf_count[c1.label]
						c1_sum = length_sum[c1.label]
						count[0] += c1_count[0]
						count[1] += c1_count[1]
						bl_sum[0] += c1_sum[0]
						bl_sum[1] += c1_sum[1]
						for c2 in node.child_nodes():
							if c1 != c2:
								c2_count = leaf_count[c2.label]
								total_sum += c1_sum[0] * leaf_count[c2.label][1] + c1_sum[1] * leaf_count[c2.label][0]
								# existing_pairs += c1_count[0] * c2_count[1]
					if not node.is_root():
						bl_sum[0] += node.edge_length * count[0]
						bl_sum[1] += node.edge_length * count[1]
						leaf_count[node.label] = count
						length_sum[node.label] = bl_sum
			if not n.is_root():
				Atd[n.label] = total_sum

		else:
			for node in refTree.traverse_postorder():
				if node.is_leaf():
					if node.label in bottom_leaves:
						leaf_count[node.label] = [1, 0]
						length_sum[node.label] = [node.edge_length, 0]
					elif node.label in all_leaves:
						leaf_count[node.label] = [0, 1]
						length_sum[node.label] = [0, node.edge_length]
					else:
						leaf_count[node.label] = [0, 0]
						length_sum[node.label] = [0, 0]
				else:
					count = [0, 0]
					bl_sum = [0, 0]
					for c1 in node.child_nodes():
						c1_count = leaf_count[c1.label]
						c1_sum = length_sum[c1.label]
						count[0] += c1_count[0]
						count[1] += c1_count[1]
						bl_sum[0] += c1_sum[0]
						bl_sum[1] += c1_sum[1]
						for c2 in node.child_nodes():
							if c1 != c2:
								c2_count = leaf_count[c2.label]
								total_sum += c1_sum[0] * leaf_count[c2.label][1] + c1_sum[1] * leaf_count[c2.label][0]
								# existing_pairs += c1_count[0] * c2_count[1]
					if not node.is_root():
						bl_sum[0] += node.edge_length * count[0]
						bl_sum[1] += node.edge_length * count[1]
						leaf_count[node.label] = count
						length_sum[node.label] = bl_sum
			if not n.is_root():
				Atd[n.label] = total_sum * true_pairs / ref_pairs
	return Atd


def compute_A_transpose_d_naive(tree, refTree):
	dist_dict = refTree.distance_matrix(leaf_labels=True)
	leaf_set_dict = get_leaf_set(tree)
	all_leaves = leaf_set_dict[tree.root.label]
	Atd = {}
	for n in tree.traverse_preorder():
		bottom_leaves = leaf_set_dict[n.label]
		total_sum = 0
		for i in all_leaves:
			for j in all_leaves:
				if i != j:
					if i in bottom_leaves and j not in bottom_leaves:
						total_sum += dist_dict[i][j]
		Atd[n.label] = total_sum
	return Atd

def get_edge_indices(tree):
	index_to_edge = [n.label for n in tree.traverse_preorder()]
	root = tree.root
	index_to_edge.remove(root.label)
	if len(root.child_nodes()) == 2:
		index_to_edge.remove(root.child_nodes()[1].label)
	edge_to_index = {index_to_edge[i]:i for i in range(len(index_to_edge))}

	return edge_to_index, index_to_edge
	
def get_matrices(tree, refTree):
	edge_to_index, index_to_edge = get_edge_indices(tree)
	n_edges = len(edge_to_index)
	edge_lengths = {n.label:n.edge_length for n in tree.traverse_preorder()}
	root = tree.root
	if len(root.child_nodes()) == 2:
		child1 = root.child_nodes()[0]
		child2 = root.child_nodes()[1]
		sum_lengths = edge_lengths[child1.label] + edge_lengths[child2.label]
		edge_lengths[child1.label] = sum_lengths
		del edge_lengths[child2.label]
	# print(edge_lengths)


	Atd = compute_A_transpose_d_combined(tree, refTree)
	AtA = compute_A_transpose_A(tree)
	# with open(output_file + ".ata", "w") as f: 
	# 	for k in AtA:
	# 		f.write(str(k) + "\t" + str(AtA[k]) + "\n")
	Atd_mat = np.zeros(n_edges)
	edge_mat = np.zeros(n_edges)
	# AtA_mat = np.zeros((n_edges, n_edges))
	AtA_mat = np.zeros((n_edges, n_edges), dtype=np.uint32)
	for i in range(n_edges):
		Atd_mat[i] = Atd[index_to_edge[i]]
		edge_mat[i] = edge_lengths[index_to_edge[i]]
		for j in range(i, n_edges):
			i_label = index_to_edge[i]
			j_label = index_to_edge[j]
			if (i_label, j_label) in AtA:
				AtA_mat[i,j] = AtA[(i_label, j_label)]
				AtA_mat[j,i] = AtA[(i_label, j_label)]
			else:
				AtA_mat[i,j] = AtA[(j_label, i_label)]
				AtA_mat[j,i] = AtA[(j_label, i_label)]

	return Atd_mat, AtA_mat, edge_mat, edge_to_index, index_to_edge

def compute_optimal_rates(tree, refTree, r = 1):
	Atd, AtA, w, edge_to_index, index_to_edge = get_matrices(tree, refTree)
	# mean_diag = np.mean(np.diagonal(AtA))
	# sum_diag = np.mean(np.diagonal(AtA))
	# let n be the number of edges
	n = len(w) # 
	# print(n**2, sum_diag)
	# x is the optimal weight parameter (not the scale)
	x = cp.Variable(n)

	objective  = cp.quad_form(x, cp.psd_wrap(AtA)) - 2 * Atd.T @ x
	threshold = 0

	#run without regularization
	prob = cp.Problem(cp.Minimize(objective),[x >= threshold])
	prob.solve(verbose = False)
	w_0 = x.value
	# print(w_0)
	initial_obj = np.fabs((w_0.T @ AtA @ w_0) - 2 * Atd.T @ w_0)

	if r > 0:
		# print("r > 0")
		w_mask = w > 1e-8
		regularization = cp.sum_squares((x[w_mask] * 1e+4)/ (w[w_mask] * 1e+4) - (cp.sum((x[w_mask] * 1e+4)/(w[w_mask] * 1e+4)) / np.sum(w_mask))) / np.sum(w_mask)
		prob = cp.Problem(cp.Minimize(objective / initial_obj + r * regularization),[x >= threshold])
		prob.solve(verbose = False, warm_start=True)

	new_tree = read_tree_newick(tree.newick())
	# rates = ['0' for i in range(n)]
	# bls = ['0' for i in range(n)]

	# is_rooted = False
	# is_balanced = True
	if len(new_tree.root.child_nodes()) == 2:
		# is_rooted = True
		sum_root_bls = new_tree.root.child_nodes()[0].edge_length + new_tree.root.child_nodes()[1].edge_length

	rates = []
	bls = []

	# for n in edge_to_index:
	# 	print(n, x.value[edge_to_index[n]])
	# print(x.value)
	if x.value is None:
		print("unsuccessful")
		for n in new_tree.traverse_postorder():
			n.edge_length = 0
			rate = 0
			rates.append(str(rate))
			bls.append(str(n.edge_length))
		return new_tree.newick(), str(objective.value), rates, bls

	dropped_index=[]
	index=-1
	for n in new_tree.traverse_postorder():
		# print(n.label)
		# if n.label in edge_to_index:
		# 	print(x.value[edge_to_index[n.label]])
		# print("="*10)
		index+=1
		if n.is_root():
			dropped_index.append(index)
			continue
		dupl = False
		old_edge = n.edge_length
		#ignore I1?
		if n.parent.is_root() and len(new_tree.root.child_nodes()) == 2:
			new_length = x.value[edge_to_index[new_tree.root.child_nodes()[0].label]]
			n.edge_length = new_length * old_edge / sum_root_bls
			if n == new_tree.root.child_nodes()[1]:
				dupl=False #tentatively not ignoring
				# dupl = True
		else:
			new_length = x.value[edge_to_index[n.label]]
			n.edge_length = new_length

		if not dupl:
			rate = n.edge_length / old_edge
			rates.append(str(rate))
			bls.append(str(new_length))
		else:
			dropped_index.append(index)
			# print("dupl:",n)
	# print(new_tree.newick())

	return new_tree.newick(), str(objective.value), rates, bls,dropped_index


def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', required=True, help="Input tree")
	parser.add_argument('-r', '--ref', required=True, help="Reference trees")
	parser.add_argument('-l', '--reg_coef', default=0.001, required=False, help="Regularization Coefficient")
	parser.add_argument('-o', '--output_file', required=True, help="Output file")
	parser.add_argument('-g', '--log_file', required=False, help="Log file")


	args = parser.parse_args()

	random.seed(a=1105)
	np.random.seed(1105)

	start = time.time()

	with open(args.input,'r') as f:
		inputTree = f.read().strip().split("\n")[0]

	with open(args.ref,'r') as f:
		refTrees = f.read().strip().split("\n")

	new_trees = []
	losses = []
	all_rates = []
	all_bls = []

	for i in range(len(refTrees)):
		# print(i+1)
		# print(refTrees[i])
		# print("="*200)
		num = str(i+1)
		r = refTrees[i]
		inputTree_obj = read_tree_newick(inputTree)
		refTree_obj = read_tree_newick(r)
		__label_tree__(inputTree_obj)
		__label_tree__(refTree_obj)
		# print(refTree_obj.newick())
		new_tree, obj, rates, bls = compute_optimal_rates(inputTree_obj, refTree_obj, r = float(args.reg_coef))
		new_trees.append(new_tree)
		losses.append([num,obj])
		all_rates.append([num] + rates)
		all_bls.append([num] + bls)

	with open(args.output_file + ".trees", 'w') as f:
		for t in new_trees:
			f.write(t + '\n')

	if args.log_file:
		with open(args.log_file, 'w') as f:
			for l in losses:
				f.write('\t'.join(l) + '\n')

	with open(args.output_file + ".rates", 'w') as f:
		for r in all_rates:
			f.write('\t'.join(r) + '\n')

	with open(args.output_file + ".branches", 'w') as f:
		for r in all_bls:
			f.write('\t'.join(r) + '\n')


	end = time.time()
	print("Runtime: ", end - start) 


if __name__ == "__main__":
	main()