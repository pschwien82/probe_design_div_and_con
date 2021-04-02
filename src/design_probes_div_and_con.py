import os

# up the recursion limit or else will run into "max recursion depth exceeded" error on k>=4
#import sys
#sys.setrecursionlimit(10000) # still not enough on my system, using iterative versions of all algorithms
#resource.setrlimit(resource.RLIMIT_STACK, [0x10000000, resource.RLIM_INFINITY]) # only available on *nix machines

import time
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaIterator
from argparse import ArgumentParser
from anytree import AnyNode, Node, RenderTree, AsciiStyle
import matplotlib.pyplot as plt
from collections import deque # more efficient stack

# comment out any amino acids that should be ignored
std_alphabet = {}
std_alphabet["G"] = "Glycine (Gly)"
std_alphabet["A"] = "Alanine (Ala)"
std_alphabet["L"] = "Leucine (Leu)"
std_alphabet["M"] = "Methionine (Met)"
std_alphabet["F"] = "Phenylalanine (Phe)"
std_alphabet["W"] = "Tryptophan (Trp)"
std_alphabet["K"] = "Lysine (Lys)"
std_alphabet["Q"] = "Glutamine (Gln)"
std_alphabet["E"] = "Glutamic Acid (Glu)"
std_alphabet["S"] = "Serine (Ser)"
std_alphabet["P"] = "Proline (Pro)"
std_alphabet["V"] = "Valine (Val)"
std_alphabet["I"] = "Isoleucine (Ile)"
std_alphabet["C"] = "Cysteine (Cys)"
std_alphabet["Y"] = "Tyrosine (Tyr)"
std_alphabet["H"] = "Histidine (His)"
std_alphabet["R"] = "Arginine (Arg)"
std_alphabet["N"] = "Asparagine (Asn)"
std_alphabet["D"] = "Aspartic Acid (Asp)"
std_alphabet["T"] = "Threonine (Thr)"
# std_alphabet["U"] = "*** Selenocysteine (Sec)"
# std_alphabet["X"] = "*** Any / Unknown (Xaa)"
# std_alphabet["B"] = "*** Aspartic acid or Asparagine (Asx)"
# std_alphabet["Z"] = "*** Glutamic acid or Glutamine (Glx)"
# std_alphabet["J"] = "*** Leucine or Isoleucine (Xle)"


def read_fasta(filename):
	"""
	Reads all protein sequences into memory
	"""
	protein_seqs = []
	protein_ids = []
	with open(filename) as handle:
		for record in FastaIterator(handle):
			protein_ids.append(str(record.id))
			protein_seqs.append(str(record.seq))
			if len(protein_ids) % 5000 == 0:
				print("Read {} protein sequences.".format(len(protein_ids)), flush=True)
	return protein_seqs, protein_ids


def determine_alphabet(protein_seqs):
	"""
	Scans the proteome for the used protein alphabet in search of non-canonical amino acids
	"""
	alphabet = {}
	for protein in protein_seqs:
		for aa in protein:
			alphabet[aa] = alphabet.get(aa, 0) + 1
	return alphabet


def count_kmers(protein_seqs, k, protein_indices=[]):
	"""
	Builds kmer index of given protein indices
	"""
	kmers = {}
	protein_index = 0
	for protein in protein_seqs:
		if not protein_indices or protein_index in protein_indices:
			if len(protein) >= k:
				for i in range(len(protein) - k + 1):
					kmer = protein[i:i+k]
					if kmer in kmers:
						if protein_index in kmers[kmer]:
							kmers[kmer][protein_index] += 1
						else:
							kmers[kmer][protein_index] = 1
					else:
						kmers[kmer] = {protein_index: 1}
		protein_index += 1
	return kmers


def curate_kmers(kmers, non_canonical_amino_acids=[]):
	"""
	Removes kmers that do not conform to canonical amino acid code
	"""
	kmer_delete_list = []
	for kmer in kmers:
		for k in kmer:
			if k in non_canonical_amino_acids:
				kmer_delete_list.append(kmer)
				break
	if kmer_delete_list:
		print("Eliminated the following {} kmers due to containing non-canonical amino acids:".format(len(kmer_delete_list)))
		print(kmer_delete_list)
		for kmer in kmer_delete_list:
			del kmers[kmer]
	return kmers


def build_kmers_from_index(proteins_lookup, protein_indices):
	"""
	Builds kmer index specifically for the given set of protein indices
	"""
	kmers = {}
	for protein_index in protein_indices:
		if protein_index in proteins_lookup: # there are rare cases in which a protein has no kmer, e.g. sp|P0DPR3|TRDD1_HUMAN
			for kmer in proteins_lookup[protein_index]:
				if kmer in kmers:
					kmers[kmer][protein_index] = proteins_lookup[protein_index][kmer]
				else:
					kmers[kmer] = {protein_index: proteins_lookup[protein_index][kmer]}
	return kmers


def report_results(protein_indices, used_kmers, used_kmer_polarity, total_protein_count, total_found, error=""):
	"""
	If the verbose flag is set this method prints our progress statistics during main execution
	"""
	string = used_kmer_polarity[0]+used_kmers[0]
	for i in range(1, len(used_kmers)):
		string += ","+used_kmer_polarity[i]+used_kmers[i]
	print("{}\t{}\t{} ({:.2f}%)\t{}".format(",".join(str(pi) for pi in protein_indices), ",".join(ukp+uk for ukp, uk in zip(used_kmer_polarity, used_kmers)), total_found, total_found/total_protein_count*100, error), flush=True)


def divide_and_conquer_iterative(total_protein_count, k, protein_lookup, verbose):
	"""
	This method does the interesting work of finding the most appropriate probes of size k in a greedy fashion who's absence/presence pattern splits 
	the subset of proteins by as close to 50% as possible. This should in theory minimize the number of probes needed to uniquely identify all
	proteins. Since each introduction of a new probe splits the subset into two, the resulting data structure is a binary tree with each internal node
	representing a probe (AA k-mer) as well as a "polarity", indicating required absence (-) or presence (+) of this probe for the subtree. Each leaf
	of the tree is one or more protein sequences (multiple if no differentiating k-mers could be found). The path from a leaf to root specifies the 
	probe set that uniquely identifies the protein.

	"""
	root = Node(name="root", children=[])
	stack = deque([{"refnode": root, "used_kmers": [], "used_kmer_polarity": [], "protein_indices": list(range(total_protein_count))}])
	
	total_found = 0 # keep track of progress and used for stats
	
	while len(stack) > 0:
		# get a work item (subtree) from the stack
		subtree_work = stack.pop()

		node = subtree_work["refnode"]
		used_kmers = subtree_work["used_kmers"]
		used_kmer_polarity = subtree_work["used_kmer_polarity"]
		protein_indices = subtree_work["protein_indices"]

		if len(protein_indices) <= 1:
			if len(protein_indices) == 1:
				total_found += len(protein_indices)
				if verbose: report_results(protein_indices, used_kmers, used_kmer_polarity, total_protein_count, total_found)
				Node(name=protein_indices, parent=node)
			else:
				print("!!! Found protein indices to be empty at node {}, this should never happen.".format(node))
			continue

		# count all kmers in the set of provided protein sequences
		kmers = build_kmers_from_index(protein_lookup, protein_indices)
		
		# record how many proteins this kmer positively IDs
		protein_indices_covered = []

		# sort the kmers so that those covering closest to 50% of the protein sequences are on the top
		sorted_kmers = sorted(kmers, key = lambda key: abs(len(kmers[key])-len(protein_indices)/2))
		
		# go through kmers and pick the one closes to 50% that has not already been used
		for kmer in sorted_kmers:
			# return if we encounter a kmer that covers 100% of the proteins as all subsequent kmers
			# will have the same coverage since they were sorted that way
			if len(kmers[kmer]) == len(protein_indices):
				total_found += len(protein_indices)
				if verbose: report_results(protein_indices, used_kmers, used_kmer_polarity, total_protein_count, total_found, "No kmers found that distinguish the protein IDs")
				Node(name=protein_indices, parent=node)
				break

			# take the kmer with the closest to 50% protein coverage that has not already been used
			if kmer not in used_kmers:
				# remember all the proteins that have been covered with this kmer
				protein_indices_covered = list(kmers[kmer].keys())
								
				# put two new work items on the stack
				# search for the next kmers maximally dividing the uncovered proteins
				# using list concatenation here to prevent passing lists by reference
				stack.append({"refnode": Node(kmer=kmer, name="-"+kmer, parent=node, children=[], polarity="-"), "used_kmers": used_kmers + [kmer], "used_kmer_polarity": used_kmer_polarity + ["-"], "protein_indices": list(set(protein_indices) - set(protein_indices_covered))})
				# and search for the next kmers maximally dividing the covered proteins
				stack.append({"refnode": Node(kmer=kmer, name="+"+kmer, parent=node, children=[], polarity="+"), "used_kmers": used_kmers + [kmer], "used_kmer_polarity": used_kmer_polarity + ["+"], "protein_indices": protein_indices_covered})
				break
		else:
			# handle case in which no kmer is left to be used
			total_found += len(protein_indices)
			if verbose: report_results(protein_indices, used_kmers, used_kmer_polarity, total_protein_count, total_found, "No more unused kmers to distinguish protein IDs")
			Node(name=protein_indices, parent=node)
			continue
				
	return root	


def create_protein_index(all_kmers):
	"""
	Build a protein index-based lookup table for fast access to kmer data
	"""
	protein_lookup = {}
	for kmer in all_kmers:
		for protein_index in all_kmers[kmer]:
			if protein_index in protein_lookup:
				protein_lookup[protein_index][kmer] = all_kmers[kmer][protein_index]
			else:
				protein_lookup[protein_index] = {kmer: all_kmers[kmer][protein_index]}
	return protein_lookup


def validate_results(root, protein_seqs):
	"""
	Double check that the absence/presence pattern of the probes is indeed correct for each protein
	by bruteforce testing the sequences
	"""	
	stack = deque()
	node = root
	stack.append(node)
	while len(stack) > 0:
		if node.children:
			stack.append(node)
			node = node.children[0]
		else:
			probe_stack = []
			leaf_node = node
			while node.parent.parent:
				for sequence_id in leaf_node.name:
					if node.parent.kmer in protein_seqs[sequence_id]:
						if node.parent.polarity == "-":
							print("Found kmer '{}' in protein '{}' that should not be present!".format(node.parent.kmer, sequence_id))
							print("Node: {}".format(leaf_node))
					elif node.parent.polarity == "+":
						print("Did not find kmer '{}' in protein '{}' that should be present!".format(node.parent.kmer, sequence_id))
						print("Node: {}".format(leaf_node))
				node = node.parent
			node = stack.pop()
			node = stack.pop()
			node = node.children[1]


def save_probes_per_protein_iterative(root, protein_ids, output_filename):
	"""
	Stores all probes neccessary to uniquely ID a protein in a tab separated file
	one line per protein
	"""
	stack = deque()
	used_kmers = {}
	unique = 0
	ambiguous = 0
	kmers_per_protein = []
	with open(output_filename, "w") as f:
		node = root
		stack.append(node)
		while len(stack) > 0:
			if node.children:
				stack.append(node)
				node = node.children[0]
			else:
				probe_stack = []
				internal_label = ';'.join([str(n) for n in node.name])
				label = ';'.join([protein_ids[n] for n in node.name])
				if len(node.name) == 1:
					unique += 1
				else:
					ambiguous += len(node.name)
				
				while node.parent.parent:
					probe_stack.append(node.parent.polarity+node.parent.kmer)
					# count occurence
					if node.parent.kmer in used_kmers:
						used_kmers[node.parent.kmer] += 1
					else:
						used_kmers[node.parent.kmer] = 1
					node = node.parent
				probe_stack.reverse()
				kmers_per_protein.append(len(probe_stack))
				f.write("{}\t{}\t{}\n".format(internal_label, label, ",".join(probe_stack)))
				node = stack.pop()
				node = stack.pop()
				node = node.children[1]
	return used_kmers, unique, ambiguous, kmers_per_protein


def get_newick_string_recursive(node):
	"""
	Simple recursive function to convert AnyNode-based tree into newick tree format
	Note: runs into recursion depth problems for k>3, so using iterative method below
	"""	
	label = children = ""
	
	# differentiate internal nodes and leaf node labels
	# the former is a concatenation of polarity and kmer
	# whereas the latter is a concatenation of protein IDs
	if node.name == "root":
		label = node.name
	elif node.children:
		label = node.polarity + node.kmer
	else:
		label = '-'.join([str(n) for n in node.name])

	# join subtree together by recursively traversing the children
	if node.children:
		children = ','.join([get_newick_string(n) for n in node.children])
	
	# wrap in parenthesis if children are not None
	if children:
		children = '(' + children + ')'

	return children + label


def get_newick_string_iterative(node):
	"""
	Iterative function to convert AnyNode-based tree into newick tree format.
	"""	
	stack = []
	res = label = ""
	while True:
		while node:
			if node.children:
				stack.append(node.children[0])
				res += '('

			stack.append(node)
			if len(node.children)==2:
				node = node.children[1]
			else:
				node = None
		node = stack.pop()
		if node.children and (stack and node.children[0] == stack[-1]):
			tmp = stack.pop()
			stack.append(node)
			if res[-1] == ')':
				res = res[:-1]
			if len(node.children)==2:
				res += ','
			node = tmp
		else:
			if node.name == "root":
				label = node.name
			elif node.children:
				label = node.polarity + node.kmer
			else:
				label = '-'.join([str(n) for n in node.name])
			res += label + ')'
			node = None
		if not stack:
			break
	res = res[:-1]
	return res + ";"


def main():
	argparser = ArgumentParser()
	argparser.add_argument('--k', type=int, help='set the probe\'s target size in amino acids', default='3')
	argparser.add_argument('proteome', type=str, help='path to fasta file with all protein sequences to design probes for (e.g. uniprot-proteome_UP000005640_canonical_plus_isoforms.fasta)')
	argparser.add_argument('--verbose', help='print progress information during divide & conquer', action="store_true")
	args = argparser.parse_args()

	# set kmer size
	k = args.k
	# show divide & conquer progress
	verbose = args.verbose
	# fasta proteome input file
	input_filename = args.proteome
	#input_filename = "C:\\Users\\psusa\\Downloads\\uniprot-proteome_UP000005640_canonical_plus_isoforms.fasta"
	#input_filename = "C:\\Users\\psusa\\Downloads\\uniprot-proteome_UP000005640_canonical.fasta"
	#input_filename = "C:\\Users\\psusa\\Downloads\\uniprot-proteome_UP000005640_canonical_mini.fasta"
	filename, file_extension = os.path.splitext(input_filename)
	
	# add date/time to output files to avoid overwrite
	timestr = time.strftime("%Y%m%d-%H%M%S")

	# setup output files
	output_filename = filename + "_" + timestr + "_{}-mer_run_results.txt".format(k)
	newick_file = filename + "_" + timestr + "_{}-mer.newick".format(k)
	tree_file = filename + "_" + timestr + "_{}-mer_tree.txt".format(k)

	# start work
	print("Reading protein file...", flush=True)
	protein_seqs, protein_ids = read_fasta(input_filename)
	print("Read {} proteins.\n".format(len(protein_seqs)), flush=True)

	print("Checking AA alphabet (non-canonical AAs will be ignored during probe design)...", flush=True)
	alphabet = determine_alphabet(protein_seqs)
	total_amino_acids = sum(alphabet.values())
	non_canonical_amino_acids = []
	for character in alphabet:
		if character in std_alphabet:
			print("{} {} - {:.4f}% ({})".format(character, std_alphabet[character], alphabet[character]/total_amino_acids*100, alphabet[character]))
		else:
			non_canonical_amino_acids.append(character)
			print("Non-canonical: {} - {:.4f}% ({})".format(character, alphabet[character]/total_amino_acids*100, alphabet[character]))
	print()

	start_time = time.monotonic()
	print(time.ctime(), " Counting kmers...", flush=True)
	kmers = count_kmers(protein_seqs, k)
	print("done in {:.1f} seconds, {} unique {}-mers found ({:.2f}% of {}^{} possible)\n".format(time.monotonic() - start_time, len(kmers), k, (len(kmers)/pow(len(alphabet), k)*100), len(alphabet), k), flush=True)

	start_time = time.monotonic()
	print(time.ctime(), " Removing non-canonical kmers...", flush=True)
	kmers = curate_kmers(kmers, non_canonical_amino_acids)
	print("done in {:.1f} seconds\n".format(time.monotonic() - start_time), flush=True)

	start_time = time.monotonic()
	print(time.ctime(), " Building protein lookup table...", flush=True)
	protein_lookup = create_protein_index(kmers)
	print("done in {:.1f} seconds\n".format(time.monotonic() - start_time), flush=True)

	del kmers # free up mem
	
	start_time = time.monotonic()
	print(time.ctime(), " Calculating probes...",  flush=True)
	root = divide_and_conquer_iterative(len(protein_seqs), k, protein_lookup, verbose)
	print("done in {:.1f} minutes\n".format((time.monotonic() - start_time)/60), flush=True)
	
	start_time = time.monotonic()
	print(time.ctime(), " Validating results...",  flush=True)
	validate_results(root, protein_seqs)
	print("done in {:.1f} seconds\n".format(time.monotonic() - start_time), flush=True)

	print("Writing results...\n", flush=True)
	used_kmers, unique, ambiguous, kmers_per_protein = save_probes_per_protein_iterative(root, protein_ids, output_filename)

	total_kmers = sum(used_kmers[kmer] for kmer in used_kmers)
	print("A total of {} unique ({} total) kmers have been used to uniquely identify {} of {} proteins ({} could not be uniquely resolved).".format(len(used_kmers), total_kmers, unique, len(protein_ids), ambiguous))
	print("Average number of probes to uniquely identify a protein is {:.2f}.\n".format(total_kmers/unique))

	try:
		print("Writing ASCII tree...\n", flush=True)
		with open(tree_file, "w") as f:
			f.write(RenderTree(root, style=AsciiStyle()).by_attr('name'))
	except RecursionError as re:
		print("Unable to render ASCII tree due to: {}".format(re.args[0]))
		print("This is expected with k>=4.\n")

	print("Writing Newick tree...\n", flush=True)
	with open(newick_file, "w") as f:
		#f.write(get_newick_string_recursive(root)) # my system runs out of steam for k>3, using iterative instead but less elegant
		f.write(get_newick_string_iterative(root))

	# matplotlib histogram of how many kmers are needed for a unique protein id
	plt.hist(kmers_per_protein, color = 'blue', align='left', edgecolor = 'blue', bins=max(max(kmers_per_protein)-min(kmers_per_protein)+1, 10))
	plt.title('Histogram of {}-mer probes needed to uniquely id proteins'.format(k))
	plt.xlabel('number of {}-mer probes'.format(k))
	plt.ylabel('number of proteins')
	plt.show()


if __name__ == '__main__':
	main()
