#!/usr/bin/env python
__author__ = 'Benjamin Bolduc'

import sys
import os
import optparse
import math
import csv
import datetime
from collections import OrderedDict
from pprint import pprint

p = optparse.OptionParser(description="""A simple script that creates a graph (*.graphml format) from a BLAST file and
applies a community detection algorithm to identify groups of highly-related sequences. Annotation information can also
be added as well, allowing high-order organization and graphics in a graphml-compatible graph visualization software.
A number of options exist to filter the BLAST HSPs going into the network, thresholds for the number of minimum members
within a group to be retained in the graph, as well as read-based analytics (using *.ace formatted files). Limited data
use visualization is also possible.""")

g = optparse.OptionGroup(p, "General Options")

g.add_option("--input_fn", dest="input_fn", metavar="FILENAME",
			help="""BLAST file in pairwise format (out_fmt 0 option).""")

g.add_option("--annotation_fn", dest="annotation_fn", metavar="FILENAME",
			help="""Tab-delimited file with header columns as arguments for annotation and the first row as the
			sequence/contig/node's name. This name MUST MATCH the sequence name in whatever fasta file was used to
			create the network from. EX:

			Name    Length  Group   Time    Spacer  Virus Type  Fraction
			Contig5 2554    8       Feb     Yes     Archaeal    17
			Contig2 3221    6       Dec     No      Bacterial   19""")

g.add_option("--output_fn", dest="output_fn", metavar="FILENAME", default="Output", help="""Output prefix.""")

p.add_option_group(g)

g = optparse.OptionGroup(p, "BLAST Options")

g.add_option("--evalue", dest="evalue", type="float", metavar="FLOAT", default=1E-30,
			help="""E-value cutoff for sequence selection (optional; default: %default).""")

g.add_option("--hsp_to_query_ratio", dest="hsp_to_query_ratio", type="float", metavar="FLOAT", default=0.0,
			help="""Percent of query that belongs to the HSP. If alignment/HSP is 200 bp and query length is 400,
			the ratio is 0.5 (optional; default: %default).""")

g.add_option("--hsp_identity", dest="hsp_identity", type="float", metavar="FLOAT", default=0.75,
			help="""Percent identity cutoff. This allows a minimum % identity that an alignment must be to be included.
			(optional; default: %default).""")

g.add_option("--hsp_length", dest="hsp_length", type="int", metavar="INTEGER", default=200,
			help="""HSP length cutoff. This allows a minimum length that an alignment must be to be included.
			(optional; default: %default).""")

p.add_option_group(g)

g = optparse.OptionGroup(p, "Network Options")

g.add_option("--num_inputs", dest="num_inputs", type="int", metavar="INTEGER",
			help="""The number of sequences in the original BLAST. This is useful for when the BLAST file is too large
			to be easily parsed.""")

g.add_option("--partitioner", dest = "partitioner", metavar = "PYTHON CODE", default = "%",
			help = """Python code to further define a separating characteristic (i.e. site, date, etc) that can be used
			during graph visualization to color-code nodes. For example, "%.split('|')[1].split('_')[1]" will identify
			0808 for contig01129|NL10_0808_vDNA_MSU_WGA; '%' will be replaced by the sequence name. An alternative is to
			use the annotation file, whose attributes will be passed to the created nodes (optional; Default: %default")""")

g.add_option("--use_reads", dest="use_reads", action="store_true", default=False,
			help="""Enables the use of read counting using 454 ACE files.""")

g.add_option("--read_fn", dest="read_fn", metavar="FILENAME",
			help="""Ace file with read-to-contig mapping information.""")

g.add_option("--read_community_cutoff", dest="read_co_cutoff", type="int", metavar="INTEGER", default=1000,
			help="""Sets a cutoff for the minimum number of reads comprising the contigs of a community must have to be retained.""")

g.add_option("--contig_community_cutoff", dest="contig_co_cutoff", type="int", metavar="INTEGER", default=50,
			help="""Sets a cutoff for the minimum number of members a community must have to be retained.""")

g.add_option("--remove_tmp", dest="tmp_remove", action="store_true", default=False,
			help="""Removes initial network created prior to any manipulation. This is useful for when memory is
			constrained and analysis crashes. This REMOVES that file at the end of a successful run.""")

p.add_option_group(g)

g = optparse.OptionGroup(p, "Post-Visualization Options")

p.add_option("--membership_per_community", dest="mpc_fig", action="store_true", default=False,
			help="Will count the members in each community and generate a figure.")

p.add_option("--data_represented", dest="dr_fig", action="store_true", default=False,
			help="Will calculate the data represented and generate a figure.")

(p, a) = p.parse_args()

def error(msg):
	print >> sys.stderr, "ERROR: {}".format(msg)
	sys.exit(1)

if (p.input_fn == None):
	error("An input sequence filename is required.")

try:
	with open(p.input_fn): pass

except IOError:
	error("Input file {} not found.".format(p.input_fn))

try:
	from blastorage import Storage

except ImportError:
	error("The BlaSTorage library is required and was not found. Available @ http://biowiki.crs4.it/biowiki/MassimilianoOrsini")

if p.mpc_fig or p.dr_fig:

	try:
		import matplotlib
	except ImportError:
		error("The matplotlib library is required and was not found. Available @ http://matplotlib.org/")

	try:
		import numpy as np
	except ImportError:
		error("The numpy library is required and was not found. Available @ ")

try:
	import networkx as nx

except ImportError:
	error("The NetworkX library is required and was not found. Available @ https://networkx.github.io/")

try:
	import community

except ImportError:
	error("The NetworkX-compatible Louvain community detection algorithm is required and was not found."
	      "It is available at https://bitbucket.org/taynaud/python-louvain")

if (p.use_reads):

	if (p.read_fn == None):
		error("An read filename is required.")

	try:
		with open(p.read_fn):
			pass

	except IOError:
		error("Read file {} not found.".format(p.read_fn))


def logger(logfile_handle, text):
	logfile_handle.write("{} {}\n".format(datetime.datetime.now(), text))
	logfile_handle.flush()
	print text

logfile_fn = "{}.log".format(p.output_fn)
logfile_fh = open(logfile_fn, 'a', 0)  # Set unbuffered, so info written to disk immediately

diG = nx.DiGraph()

node_annotation_dict = {}

if p.annotation_fn:
	with open(p.annotation_fn, 'rU') as annotation_fh:

		logger(logfile_fh, "Reading Annotation File.")

		line_count = 0
		header = []

		for lines_of_data in annotation_fh:

			if line_count == 0:  # Hack to separate 1st line of file (annotation names) from rest of file
				header = lines_of_data.strip().split('\t')
				line_count += 1

			else:
				data = lines_of_data.strip().split('\t')
				node_annotation_dict[data[0]] = dict(zip(header[1:], data[1:]))  # With header columns as keys, rows as values
				# Header should stay the same, with lengths, groups, etc always being keys and the data (changing) as values

db = p.input_fn.rsplit('.', 1)[0]+'.db'
infile = p.input_fn

graph_tmp_out = "{0}.tmp.graphml".format(p.output_fn)
uGraph = ''
input_queries = 0
non_singleton_graph_node_count = 0

NoMatches = 0
NoMatches_fh = open(p.output_fn + '-NoMatches.list', 'a', 0)

if os.path.isfile("{}.graphml".format(p.output_fn)):
	print "End file {}.graphml already exists. It is safe to delete this file.".format(p.output_fn)

else:
	if os.path.isfile(graph_tmp_out):

		logger(logfile_fh, "Graph file found! Jumpstarting Analysis.")
		uGraph = nx.read_graphml(graph_tmp_out)
		logger(logfile_fh, "Graph file loaded.")

		non_singleton_graph_node_count = uGraph.number_of_nodes()  # Update nodes if reloading graph file

	if not os.path.isfile(graph_tmp_out):

		logger(logfile_fh, "Processing main BLAST file.")

		st = Storage(db, infile)
		st.store()
		api = st()

		logger(logfile_fh, "BLAST file conversion finished. Parsing for network.")

		blast_queries = api.getQueriesNumber()
		input_queries = blast_queries  # Number of reference sequences should/cold be subtracted

		logger(logfile_fh, "There are {} blast queries to parse.".format(blast_queries))

		n = 0

		for blast in api:  # One Blast result for each query

			query = str(blast.getQueryTitle()).strip()
			query_len = int(blast.getQueryLength())

			query_dict = {
				"id": query,
				"length": query_len,
				"size": 1,
			}

			if query in node_annotation_dict:
				query_dict.update(node_annotation_dict[query])

			# Round result, multiple objects for psiblast... apparently only 1 (stills needs iterated through though)
			for iteration in blast:

				if not iteration.found:  # Another way of saying "if cond is False"

					NoMatches += 1
					NoMatches_fh.write(query + '\n')
					NoMatches_fh.flush()

				else:

					# All sequences will get 'a chance.' If good match, then edge will connect.
					# If no good match (i.e. fails >= 1 conditions), will be removed during singleton analysis
					diG.add_node(query, query_dict)  # Pass dict to networkx

					for alignment in iteration:  # Alignment, one for each subject

						target = str(alignment.getSubjectTitle()).strip()
						target_len = int(alignment.getSubjectLength())  # True target length, i.e. the contig size

						target_dict = {
							"id": target,
							"length": query_len,
							"size": 1,
						}

						if target in node_annotation_dict:
							target_dict.update(node_annotation_dict[target])

						for hsp in alignment:

							hsp.setGradesDict()

							evalue = hsp.getEvalueAsFloat()
							identity = hsp.getIdentitiesAsPercentage()
							identities = hsp.getIdentitiesAsAbsoluteValue()  # Int
							hsp_length = int(hsp.getIdentitiesExpression().split(' ')[0].split('/')[1])  # Hack to get hsp length, Int
							hsp_identity = float(identities) / hsp_length  # Percent identity over HSP length
							total_alignment_length = hsp_length / target_len  # Alignment length over target
							hsp_to_query_ratio = hsp_length / (query_len + .0)  # Can be greater than 1

							if query == target:
								continue  # At the very least, it's a self-match, at worst, there's already a node in the network

							if (evalue <= p.evalue) and (hsp_to_query_ratio >= p.hsp_to_query_ratio) \
								and (hsp_identity >= p.hsp_identity) and (hsp_length >= p.hsp_length):

								# Going to add a node anyway, might as well add it once
								diG.add_node(target, target_dict)

								if evalue == 0:
									diG.add_edge(query, target, weight=500.0, length=hsp_length, identity=hsp_identity, label='BLAST')

								if evalue != 0:
									expect = math.log(evalue) * -1
									diG.add_edge(query, target, weight=expect, length=hsp_length, identity=hsp_identity, label='BLAST')
			n += 1

			sys.stdout.write("\r{0:.4f}".format((float(n) / blast_queries) * 100))
			sys.stdout.flush()

		logger(logfile_fh, "\nThere were {} queries without a match.".format(NoMatches))
		logger(logfile_fh, "Finished processing BLAST file.")

		# Remove nodes with fewer than 1 edge in a directed graph
		# Following speeds things up by removing nodes that are not relevant
		# Any nodes with matches but not matches of sufficient strength
		singles = [node for node, degree in diG.degree_iter() if degree == 0]
		logger(logfile_fh, "There were {} queries that were singletons.".format(len(singles)))
		diG.remove_nodes_from(singles)

		Insufficient_Match_fh = open(p.output_fn + '-PoorMatches.list', 'w')

		formatted_singles = "\n".join(singles)
		print >> Insufficient_Match_fh, formatted_singles
		Insufficient_Match_fh.close()

		logger(logfile_fh, "Removed singletons.")
		logger(logfile_fh, "Converting directed graph to undirected.")

		# Convert graph to undirected - can't run community analysis without
		uGraph = diG.to_undirected(diG)
		non_singleton_graph_node_count = uGraph.number_of_nodes()

		logger(logfile_fh, "Writing raw graphml file to disk.")

		# In case downstream analysis fails... wont have to repeat
		nx.write_graphml(uGraph, "{0}.tmp.graphml".format(p.output_fn))

		logger(logfile_fh, "Attempting to close and remove BlaSTorage DB. If it fails, kill job and restart. It will resume."
							" This is likely to happen if your BLAST file is >= 2 GB.")

		st.close()
		os.remove(db)

		logger(logfile_fh, "Closing finished.")

	# Sadly needed if graph BlaSTorage fails to close, halting the program. No other way to obtain the input queries unless
	# the blast file is reloaded... so could either force user to supply the number, the fasta file used for the blast, or
	# just grep -c '>' on the original fasta file and be done with it
	if (p.num_inputs):
		input_queries = int(p.num_inputs)
	else:
		if input_queries == 0:
			error("There are no input queries or the number of input queries has not been specified. This can happen if"
			      " the BlaSTorage library failed previously.")

	def find_partition(graph):

		g = graph
		partition = community.best_partition(g)
		nx.set_node_attributes(g, 'orig_partition', partition)
		return g, partition

	logger(logfile_fh, "Identifying partitions...")

	# Find communities
	coGraph, parts = find_partition(uGraph)

	new_parts = {}
	cluster_membership = {}

	out_part_mem_fh = open(p.output_fn + '-Partition-Membership.tab', 'w')

	for sequences, partitions in parts.iteritems():

		print >> out_part_mem_fh, '{0}\t{1}'.format(sequences, partitions)

		if not partitions in new_parts.keys():  # The keys are partition numbers
			new_parts[partitions] = 1

		else:
			new_parts[partitions] += 1

		if not partitions in cluster_membership:
			cluster_membership[partitions] = [sequences]
		else:
			cluster_membership[partitions].append(sequences)

	out_part_mem_fh.close()

	if p.partitioner != "%":

		logger(logfile_fh, "Partitioner enabled due to non-default options.")

		logger(logfile_fh, "Identifying fraction information...")

		part_site_comp_fh = csv.writer(open(p.output_fn + '-Partition-Site-Composition.csv', 'w'))

		partition_sites = {}

		for contigs, partitions in parts.iteritems():

			partition_num = str(partitions)

			get_characteristic = eval("lambda x: " + p.partitioner.replace('%', 'x').replace("\\x", '%'))

			characteristic = get_characteristic(contigs)

			if not partition_num in partition_sites:
				partition_sites[partition_num] = {}

			if not characteristic in partition_sites[partition_num]:
				partition_sites[partition_num][characteristic] = 0

			partition_sites[partition_num][characteristic] += 1

		characteristics = {}

		for cluster_id in partition_sites:

			for characteristic in partition_sites[cluster_id]:
				characteristics[characteristic] = None

		characteristics = sorted(characteristics)

		part_site_comp_fh.writerow(["Cluster"] + characteristics + ["Total Contigs"])

		for cluster_id in sorted(partition_sites):

			row = [cluster_id]

			total_contigs = sum(partition_sites[cluster_id].values())

			for characteristic in characteristics:
				row.append(partition_sites[cluster_id].get(characteristic, 0))

				# If reference sequence present, remove count from total
				if 'REF' in characteristic:
					total_contigs -= partition_sites[cluster_id].get(characteristic, 0)

			row.append(total_contigs)

			part_site_comp_fh.writerow(row)

	cutoff = []

	if (p.use_reads):

		print "Working on ace file {}".format(p.read_fn)

		contig_read_dict = {}
		contig_read_len_dict = {}

		from Bio.Sequencing import Ace

		with open(p.use_reads, 'rU') as ace_fh:

			for contig in Ace.parse(ace_fh):
				"""rd (reads) - read with name, sequence, etc
				qa (read qual) - which parts used as consensus
				ds - file name of read's chromatogram file
				af - loc of read within contig
				bs (base segment) - which read chosen at consensus at each pos
				rt (transient read tags) - generated by crossmatch and phrap
				ct (consensus tag)
				wa (whole assembly tag) - hosts assembly program name, version, etc
				wr
				reads - info about read supporting ace contig
				contig - holds info about contig from ace record"""

				contig_name = "{}".format(contig.name)  # contig00001

				if not contig_name in contig_read_dict:
					contig_read_dict[contig_name] = []
				if not contig_name in contig_read_len_dict:
					contig_read_len_dict[contig_name] = []

				for read_id, read in enumerate(contig.reads):
					# Using enumerate to give int because parsing through list of contig's reads. The list isn't necessarily in order
					# meaning that parsing in one block won't correspond to to its location (in the list) in another block
					read_name = read.rd.name.split('.')[0]  # Functionally equivalent to contig.af[read_id].name
					read_length = len(read.rd.sequence)

					if not read_name in contig_read_dict[contig_name]:
						contig_read_dict[contig_name].append(read_name)
						contig_read_len_dict[contig_name].append(read_length)

		cluster_membership = OrderedDict(sorted(cluster_membership.items(), key=lambda t: len(t[1]), reverse=True))  # Sorted by value

		cluster_read_abundances = OrderedDict()
		cluster_read_length_abundances = OrderedDict()

		for cluster, contigs in cluster_membership.iteritems():  # Ordered above, but need another dict to keep it that way

			# Might as well filter the list here and now. Later the reads *would of* been filtered, but that's pointless since
			# N reads will never be sufficient to generate N contigs in order to be included in network

			for contig in contigs:

				try:
					num_reads = len(contig_read_dict[contig])  # Get length of contig's read list
					read_lengths = sum(contig_read_len_dict[contig])  # Get sum of all read lengths

				except KeyError:  # References will not be in the metagenome's ace files, obviously
					print "Skipping {}".format(contig)
					break

				if not cluster in cluster_read_abundances:
					cluster_read_abundances[cluster] = 0
				if not cluster in cluster_read_length_abundances:
					cluster_read_length_abundances[cluster] = 0

				cluster_read_abundances[cluster] += num_reads
				cluster_read_length_abundances[cluster] += read_lengths

		out_part_size_fh = open(p.output_fn + '-Partition-Contig-Read-Sizes.tab', 'w')

		for keys, values in new_parts.iteritems():  # partition: count

			try:
				cluster_reads = cluster_read_abundances[keys]

				print >> out_part_size_fh, '{0}\t{1}\t{2}'.format(keys, values, cluster_reads)  # print partition \t count

				if values >= p.read_co_cutoff:
					cutoff.append(keys)

			except KeyError:  # All clusters are in new_parts, not all clusters have sufficient
				print "Unable to find {}. It has {} contigs".format(keys, values)

		out_part_size_fh.close()

		if p.partitioner != "%":

			part_site_read_abund_fh = csv.writer(open(p.output_fn + '-Partition-Site-Read-Abundances.csv', 'w'))

			partition_read_abundance_sites = {}

			for contigs, partitions in parts.iteritems():

				partition_num = str(partitions)

				try:
					num_reads = len(contig_read_dict[contigs])
					#print contigs, partitions, num_reads
				except KeyError:  # References will not be in the metagenome's ace files, obviously
					print "Skipping {}".format(contigs)
					num_reads = 0

				get_characteristic = eval("lambda x: " + p.partitioner.replace('%', 'x').replace("\\x", '%'))

				characteristic = get_characteristic(contigs)

				if not partition_num in partition_read_abundance_sites:
					partition_read_abundance_sites[partition_num] = {}

				if not characteristic in partition_read_abundance_sites[partition_num]:
					partition_read_abundance_sites[partition_num][characteristic] = 0

				partition_read_abundance_sites[partition_num][characteristic] += num_reads

			characteristics = {}

			for cluster_id in partition_read_abundance_sites:

				for characteristic in partition_read_abundance_sites[cluster_id]:
					characteristics[characteristic] = None

			characteristics = sorted(characteristics)

			part_site_read_abund_fh.writerow(["Cluster"] + characteristics + ["Total Reads", "Total Contigs"])

			read_composition_array = []

			for cluster_id in sorted(partition_read_abundance_sites):

				row = [cluster_id]

				total_reads = sum(partition_read_abundance_sites[cluster_id].values())
				total_contigs = new_parts[int(cluster_id)]
				list_to_append = []

				for characteristic in characteristics:
					row.append(partition_read_abundance_sites[cluster_id].get(characteristic, 0))

					if total_contigs >= p.co_cutoff:
						list_to_append.append(partition_read_abundance_sites[cluster_id].get(characteristic, 0))

				row.append(total_reads)
				row.append(total_contigs)

				if list_to_append:  # Dont want to add empty arrays
					read_composition_array.append(list_to_append)

				part_site_read_abund_fh.writerow(row)

			# http://stackoverflow.com/questions/19872530/python-sort-lists-based-on-their-sum
			read_composition_array.sort(key=sum)  # So easy it makes me want to cry

	# If not using reads
	else:

		out_part_size_fh = open(p.output_fn + '-Partition-Contig-Sizes.tab', 'w')

		for keys, values in new_parts.iteritems():  # partition: count

			print >> out_part_size_fh, '{0}\t{1}'.format(keys, values)  # print partition \t count

			if values >= p.contig_co_cutoff:
				cutoff.append(keys)

	logger(logfile_fh, "Removing clusters containing less than {} members".format(p.contig_co_cutoff))

	# Create a subgraph that includes only those above a specified cutoff value
	subGraph = coGraph.subgraph([n for n, attrdict in coGraph.node.items() if attrdict['orig_partition'] in cutoff])

	logger(logfile_fh, "Writing final network to disk...")

	nx.write_graphml(subGraph, "{0}.graphml".format(p.output_fn))

	logger(logfile_fh, "Network written to disk")

	if p.tmp_remove:

		os.remove(graph_tmp_out)
		logger(logfile_fh, "Removing temporary network file...")

	if p.mpc_fig or p.dr_fig:

		from mpl_toolkits.mplot3d import Axes3D
		import matplotlib.pyplot as plt
		import numpy as np
		from matplotlib import rcParams

		if p.mpc_fig:

			logger(logfile_fh, "Generating Membership Per Community Figure")

			rcParams['ytick.direction'] = 'out'

			font = {'family': 'sans-serif',
					'weight': 'normal',
					'size': 10}
			lines = {'linewidth': 0.5}

			Data = [value for (key, value) in new_parts.iteritems() if key in cutoff]
			Data.sort(reverse=True)

			# Removes the "long tail" often associated with rare species
			Data_NoFew = filter(lambda few: few >= 2, Data)

			fig = plt.figure(figsize=(8, 5))
			ax = fig.add_subplot(1, 1, 1)

			N = len(Data_NoFew)
			ind = np.arange(N)
			width = 0.8
			margin = 0.1

			colors = []
			for data in Data_NoFew:
				if data >= p.contig_co_cutoff:
					colors.append('b')
				else:
					colors.append('r')

			ax.bar(ind + margin, Data_NoFew, width, color=colors, align='edge', alpha=.8, zorder=1, edgecolor="none", linewidth=0)  # Removed label=

			ax.set_ylabel('Members/Community')
			ax.set_xlabel('Viral Group Rank')
			ax.legend(loc='upper left', frameon=False)

			fig.savefig(p.output_fn + "-Members-in-Each-Community.png", dpi=600, bbox_inches='tight')

		if p.dr_fig:

			logger(logfile_fh, "Generating NonSingleton Data Represented Figure.")

			import matplotlib.cm as cm
			from matplotlib.colors import Normalize

			font = {'family': 'sans-serif',
					'weight': 'normal',
					'size': 10}
			lines = {'linewidth': 0.5}

			matplotlib.rc('font', **font)
			matplotlib.rc('lines', **lines)

			fig = plt.figure(figsize=(7, 4))
			ax = fig.add_subplot(1, 1, 1)

			# new_parts = partition#, counts
			Data = [value for (key, value) in new_parts.iteritems()]

			N = 15  # Change to 16 if including 16 cutoff range
			ind = np.arange(N)
			width = 0.8

			Total_Contigs = int(input_queries)
			# Will create a list with values of the counts of each of the partitions (done above)
			Total_Contigs_In_Network_list = [value for value in Data if value >= 2]
			Total_Contigs_In_Network = sum(Total_Contigs_In_Network_list)
			Total_Contigs_In_Gt1 = non_singleton_graph_node_count

			data_represented_membership_size_cutoffs = [2, 3, 4, 5, 10, 15, 20, 30, 40, 50, 100, 200, 300, 400, 500]

			Data_Represented = []

			for cutoff in data_represented_membership_size_cutoffs:

				Data_prerepresented = []

				for community_sizes in Data:

					if community_sizes >= cutoff:
						Data_prerepresented.append(community_sizes)

				# Calculate % of data represented by the new communities that are above the cutoff
				to_append = (sum(Data_prerepresented) / (Total_Contigs_In_Gt1 + .0)) * 100

				Data_Represented.append(to_append)

			my_cmap = cm.get_cmap('jet')
			my_norm = Normalize(vmin=min(Data_Represented), vmax=max(Data_Represented))

			ax.bar(ind, Data_Represented, width, color=my_cmap(my_norm(Data_Represented)), align='edge', alpha=1, zorder=1)

			xlabels = ["{0}".format(x) for x in data_represented_membership_size_cutoffs]

			ax.set_xlabel('Members/Group')
			ax.set_ylabel('Contigs Represented (%)')
			ax.set_xlim(0)
			ax.set_ylim(0)
			ax.set_xticks(ind + (width / 2)) # ind+width
			ax.set_xticklabels((xlabels))
			xtickNames = ax.set_xticklabels(xlabels)
			plt.setp(xtickNames, rotation=0, fontsize=10)

			sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=Normalize(vmin=min(Data_Represented), vmax=max(Data_Represented)))
			sm._A = []
			plt.colorbar(sm)

			fig.savefig(p.output_fn + "-NoSingles-Contig-Data-Represented.png", dpi=600, bbox_inches='tight')

		if p.dr_fig:

			logger(logfile_fh, "Generating Total Data Represented Figure.")

			import matplotlib.cm as cm
			from matplotlib.colors import Normalize

			font = {'family': 'sans-serif',
					'weight': 'normal',
					'size': 10}
			lines = {'linewidth': 0.5}

			matplotlib.rc('font', **font)
			matplotlib.rc('lines', **lines)

			fig = plt.figure(figsize=(7, 4))
			ax = fig.add_subplot(1, 1, 1)

			Data = [value for (key, value) in new_parts.iteritems()]

			N = 15  # Change to 16 if including 16 cutoff range
			ind = np.arange(N)
			width = 0.8

			Total_Contigs = int(input_queries) + .0
			# Will create a list with values of the counts of each of the partitions (done above)
			Total_Contigs_In_Network_list = [value for value in Data if value >= 1]
			Total_Contigs_In_Network = sum(Total_Contigs_In_Network_list)

			data_represented_membership_size_cutoffs = [2, 3, 4, 5, 10, 15, 20, 30, 40, 50, 100, 200, 300, 400, 500]

			Data_Represented = []

			for cutoff in data_represented_membership_size_cutoffs:

				Data_prerepresented = []

				for community_sizes in Data:

					if community_sizes >= cutoff:
						Data_prerepresented.append(community_sizes)

				# Calculate % of data represented by the new communities that are above the cutoff
				to_append = (sum(Data_prerepresented) / (Total_Contigs + .0)) * 100

				Data_Represented.append(to_append)

			my_cmap = cm.get_cmap('jet')
			my_norm = Normalize(vmin=min(Data_Represented), vmax=max(Data_Represented))

			ax.bar(ind, Data_Represented, width, color=my_cmap(my_norm(Data_Represented)),
					align='edge', alpha=1, zorder=1)

			xlabels = ["{0}".format(x) for x in data_represented_membership_size_cutoffs]

			ax.set_xlabel('Members/Group')
			ax.set_ylabel('Contigs Represented (%)')
			ax.set_xlim(0)
			ax.set_ylim(0)
			ax.set_xticks(ind + (width / 2))  # ind+width
			ax.set_xticklabels(xlabels)
			xtickNames = ax.set_xticklabels(xlabels)
			plt.setp(xtickNames, rotation=0, fontsize=10)

			sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=Normalize(vmin=min(Data_Represented), vmax=max(Data_Represented)))
			sm._A = []
			plt.colorbar(sm)

			fig.savefig(p.output_fn + "-Total-Data-Represented.png", dpi=600, bbox_inches='tight')

logger(logfile_fh, "Program Complete.")
logfile_fh.close()