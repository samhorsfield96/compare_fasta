import networkx as nx
import numpy as np
import scipy
from scipy import stats
import argparse


def get_options():
	parser = argparse.ArgumentParser(description='Print cluster size information from gml file', prog='python parse_gml.py')

	# input options
	parser.add_argument('--infile',
						required=True,
						help='.gml file to analysis')
	parser.add_argument('--outpref',
						default="parsed_graph",
						help='Output prefix. Default = parsed_graph')

	return parser.parse_args()

def get_node_stats(G):
	stat_list_IQR = []
	stat_list_prop_IQR = []
	stat_list_MAD = []
	stat_list_prop_MAD = []
	stat_list_prop_stdev = []
	stat_list_stdev = []
	stat_size_cluster = []
	for n in G.nodes():
		node_info = G.nodes[n]
		ORF_sizes = node_info['lengths']
		ORF_sizes = np.array(ORF_sizes).astype(int)
		stat_size_cluster.append(ORF_sizes.size)
		q75, q50, q25 = np.percentile(ORF_sizes, [75, 50, 25])

		if ORF_sizes.size == 1:
			stat_list_MAD.append(0.0)
			stat_list_prop_MAD.append(0.0)
		else:
			test = scipy.stats.median_abs_deviation(ORF_sizes)
			stat_list_MAD.append(scipy.stats.median_abs_deviation(ORF_sizes))
			stat_list_prop_MAD.append(scipy.stats.median_abs_deviation(ORF_sizes) / q50)

		stat_list_IQR.append(float(q75 - q25))
		stat_list_prop_IQR.append(float((q75 - q25) / q50))

		stat_list_prop_stdev.append(float(np.std(ORF_sizes)) / float(np.mean(ORF_sizes)))
		stat_list_stdev.append(float(np.std(ORF_sizes)))

	stat_list_IQR = np.array(stat_list_IQR)
	stat_list_prop_IQR = np.array(stat_list_prop_IQR)

	stat_list_MAD = np.array(stat_list_MAD)
	stat_list_prop_MAD = np.array(stat_list_prop_MAD)

	stat_list_stdev = np.array(stat_list_stdev)
	stat_list_prop_stdev = np.array(stat_list_prop_stdev)

	stat_size_cluster = np.array(stat_size_cluster)

	return stat_list_IQR, stat_list_prop_IQR, stat_list_MAD, stat_list_prop_MAD, stat_list_prop_stdev, stat_list_stdev, stat_size_cluster


if __name__ == "__main__":
	options = get_options()
	infile = options.infile
	outpref = options.outpref

	G = nx.read_gml(infile)

	stat_list_IQR, stat_list_prop_IQR, stat_list_MAD, stat_list_prop_MAD, stat_list_prop_stdev, stat_list_stdev, stat_size_cluster = get_node_stats(G)

	combined_stdev = np.column_stack((stat_list_prop_stdev, stat_list_stdev))
	combined_IQR = np.column_stack((stat_list_prop_IQR, stat_list_IQR))
	combined_MAD = np.column_stack((stat_list_prop_MAD, stat_list_MAD))

	np.savetxt(outpref + "_stddev.txt", combined_stdev, delimiter=",")
	np.savetxt(outpref + "_IQR.txt", combined_IQR, delimiter=",")
	np.savetxt(outpref + "_MAD.txt", combined_MAD, delimiter=",")
	np.savetxt(outpref + "_sizes.txt", stat_size_cluster, delimiter=",")

	print("Average prop. standard deviation: {}".format(np.mean(stat_list_prop_stdev)))
	print("Stdev prop. standard deviation: {}".format(np.std(stat_list_prop_stdev)))
	print("Average standard deviation: {}".format(np.mean(stat_list_stdev)))
	print("Stdev standard deviation: {}".format(np.std(stat_list_stdev)))

	print("Average prop. IQR: {}".format(np.mean(stat_list_prop_IQR)))
	print("Stdev prop. IQR: {}".format(np.std(stat_list_prop_IQR)))
	print("Average IQR: {}".format(np.mean(stat_list_IQR)))
	print("Stdev IQR: {}".format(np.std(stat_list_IQR)))

	print("Average prop. MAD: {}".format(np.mean(stat_list_prop_MAD)))
	print("Stdev prop. MAD: {}".format(np.std(stat_list_prop_MAD)))
	print("Average MAD: {}".format(np.mean(stat_list_MAD)))
	print("Stdev MAD: {}".format(np.std(stat_list_MAD)))