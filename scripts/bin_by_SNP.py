import argparse
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument("file_in", help="input file")
parser.add_argument("file_out", help="output file")
parser.add_argument("bin_size", help='snps per bin')
parser.add_argument("slide_rate", help='bin slide rate')
args = parser.parse_args()

phial_in = args.file_in
phial_out = args.file_out
bin_size = int(args.bin_size)
slide_rate = int(args.slide_rate)



def parse_Snps(phial_in):
	comparison_dict = defaultdict(list)
	parser = open(phial_in, 'r')
	for line in parser.readlines():
		split_line = line.split('\t')
		comparison_dict[split_line[0]].append([int(split_line[1]), float(split_line[3])])
	parser.close()
	return comparison_dict


#modular window reporting function
def window_aggregator(window_in):
	return sum([w[1] for w in window_in])

def window_slider(comparison_dict):
	windows_list = []
	contig_list = comparison_dict.keys()
	for contig in contig_list:
		snps = sorted(comparison_dict[contig])
		for i in range(0,len(snps)-bin_size+slide_rate,slide_rate):
			window = snps[i:i+bin_size]
			start = window[0][0]
			end = window[-1][0]
			stat = window_aggregator(window)
			windows_list.append([contig, start, end, stat])
	return windows_list


def file_writer(win, phial):
	output = open(phial, 'w')
	for lyne in win:
		output.write('%s\t%s\t%s\t%s\n' % tuple(lyne))
	output.close()

comparison_dict = parse_Snps(phial_in)
windows = window_slider(comparison_dict)
file_writer(windows, phial_out)

