import argparse
from string import upper
import re

parser = argparse.ArgumentParser()
parser.add_argument("hybrid_pileup", help="pileup file for hybrid")
parser.add_argument("parent_pileup", help="pileup file for parent")
parser.add_argument("output", help="file to write shared SNPs to")
parser.add_argument("-v", "--verbose", action="store_true", help="verbose reporting")
args = parser.parse_args()

indel_regex = re.compile('[\+\-][0-9]+[ACGTNacgtn]+')
mismatch_regex = re.compile('[ACGTacgt]')


def pileup_scanner(file_in):
	# 	Opens a pileup $file_in;
	# 	records all the places where a read base disagrees with the ref
	# 	Outputs as a nested dictionary
	# 	{contig_name: {contig_metadata;
	# 		{ postion: {site_ref, site_cov, {base counts} } }} }
	# 	for pileup specs:
	# 	https://en.wikipedia.org/wiki/Pileup_format#Column_5:_The_bases_string
	# 	or http://samtools.sourceforge.net/pileup.shtml (deprecated)
	mpile_file = open(file_in, "r")
	contig_dict = {}
	for line in (mpile_file):

		try:
			contig, position, ref_base, coverage, read_base, qual = line.split()
		except ValueError:  # zero-cov sites don't report a readbase/qual Dx
			contig, position, ref_base, coverage = line.split()
			read_base, qual = "N", "I"

		if contig not in contig_dict.keys():
			# intialize the dictionary for each new contig
			contig_dict[contig] = dict((
										('mini_counter', 0),
										('sum_coverage', 0),
										('mm_count', 0),
										('mini_counter', 0),
										('position_dict', {})))
		contig_dict[contig]['mini_counter'] += 1
		contig_dict[contig]['sum_coverage'] += int(coverage)
		if (
			int(coverage) > 1 and "*" not in read_base and not
			indel_regex.search(read_base) and upper(ref_base) != "N" and
			mismatch_regex.search(read_base)):
			# if cov >1 and read_base is not an indel and the reference base
			# is not null and read bases contain mismatches
			# ... then tally the bases at this site
			read_base_dict = {
								"A": read_base.upper().count("A"),
								"T": read_base.upper().count("T"),
								"C": read_base.upper().count("C"),
								"G": read_base.upper().count("G")}
			read_base_dict[upper(ref_base)] += (
												read_base.upper().count(",") +
												read_base.upper().count("."))
			#  incremenet when the read base is the reference base
			contig_dict[contig]['mm_count'] += 1  # increment mismatch count
			contig_dict[contig]['position_dict'][int(position)] = {'ref_base':upper(ref_base), 'pileup_cov':int(coverage), 'base_dict':read_base_dict} 

		if int(position) % 1000 == 0 and args.verbose:
			print "Contig %s: %s KiloBases scanned!" % tuple([contig, int(position)/1000])

	mpile_file.close()
	return contig_dict

def contig_report(contig_dict):
	for contig in contig_dict.keys():
		if contig_dict[contig]['mm_count'] > 0:
			mismatch_warn_string = '.'
		else:
			mismatch_warn_string = ' (only Ns in the pileup reference base column?)'
		print "contig %s had an average coverage depth of %s reads and a raw mismatch count of %s%s" % tuple([contig, contig_dict[contig]['sum_coverage']/contig_dict[contig]['mini_counter'],  contig_dict[contig]['mm_count'], mismatch_warn_string])


def mismatch_chooser(site_dict):
	#{'base_dict': {
	#	'A': 0, 'C': 0, 'T': 0, 'G': 28}, 'ref_base': 'A', 'pileup_cov': 28}}
	choice = 'N'
	meta = None # 	We may want to later report information on polymorphism
	for bass in site_dict['base_dict'].keys():
		if 	site_dict['base_dict'][bass]/float(site_dict['pileup_cov']) >= 0.9:
			choice = bass
	return choice, meta

def contig_dict_comparator(parent_dict, hybrid_dict):
	shared_contig_list = list(set(hybrid_dict.keys()).intersection(set(parent_dict.keys())))
	comparison_dict = dict.fromkeys(shared_contig_list, [])
	for contig in shared_contig_list:
		minicount = 0
		for parent_pos in parent_dict[contig]['position_dict'].keys():
			minicount += 1
			total_count = len(parent_dict[contig]['position_dict'].keys())
			if parent_pos in hybrid_dict[contig]['position_dict'].keys():
				# 	If the parent variant site is variant in the hybrid...
				parent_minidict = parent_dict[contig]['position_dict'][parent_pos]
				hybrid_minidict = hybrid_dict[contig]['position_dict'][parent_pos]
				hyb_var, hyb_meta = mismatch_chooser(hybrid_minidict)
				par_var, par_meta = mismatch_chooser(parent_minidict)
				if hybrid_minidict[parent_pos]['ref_base'] !=  parent_minidict['ref_base']:
					# If this happens... something, somewhere has gone terribly wrong x_x
					print "WARNING: reference sequences disagree on contig %s, position %s !!!" % tuple([contig, parent_pos])
				elif par_var == parent_minidict['ref_base']:
					#if the site isn't actually variable in the parent, ignore it
					pass					
				elif (hyb_var == par_var and hybrid_minidict['ref_base'] == parent_minidict['ref_base']):
					#if the site has the same variant in hybrid and parent, record it as parent-derived
					comparison_dict[contig].append([parent_pos, 1])
				else:
					#if the parent site is the wrong variant in the hybrid, it's not parent-derived
					comparison_dict[contig].append([parent_pos, 0])
			else: 
			# 	If the parent variant site isn't variant in the hybrid, the hybrid site isn't parent-derived.
				comparison_dict[contig].append([parent_pos, 0])
			if minicount % 10000 == 0 and args.verbose:
				print "%s parent mismatch sites investigated of %s!" % tuple([minicount, total_count])
		print "Contig %s of %s compared..." % tuple([shared_contig_list.index(contig), len(shared_contig_list)])
		print
	return comparison_dict

def comparison_writer(comp_dict, file_out):
	write_file = open(file_out, "w")

	for contig in comp_dict.keys():
		for coord in comp_dict[contig]:
			write_file.write("%s\t%s\t%s\n" % tuple([contig, coord[0], coord[1]]) )
	write_file.close()

print "gathering hybrid contigs...".upper()
hyb_contig_dict = pileup_scanner(args.hybrid_pileup)
contig_report(hyb_contig_dict)

print
print "gathering parent contigs...".upper()
par_contig_dict = pileup_scanner(args.parent_pileup)
contig_report(par_contig_dict)

print
print "comparing parent and offspring...".upper()
comparison = contig_dict_comparator(par_contig_dict, hyb_contig_dict)

print "writing comparison to a file....".upper()
comparison_writer(comparison, args.output)

print "DONE!"

