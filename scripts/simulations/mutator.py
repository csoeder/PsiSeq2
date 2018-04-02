from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint, choice
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("reference_genome", help="reference genome to be mutated")
parser.add_argument("custom_genome", help="output genome, with SNPs added")
parser.add_argument("-m", "--mutation_rate", help="SNP rate, per bp; default 1%% ", type=float, default=0.01 )

args = parser.parse_args()



fasta_in = args.reference_genome
fasta_out = args.custom_genome
mut_rate = args.mutation_rate



#samtools faidx /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_melanogaster/dm6.fa chr2L > dm6.chr2L.fa

chroms_out = []

for std_chrom in SeqIO.parse(fasta_in, "fasta"):

	std_chrom = std_chrom.seq.tostring()
	std_id = std_chrom.id
	std_name = std_chrom.name

	snp_sites = []

	for i in range(1, int(len(std_chrom)*mut_rate)): #	calculate number of sites to perturb
		snp_sites.append(randint(1, len(std_chrom)))

	snp_sites = list(set(snp_sites))
	snp_sites.sort()
	snp_sites.reverse()

	nu_start=0
	nu_chrom=''
	while len(snp_sites) > 0:
		old_start = nu_start
		nu_start = snp_sites.pop()
		nukes = ["A","T","C","G"]
		try:
			nukes.remove(std_chrom[nu_start].upper())
		except ValueError:
			#got one of these; if the nucleotide in the standard chromosome isn't in the replacement list
			#eg, N
			#don't worry about it
			pass
		nu_chrom = "%s%s%s" % tuple([nu_chrom, std_chrom[old_start:nu_start-1], choice(nukes)])
	nu_chrom = "%s%s" % tuple([nu_chrom, std_chrom[nu_start:]])

	chroms_out.append(SeqRecord(Seq(nu_chrom), id=std_id, name=std_name))


SeqIO.write(chroms_out, fasta_out, "fasta")
