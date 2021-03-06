import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint, choice, seed

parser = argparse.ArgumentParser()
parser.add_argument("parent1", help="first parent")
parser.add_argument("parent2", help="second parent")
parser.add_argument("hybrid_fasta", help="hybrid FASTA out")
parser.add_argument("collage", help="BED file of hybrid ancestry out")
parser.add_argument("--introgression_number", help="number of introgression sites, total", type=int, default=10 )
parser.add_argument("--introgression_size", help="size of introgression sites, in bp", type=int, default=5000 )
parser.add_argument("-s", "--set_seed", help="set RNG seed; default  ", type=int, default=None )
args = parser.parse_args()

introg_num = args.introgression_number
introg_size = args.introgression_size
introg_points = []
seed(args.set_seed)

p1 = next(SeqIO.parse(args.parent1, "fasta"))
chrom_name = p1.id
p1 = p1.seq.tostring()
p2 = next(SeqIO.parse(args.parent2, "fasta")).seq.tostring()
chrom_len = len(p1)
# add multichrom function?


for i in range(0, introg_num):
	introg_points.append(randint(1,chrom_len))


introg_points = list(set(introg_points))
introg_points.sort()
introg_points.reverse()

nu_start = 0
nu_chrom=''
collage_file = open(args.collage, 'w')

while len(introg_points) > 0:
	old_start = nu_start
	nu_start = introg_points.pop()
	nu_chrom = "%s%s%s" % tuple([nu_chrom, p1[old_start:nu_start], p2[nu_start:nu_start+introg_size]])
	collage_file.write('%s\t%s\t%s\tP1\n' % tuple([chrom_name, old_start, nu_start]))
	collage_file.write('%s\t%s\t%s\tP2\n' % tuple([chrom_name, nu_start, nu_start+introg_size]) )
	nu_start += introg_size
nu_chrom = "%s%s" % tuple([nu_chrom, p1[nu_start:]])
collage_file.write('%s\t%s\t%s\tP1\n' % tuple([chrom_name, nu_start, chrom_len]))


record=SeqRecord(Seq(nu_chrom),
                   id=chrom_name, name=chrom_name)
#print len(nu_chrom)
#print len(p1)
#print len(p2)
SeqIO.write(record, args.hybrid_fasta, "fasta")