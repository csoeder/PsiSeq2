from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random import randint, choice


#samtools faidx /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_melanogaster/dm6.fa chr2L > dm6.chr2L.fa

mut_rate=0.001 # 	SNP frequency, per bp

std_chrom = SeqIO.parse("dm6.chr2L.fa", "fasta").next().seq.tostring()

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
	nukes.remove(std_chrom[nu_start].upper())
	nu_chrom = "%s%s%s" % tuple([nu_chrom, std_chrom[old_start:nu_start-1], choice(nukes)])
nu_chrom = "%s%s" % tuple([nu_chrom, std_chrom[nu_start:]])


record=SeqRecord(Seq(nu_chrom),
                   id="chr2L", name="chr2L")
print len(nu_chrom)
print len(std_chrom)

SeqIO.write(record, "parent.fasta", "fasta")
