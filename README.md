# PsiSeq2
An update to PSISeq (Earley &amp; Jones 2011)


Software for mapping of complex traits through selection, introgression, and sequencing.

inputs:
	reference parent genome(s), FASTA
	sequenced parent DNA, FASTQ
		(optional, otherwise may be generated from reference genome)
	sequenced offpsring DNA, FASTQ

outputs:
	Analysis based on Jones & Early 2011 reimplementation
	Analysis based on variant-callers?
	R markdown summary




# Installation, Dependencies

bwa
samtools
bedtools
ngm
art
picard
vcftools
bedops

# Components, Default Structure
```
├── .gitignore
├── FASTQs/ 
		a folder containing the sequenced reads to be analyzed
├── README.md
├── #LICENSE.md
├── config.yaml
		snakemake config file, used to define pipeline parameters
├── scripts/
│   ├── bwa_pe.sh
			paired-end BWA aligner script
│   ├── bwa_se.sh
			single-end BWA aligner script
│   ├── out2bed.sh
			converts shared_snp comparison output to a BED file
│   ├── shared_snps.py
			compares two pileups for simliarities in one seen in the other
│   ├── simulations/
			scripts for proof-of-concept simulation
│   │   ├── mutator.py
│   └── ....................................
├── #envs
│   └── #myenv.yaml
└── Snakefile
		pipeline rules and processes
```
# Use, Default Process

* crossover simulation as a default use case example/proof of concept
	* also demonstrates further utility for the crossover/gene conversion applications?

