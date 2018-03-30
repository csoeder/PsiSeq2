# PsiSeq2
An update to PSISeq (Earley &amp; Jones 2011)


# Installation, Dependencies

# Components, Default Structure

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
│   ├── simulations/
			scripts for proof-of-concept simulation
│   │   ├── mutator.py
│   └── shared_snps.py
			compares two pileups for simliarities in one seen in the other
├── #envs
│   └── #myenv.yaml
└── Snakefile
		pipeline rules and processes

# Use, Default Process



