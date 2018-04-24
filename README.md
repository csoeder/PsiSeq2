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

snakemake --restart-times 6 --latency-wait 60 --jobs 24 -p --cluster "sbatch --time=24:00:00 -n 4 --mem=32G "


# Installation, Dependencies

bwa
samtools
bedtools
ngm
art
picard
vcftools
bedops
UCSC LiftOver
freebayes
R with markdown dependencies


# Components, Default Structure
```
├── .gitignore
├── Build_PsiSeq2_Results.snakefile 
		Rules file
├── FASTQs/ 
		a folder containing the sequenced reads to be analyzed
├── README.md
├── #LICENSE.md
├── PsiSeq2.yaml
		snakemake config file for standard PsiSeq2 analysis
├── PsiSeq2.noRich.yaml
		snakemake config file for standard PsiSeq2 analysis, without Rich's data
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

Analysis of simulans/sechellia introgression, including Rich's sugarflies:

```
snakemake --snakefile Build_PsiSeq2_Results.snakefile --config PsiSeq2.yaml --cluster "sbatch --time={params.runtime} -n 4 --mem={params.runmem_gb}G "  -p PsiSeq2.pdf 
```



* crossover simulation as a default use case example/proof of concept
	* also demonstrates further utility for the crossover/gene conversion applications?

