

# borrowed from examples at https://bitbucket.org/holtgrewe/snakemake/

configfile: 'config.yaml'



#depends: bedtools samtools bwa vcftools freebayes


sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}
chain_dict_by_destination = config['lift_genomes']
samps2process = [c['name'] for c in config['data_sets'] if c['pedigree'] == 'offspring' ]


def get_input_files(sample_name, treatment):
	#collect files
	result = []
	directory = 'FASTQs/%s/%s/' % tuple([treatment, sample_name])
	for fname in os.listdir(directory):
		if fname.endswith('.fastq') or fname.endswith('.fq'):
			result.append(os.path.join(directory, fname))
	result = list(sorted(result))
	#sanity chex go here
	return result


rule all:
	input:
		expand("analysis_out/{sample}_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed", sample= [k for k,v in sample_by_name.items() if v['pedigree'] == 'offspring']),
		expand("analysis_out/{sample}_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed", sample= [k for k,v in sample_by_name.items() if v['pedigree'] == 'offspring'])

rule synthetic_reads_se:
	output: 
		reads='FASTQs/{treatment}/{sample}/{sample}.fq'
	params:
		platform='HS25',
		read_len=100,
		cov_depth=20
	run:
		ref_genome = ref_genome_by_name[sample_by_name[wildcards.sample]['pedigree']]['path']
		if not sample_by_name[wildcards.sample]['paired']:
			shell(
				"~/modules/art_bin_MountRainier/art_illumina -na -ss '{params.platform}' -i {ref_genome} -l '{params.read_len}' -f '{params.cov_depth}' -o FASTQs/{wildcards.treatment}/{wildcards.sample}/{wildcards.sample}"
				)
		else:
			shell("echo 'THESE READS ARE SUPPOSED TO BE PAIRED'")


rule synthetic_reads_pe:
	output: 
		reads1='FASTQs/{treatment}/{sample}/{sample}_1.fq',
		reads2='FASTQs/{treatment}/{sample}/{sample}_2.fq'
	params:
		platform='HS25',
		read_len=100,
		insert_mean=200,
		insert_sd=10,
		cov_depth=20
	run:
		ref_genome = ref_genome_by_name[sample_by_name[wildcards.sample]['pedigree']]['path']
		if sample_by_name[wildcards.sample]['paired']:
			shell(
				"~/modules/art_bin_MountRainier/art_illumina -p -na -ss '{params.platform}'  -i {ref_genome} -l '{params.read_len}' -f '{params.cov_depth}' -m '{params.insert_mean}' -s '{params.insert_sd}' -o FASTQs/{wildcards.treatment}/{wildcards.sample}/{wildcards.sample}_"
				)
		else:
			shell("echo 'THESE READS ARE NOT SUPPOSED TO PAIRED'")

ruleorder: synthetic_reads_pe > synthetic_reads_se

#https://bitbucket.org/snakemake/snakemake/issues/37/add-complex-conditional-file-dependency
def check_for_reads(wildcards):
	if sample_by_name[wildcards.sample]['type'] == 'synthetic':
		if not sample_by_name[wildcards.sample]['paired']:
			return expand("FASTQs/{treatment}/{sample}/{sample}.fq", sample=[wildcards.sample], treatment=sample_by_name[wildcards.sample]['treatment'] )
		else:
			return expand("FASTQs/{treatment}/{sample}/{sample}_{readnum}.fq", sample=[wildcards.sample], treatment=sample_by_name[wildcards.sample]['treatment'], readnum=[1,2])
	return get_input_files(wildcards.sample, sample_by_name[wildcards.sample]['treatment'])


rule bwa_sam:
	output:
		bam='mapped_reads/{sample}/{sample}_vs_{parent}.sort.bam',
	input:
		check_for_reads
	run:
		read_files = get_input_files(wildcards.sample, sample_by_name[wildcards.sample]['treatment'])
		paired_reads = sample_by_name[wildcards.sample]['paired']
		ref_genome = ref_genome_by_name[wildcards.parent]['path']		
		if paired_reads:
			shell(
				"sh scripts/bwa_pe.sh mapped_reads/{wildcards.sample}/{wildcards.sample}_vs_{wildcards.parent} {read_files[0]} {read_files[1]} {ref_genome}"
				)
		else:
			shell(
				"sh scripts/bwa_se.sh mapped_reads/{wildcards.sample}/{wildcards.sample}_vs_{wildcards.parent} {read_files[0]} {ref_genome}"
				)


rule mpiler:
	input:
		sorted_bam='mapped_reads/{sample}/{sample}_vs_{parent}.sort.bam'
	output:
		mpile='mapped_reads/{sample}/{sample}_vs_{parent}.mpileup'
	run:
		ref_genome = ref_genome_by_name[wildcards.parent]['path']		
		shell(
		"samtools mpileup -Bf {ref_genome} mapped_reads/{wildcards.sample}/{wildcards.sample}_vs_{wildcards.parent}.sort.bam > mapped_reads/{wildcards.sample}/{wildcards.sample}_vs_{wildcards.parent}.mpileup"
		)


rule shared_snipper:
	input:
		offspring_mpile='mapped_reads/{sample}/{sample}_vs_{parent}.mpileup',
		otherparent_mpile='mapped_reads/{compare}/{compare}_vs_{parent}.mpileup'
	output:
		shared_snps='variant_comparisons/{sample}_and_{compare}_vs_{parent}.sharedSnps.out'
	params:
		minCov=5,
		minFrac=0.9
	shell:
		"python scripts/shared_snps.py -c {params.minCov} -f {params.minFrac} {input.offspring_mpile} {input.otherparent_mpile} {output.shared_snps}"		

rule out2bed:
	input:
		snps_out='variant_comparisons/{sample}_and_{compare}_vs_{parent}.sharedSnps.out'
	output:
		snps_bed='variant_comparisons/{sample}_and_{compare}_vs_{parent}.sharedSnps.bed'
	shell:
		'sh scripts/out2bed.sh {input.snps_out} {output.snps_bed}'

rule lifter:
	input:
		unlifted='variant_comparisons/{sample}_and_{compare}_vs_{parent}.sharedSnps.bed'
	output:
		lifted='variant_comparisons/{sample}_and_{compare}_vs_{parent}.lift2{lift_genome}.sharedSnps.bed',
		too_heavy='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{lift_genome}.sharedSnps.tooHeavy'
	run:
		chain = chain_dict_by_destination[wildcards.lift_genome][wildcards.parent]
		shell(
			'~/modules/UCSC_utils/liftOver {input.unlifted} {chain} {output.lifted}.tmp {output.too_heavy}'
		)
		shell(
			'bedtools sort -i {output.lifted}.tmp > {output.lifted}'
		)
		shell(
			'rm {output.lifted}.tmp'
		)

ruleorder: lifter > out2bed

rule heavy2bed:
	input:
		heavy='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{lift_genome}.sharedSnps.tooHeavy'
	output:
		heavyBed='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{lift_genome}.sharedSnps.tooHeavy.bed'
	shell:
		"cat {input.heavy} |  grep -v '#' > {output.heavyBed}"

rule window_maker:
	output:
		windowed='{ref_genome}_w{window_size}_s{slide_rate}.windows.bed'
	run:
		fai_path = ref_genome_by_name[wildcards.ref_genome]['fai'],
		shell(
			'bedtools makewindows -w {wildcards.window_size} -s {wildcards.slide_rate} -g {fai_path} | bedtools sort -i - > {output.windowed}'
		)

rule window_counter:
	input:
		windows='{window_prefix}.windows.bed',
		snps='variant_comparisons/{snp_prefix}.sharedSnps.bed'
	output:
		window_counts='analysis_out/{snp_prefix}.{window_prefix}.windowCounts.bed'
	shell:
		'bedtools map -c 5,5 -o sum,count -null 0 -a {input.windows} -b {input.snps} > {output.window_counts}'



#rule generate_report:
#	input: expand("analysis_out/{samp}_and_{comp}_vs_{ref}.lift2dm6.dm6_w100000_s10000.windowCounts.bed", samp=samps2process, comp=['synthSec','synthSim'], ref=['droSim1', 'droSec1'])



vcf_subset = ["10A", "SucSec", "SRR869587", "SRR6426002", "SRR5860570"]
vcf_subset.sort()

ruleorder: RGfix_single > bwa_sam

# rule RGfix_single:
# 	#http://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/addorreplacereadgroups.html
#     input:
#         "mapped_reads/{sample}/{sample}_vs_{parent}.sort.bam"
#     output:
#         "mapped_reads/{sample}/{sample}_vs_{parent}.RG{id}.sort.bam"
#     log:
#         "logs/picard/replace_rg/{sample}.log"
#     params:
#         "RGLB=lib1 RGPL=illumina RGPU={sample} RGSM={sample} RGID={id}"
#     wrapper:
#         "0.22.0/bio/picard/addorreplacereadgroups"

rule RGfix_single:
	input:
		"mapped_reads/{sample}/{sample}_vs_{parent}.sort.bam"
	output:
		"mapped_reads/{sample}/{sample}_vs_{parent}.RG{id}.sort.bam"
	params:
		rg="RGLB=lib1 RGPL=illumina RGPU={sample} RGSM={sample} RGID={id}"
	shell:
		"java -jar ~/modules/picard/build/libs/picard.jar AddOrReplaceReadGroups I={input} O={output} {params.rg} VALIDATION_STRINGENCY=LENIENT"

rule RGfix_all:
	input:
		["mapped_reads/%s/%s_vs_droSim1.RG%s.sort.bam" % pear for pear in [ tuple([v, v, str(vcf_subset.index(v))]) for v in vcf_subset ] ]
	output:
		"all_RG_fixed.flag"
	shell:
		"touch {output}"


rule vcf_indiv:
	input:
		bam = lambda wildcards: "mapped_reads/%s/%s.RG%s.sort.bam" % tuple([ wildcards.sample, wildcards.prefix, str(vcf_subset.index(wildcards.sample)) ])
	output:
		vcf_out = "variant_comparisons/{prefix}.indiv.{sample}.vcf"
	params:
		freebayes="--standard-filters --min-coverage 4"
	shell:
		"freebayes {params.freebayes} -f /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_simulans/droSim1.fa {input.bam} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.vcf_out}"


rule vcf_joint:
	input: rules.RGfix_all.input
	output:
		joint_vcf="variant_comparisons/allTogetherNow_vs_droSim1.vcf"
#	cluster:--mem=16G -n 4 -t 2:00:00 
	params:
		freebayes="--standard-filters --min-coverage 4"
	shell:
		"freebayes {params.freebayes} -f /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_simulans/droSim1.fa {input} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.joint_vcf}"

rule call_indiv_from_joint_vcf:
	input: rules.vcf_joint.output
	output:
		vcf_out = "variant_comparisons/{prefix}.joint.{sample}.vcf"	
	shell:
		"vcftools --indv {wildcards.sample} --vcf {input} --recode --recode-INFO-all --stdout > {output.vcf_out}"



# 	# run:
# 	# 	"sh scripts/vcf_zyggy.sh {input.vcf_in} {input.windoze} {output.bed_out} {wildcards.who}"


# rule indiv_zyggy:
# 	input:
# 		vcf_in="variant_comparisons/{vcf_prefix}.indiv.{who}.vcf",
# 		windoze="{windows_prefix}.windows.bed"
# #		windoze="droSim1_w100000_s10000.windows.bed"
# 	output:
# 		bed_out="{vcf_prefix}.{windows_prefix}.windowCounts.zygosity.bed"
# 	run:
# 		"sh scripts/vcf_zyggy.sh {input.vcf_in} {input.windoze} {output.bed_out} {wildcards.who}"




# 	for i in indivs:
# 		cat  all_vs_droSim1.filt.joint.vcf | vcftools --indv "$i"  --remove-indels --vcf - --recode --recode-INFO-all --stdout  | vcf2bed | grep  '1/1' | bedtools map -c 2 -o count -a droSim1_w100000_s10000.windows.bed -b - | awk '{print$0"\tJOINT\tHOM"}' >> "$i".vsSim.droSim1_w100000_s10000.windowCounts.bed;
# 		cat  all_vs_droSim1.filt.joint.vcf | vcftools --indv "$i"  --remove-indels --vcf - --recode --recode-INFO-all --stdout  | vcf2bed | grep  '0/1' | bedtools map -c 2 -o count -a droSim1_w100000_s10000.windows.bed -b - | awk '{print$0"\tJOINT\tHET"}' >> "$i".vsSim.droSim1_w100000_s10000.windowCounts.bed;


# rule compare_VCFs:
# 	vcftools --vcf (parent).vcf --diff (offspring).vcf --diff-site   --out vcfToolsTest

# 	cat vcfToolsTest.diff.sites_in_files  | if $4 == "B" | if $7==$8  > vcf_shared_SNPS




#sbatch --mem=16G -n 4 -t 2:00:00 -J bam2vcf_SRR869587 -o bam2vcf_SRR869587.slurm.out --wrap='freebayes --standard-filters --min-coverage 4 -f /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_simulans/droSim1.fa SRR869587_vs_droSim1.sort.bam | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout | ~/modules/vcflib/bin/vcffilter -g "GT = 1|1" > SRR869587_vs_droSim1.filt.vcf'

#for ess in $(ls | grep SRR); do 
#> cat "$ess"/"$ess"_vs_droSim1.filt.vcf | vcf2bed > "$ess"/"$ess"_vs_droSim1.filt.vcf.bed; 
#> echo $ess;
#> done

#bedtools intersect -a mapped_reads/SRR5860570/SRR5860570_vs_droSim1.filt.vcf.bed -b mapped_reads/SRR6426002/SRR6426002_vs_droSim1filt.vcf.bed | bedtools intersect -a - -b mapped_reads/SRR869587/SRR869587_vs_droSim1.filt.vcf.bed | cut -f 1-7 > intersection_of_all_SRR_vs_droSim1.vcf.bed



#sbatch -n 3 -t 6:00:00 -o freebayes_everything.slurm.out "freebayes --standard-filters --min-coverage 4 -f /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_simulans/droSim1.fa mapped_reads/C*/*sort.bam mapped_reads/PSI*/*sort.bam mapped^C vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.joint_vcf}"




