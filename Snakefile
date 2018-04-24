

# borrowed from examples at https://bitbucket.org/holtgrewe/snakemake/

configfile: 'config.yaml'



#depends: bedtools samtools bwa vcftools freebayes


sample_by_name = {c['name'] : c for c in config['data_sets']}
ref_genome_by_name = { g['name'] : g for g in config['reference_genomes']}
chain_dict_by_destination = config['lift_genomes']
#sea_dubya_dee=config['cwd']
samps2process = [c['name'] for c in config['data_sets'] if c['pedigree'] == 'offspring' ]



def get_input_files(directory):
	#collect files
	result = []
	#directory = 'FASTQs/%s/%s/' % tuple([treatment, sample_name])
	for fname in os.listdir(directory):
		if fname.endswith('.fastq') or fname.endswith('.fq'):
			result.append(os.path.join(directory, fname))
	result = list(sorted(result))
	#sanity chex go here
	return result

def samps_by_group(grup):
	subset_out = []
	grupt = {}
	for c in config['data_sets']:
		try:
			grupt[c['name']] = c['group']
		except KeyError:
			pass
	for k in grupt.keys():
		if grupt[k] == grup:
			subset_out += [k]
	subset_out.sort()
	return subset_out

def return_group_by_samp(wildcards):
	try:
		return sample_by_name[wildcards.sample]['group']
	except KeyError:
		raise KeyError("%s has no group listed - edit the config file to add one" % tuple([sample]))


# rule all:
# 	input:
# 		expand("analysis_out/{sample}_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed", sample= [k for k,v in sample_by_name.items() if v['pedigree'] == 'offspring']),
# 		expand("analysis_out/{sample}_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed", sample= [k for k,v in sample_by_name.items() if v['pedigree'] == 'offspring'])

rule synthetic_reads_se:
	output: 
		reads='FASTQs/{treatment}/{sample}/{sample}.fq'
	params:
		runmem_gb=8,
		runtime="12:00:00",
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
		cov_depth=20,
		runmem_gb=8,
		runtime="12:00:00"
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
	direc = sample_by_name[wildcards.sample]['path']
	if sample_by_name[wildcards.sample]['type'] == 'synthetic':
		if not sample_by_name[wildcards.sample]['paired']:
			return expand("{directory}{sample}.fq", directory=[direc], sample=[wildcards.sample], treatment=sample_by_name[wildcards.sample]['treatment'] )
		else:
			return expand("{directory}{sample}_{readnum}.fq", directory=[direc], sample=[wildcards.sample], treatment=sample_by_name[wildcards.sample]['treatment'], readnum=[1,2])
	return get_input_files(direc)


rule bwa_sam:
	output:
		bam='mapped_reads/{sample}/{sample}_vs_{parent}.bwa.sort.bam',
	input:
		check_for_reads
	threads: 4
	params:
		runmem_gb=8,
		runtime="12:00:00"	
	run:
		read_files = get_input_files(sample_by_name[wildcards.sample]['path'])
		paired_reads = sample_by_name[wildcards.sample]['paired']
		ref_genome = ref_genome_by_name[wildcards.parent]['path']		
		if paired_reads:
			shell(
				"sh scripts/bwa_pe.sh mapped_reads/{wildcards.sample}/{wildcards.sample}_vs_{wildcards.parent}.bwa {read_files[0]} {read_files[1]} {ref_genome}"
				)
		else:
			shell(
				"sh scripts/bwa_se.sh mapped_reads/{wildcards.sample}/{wildcards.sample}_vs_{wildcards.parent}.bwa {read_files[0]} {ref_genome}"
				)

#https://www.biostars.org/p/56246/
rule bwa_uniqueUpOnIt:
	output:
		bam_out='mapped_reads/{sample}/{sample}_vs_{parent}.bwaUniq.sort.bam'
	input:
		bam_in='mapped_reads/{sample}/{sample}_vs_{parent}.bwa.sort.bam'
	params:
		quality="-q 20 -F 0x0100 -F 0x0200 -F 0x0300 -F 0x04",
		uniqueness="XT:A:U.*X0:i:1.*X1:i:0"
	threads: 4
	params:
		runmem_gb=8,
		runtime="3:00:00"
	run:
		ref_genome = ref_genome_by_name[wildcards.parent]['path']	
		shell('samtools view {params.quality} {input.bam_in} | grep -E {params.uniqueness} | samtools view -bS -T {ref_genome} - | samtools sort -o {output.bam_out} - ')
#		'samtools view {params.quality} mapped_reads/SucSec/SucSec_vs_droSim1.sort.bam | grep "XT:A:U" | grep  "X0:i:1" | grep "X1:i:0" samtools view -bS -T referenceSequence.fa - > reads.uniqueMap.bam'


rule NGM_sam:
	output:
		bam_out='mapped_reads/{sample}/{sample}_vs_{parent}.ngm.sort.bam'
	input:
		reads=check_for_reads
	params:
		strat=1,
		topn=1,
		runmem_gb=8,
		runtime="12:00:00"
	run:
		ref_genome = ref_genome_by_name[wildcards.parent]['path']
		paired_reads = sample_by_name[wildcards.sample]['paired']
		sam="%s.sam" % '.'.join(output.bam_out.rsplit('.')[:-2])
		if paired_reads:
			shell('/nas/longleaf/home/csoeder/modules/NextGenMap-0.5.0/bin/ngm-0.5.0/ngm -1 {input.reads[0]} -2 {input.reads[1]} -r {ref_genome} -o {sam} -b -n {params.topn} --no-unal --strata {params.strat}' )
		else:
			shell('/nas/longleaf/home/csoeder/modules/NextGenMap-0.5.0/bin/ngm-0.5.0/ngm -q {input.reads[0]}  -r {ref_genome} -o {sam} -b -n {params.topn} --no-unal --strata {params.strat}' )
		shell('samtools sort -o {output.bam_out} {sam}')


rule mpiler:
	input:
		sorted_bam='mapped_reads/{sample}/{sample}_vs_{parent}.{aligner}.sort.bam'
	output:
		mpile='mapped_reads/{sample}/{sample}_vs_{parent}.{aligner}.mpileup'
	params:
		runmem_gb=8,
		runtime="3:00:00"
	run:
		ref_genome = ref_genome_by_name[wildcards.parent]['path']		
		shell(
		"samtools mpileup -Bf {ref_genome} {input.sorted_bam} > {output.mpile}"
		)


rule shared_snipper:
	input:
		offspring_mpile='mapped_reads/{sample}/{sample}_vs_{parent}.{aligner}.mpileup',
		otherparent_mpile='mapped_reads/{compare}/{compare}_vs_{parent}.{aligner}.mpileup'
	output:
		shared_snps='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{aligner}.sharedSnps.out'
	params:
		minCov=5,
		minFrac=0.9,
		runmem_gb=32,
		runtime="48:00:00"
	shell:
		"python scripts/shared_snps.py -c {params.minCov} -f {params.minFrac} {input.offspring_mpile} {input.otherparent_mpile} {output.shared_snps}"		


rule out2bed:
	input:
		snps_out='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{aligner}.sharedSnps.out'
	output:
		snps_bed='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{aligner}.sharedSnps.bed'
	params:
		runmem_gb=8,
		runtime="3:00:00"
	shell:
		'sh scripts/out2bed.sh {input.snps_out} {output.snps_bed}'

rule lifter:
	input:
		unlifted='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{aligner}.sharedSnps.bed'
	output:
		lifted='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{aligner}.lift2{lift_genome}.sharedSnps.bed',
		too_heavy='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{aligner}.{lift_genome}.sharedSnps.tooHeavy'
	params:
		runmem_gb=8,
		runtime="3:00:00"
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
		heavy='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{aligner}.{lift_genome}.sharedSnps.tooHeavy'
	output:
		heavyBed='variant_comparisons/{sample}_and_{compare}_vs_{parent}.{aligner}.{lift_genome}.sharedSnps.tooHeavy.bed' 
	params:
		runmem_gb=8,
		runtime="3:00:00"
	shell:
		"cat {input.heavy} |  grep -v '#' > {output.heavyBed}"

rule window_maker:
	output:
		windowed='{ref_genome}_w{window_size}_s{slide_rate}.windows.bed'
	params:
		runmem_gb=8,
		runtime="1:00:00"
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
	params:
		runmem_gb=8,
		runtime="1:00:00"

	shell:
		'bedtools map -c 5,5 -o sum,count -null 0 -a {input.windows} -b {input.snps} > {output.window_counts}'



# vcf_subset_psiseq2 = ["10A", "SucSec", "SRR869587", "SRR6426002", "SRR5860570"]
# vcf_subset_psiseq2.sort()
# vcf_subset_nicole = ["PARC1", "PARG1", "REC7"]
# vcf_subset_nicole.sort()

ruleorder: RGfix_single > bwa_sam



rule RGfix_single:
	input:
		"mapped_reads/{sample}/{sample}_vs_{parent}.{aligner}.sort.bam"
	output:
		"mapped_reads/{sample}/{sample}_vs_{parent}.RG{id}.{aligner}.sort.bam"
	params:
		rg="RGLB=lib1 RGPL=illumina RGPU={sample} RGSM={sample} RGID={id}",
		runmem_gb=8,
		runtime="6:00:00"
	shell:
		"java -jar ~/modules/picard/build/libs/picard.jar AddOrReplaceReadGroups I={input} O={output} {params.rg} VALIDATION_STRINGENCY=LENIENT"


def demand_rgFix_by_group(wildcards):
	return [ "mapped_reads/%s/%s_vs_%s.RG%s.%s.sort.bam" % pear for pear in [ tuple([v, v, wildcards.ref_genome, str(samps_by_group(wildcards.group).index(v)), 'bwa']) for v in samps_by_group(wildcards.group) ] ]

rule RGfix_group:
	input:
#		psiseq2=["mapped_reads/%s/%s_vs_droSim1.RG%s.%s.sort.bam" % pear for pear in [ tuple([v, v, str(vcf_subset_psiseq2.index(v)), 'bwa']) for v in vcf_subset_psiseq2 ] ],
#		nicole=["mapped_reads/%s/%s_vs_dm6.RG%s.%s.sort.bam" % pear for pear in [ tuple([v, v, str(vcf_subset_nicole.index(v)), 'bwa']) for v in vcf_subset_nicole ] ]
#		[lambda wildcards: "mapped_reads/%s/%s_vs_%s.RG%s.%s.sort.bam" % pear for pear in [ tuple([v, v, wildcards.ref_genome, str(samps_by_group(wildcards.group).index(v)), 'bwa']) for v in samps_by_group(wildcards.group) ] ]
		potato = demand_rgFix_by_group,
	params:
		runmem_gb=8,
		runtime="6:00:00"
	output:
		"{group}_vs_{ref_genome}_RG_fixed.flag"
	shell:
		"touch {output}"

rule vcf_indiv:
	input:
		#bam = lambda wildcards: "mapped_reads/%s/%s_vs_%s.RG%s.{aligner}.sort.bam" % tuple([ wildcards.sample, wildcards.prefix, wildcards.ref_genome, str(vcf_subset.index(wildcards.sample)) ])
		bam = lambda wildcards: "mapped_reads/%s/%s_vs_%s.%s.sort.bam" % tuple([ wildcards.sample, wildcards.sample, wildcards.ref_genome, wildcards.aligner ])
	output:
		vcf_out = "variant_comparisons/{prefix}_vs_{ref_genome}.{aligner}.indiv.{sample}.vcf"
	params:
		freebayes="--standard-filters --min-coverage 4",
		runmem_gb=8,
		runtime="12:00:00"
	run:
		path2ref = ref_genome_by_name[wildcards.ref_genome]['path']
		shell("freebayes {params.freebayes} -f {path2ref} {input.bam} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.vcf_out}")


ruleorder: call_indiv_from_joint_vcf > vcf_joint


rule vcf_joint:
	input: 
		flag_in="{group}_vs_{ref_genome}_RG_fixed.flag"
	output:
		joint_vcf="variant_comparisons/{group}_vs_{ref_genome}.vcf"
#	cluster:--mem=16G -n 4 -t 2:00:00 
	params:
		freebayes="--standard-filters --min-coverage 4",
		runmem_gb=8,
		runtime="12:00:00"
	run:
		path2ref = ref_genome_by_name[wildcards.ref_genome]['path']
		run("freebayes {params.freebayes} -f {path2ref} {input} | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.joint_vcf}")

rule call_indiv_from_joint_vcf:
	input:
		lambda wildcards: "variant_comparisons/%s_vs_%s.vcf" % tuple([return_group_by_samp(wildcards), wildcards.ref_genome])
	output:
		vcf_out = "variant_comparisons/{prefix}_vs_{ref_genome}.{aligner}.joint.{sample}.vcf"	
	params:
		runmem_gb=8,
		runtime="1:00:00"
	shell:
		"vcftools --indv {wildcards.sample} --vcf {input} --recode --recode-INFO-all --stdout > {output.vcf_out}"


rule zyggy_stardust:
	input: 
		#vcf_in = "variant_comparisons/{vcf_prefix}.{call_type}.{sample}.vcf",
		vcf_in = "variant_comparisons/{vcf_prefix}.vcf",
		windoze = "{window_prefix}.windows.bed",
	output:
		#bed_out = "{vcf_prefix}.{window_prefix}.{who}.windowedZygosity.bed"
		bed_out = "analysis_out/{vcf_prefix}.{window_prefix}.windowedZygosity.bed"
	params:
		runmem_gb=8,
		runtime="6:00:00"
	shell:
		#"sh scripts/vcf_zyggy.sh {input.vcf_in} {input.windoze} {output.bed_out} {wildcards.who}"
		"sh scripts/vcf_zyggy.sh {input.vcf_in} {input.windoze} {output.bed_out}"

#rule take_off_every_zyg:



#https://www.biostars.org/p/178146/
# hets & homs treated as equivalent
rule naive_vcf_compare:
	input:
		vcf_parent = "variant_comparisons/{sample_p}_vs_{genome}.{mapper_p}.{caller_p}.{sample_p}.vcf",
		vcf_offspring = "variant_comparisons/{sample_o}_vs_{genome}.{mapper_o}.{caller_o}.{sample_o}.vcf",
#		vcf_parent = "variant_comparisons/{vcf_prefix_p}.{call_type_p}.{sample_p}.vcf",
#		vcf_offspring = "variant_comparisons/{vcf_prefix_o}.{call_type_o}.{sample_o}.vcf",
#		windoze = "{window_prefix}.windows.bed",
	params:
		runmem_gb=8,
		runtime="12:00:00"
	output:
		vcf_shared_snps_bed = "variant_comparisons/{sample_o}.{mapper_o}.{caller_o}.sharedWith.{sample_p}.{mapper_p}.{caller_p}.vs_{genome}.vcfNaive.sharedSnps.bed"
	shell:
		"sh scripts/vcf_naive.sh {input.vcf_parent} {input.vcf_offspring} {output.vcf_shared_snps_bed}"

# shared vars
#vcftools --vcf SRR5860570_vs_droSim1.filt.indiv.vcf --diff  SRR6426002_vs_droSim1.filt.indiv.vcf  --diff-site  --stdout | awk '{if(($4 == "B") && ($7 == $8))print $1, $2-1, $2;}' 
# unshared parent vars 
#vcftools --vcf {input.vcf_parent} --diff  {input.vcf_offspring}  --diff-site  --stdout | awk '{if(($4 == "1"))print $1, $2-1, $2, 0;}'



#leverage trio information
#rule nuanced_vcf_compare:

#http://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#bedops





rule nicole_flies_reanalysis:
	input:
		fly_finals = ["analysis_out/REC7_and_PARG1_vs_dm6.bwa.dm6_w100000_s10000.windowCounts.bed", "analysis_out/REC7_and_PARC1_vs_dm6.bwa.dm6_w100000_s10000.windowCounts.bed", "analysis_out/REC7.bwa.joint.sharedWith.PARC1.bwa.joint.vs_dm6.vcfNaive.dm6_w100000_s10000.windowCounts.bed", "analysis_out/REC7.bwa.joint.sharedWith.PARG1.bwa.joint.vs_dm6.vcfNaive.dm6_w100000_s10000.windowCounts.bed"]
	output:
		flg = "nicoleFlies.flag"
	params:
		runmem_gb=1,
		runtime="48:00:00"
	shell:
		"touch {output.flg}"



rule build_PsiSeq2_analysis_withRich:
	input:
		lift_tests = ["analysis_out/10A_and_SynthSec_vs_droSim1.droSim1_w100000_s10000.windowCounts.bed","analysis_out/10A_and_SynthSim_vs_droSec1.droSec1_w100000_s10000.windowCounts.bed","variant_comparisons/10A_and_SynthSec_vs_droSim1.dm6.sharedSnps.tooHeavy.bed"],
		psiSeq_lines = ["analysis_out/10A_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/10B_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/13A_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/13B_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/17A_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/17B_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR303333_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/10A_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/10B_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/13A_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/13B_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/17A_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/17B_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR303333_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed"],
		check_Flybase_genomes=['analysis_out/SucSec_and_SynthSecFB_vs_dsimr2.bwa.dsimr2_w100000_s10000.windowCounts.bed','analysis_out/SRR869587_and_SynthSecFB_vs_dsimr2.bwa.dsimr2_w100000_s10000.windowCounts.bed','analysis_out/SRR6426002_and_SynthSecFB_vs_dsimr2.bwa.dsimr2_w100000_s10000.windowCounts.bed','analysis_out/SRR5860570_and_SynthSecFB_vs_dsimr2.bwa.dsimr2_w100000_s10000.windowCounts.bed','analysis_out/SucSec_and_SynthSimFB_vs_dsecr1.bwa.dsecr1_w100000_s10000.windowCounts.bed','analysis_out/SRR869587_and_SynthSimFB_vs_dsecr1.bwa.dsecr1_w100000_s10000.windowCounts.bed','analysis_out/SRR6426002_and_SynthSimFB_vs_dsecr1.bwa.dsecr1_w100000_s10000.windowCounts.bed','analysis_out/SRR5860570_and_SynthSimFB_vs_dsecr1.bwa.dsecr1_w100000_s10000.windowCounts.bed'],
		control_sec_lines = ["analysis_out/SucSec_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SucSec_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR869587_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR869587_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR6426002_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR6426002_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR5860570_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR5860570_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed"],
		different_parent_lines=["analysis_out/10A_and_SRR869587_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/10A_and_SRR6426002_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/10A_and_SRR5860570_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SucSec_and_SRR869587_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SucSec_and_SRR6426002_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SucSec_and_SRR5860570_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed"],
		different_mapping_strategies = ["analysis_out/SucSec_and_SynthSec_vs_droSim1.bwa.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SucSec_and_SynthSec_vs_droSim1.bwaUniq.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SucSec_and_SynthSec_vs_droSim1.ngm.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SucSec_and_SynthSim_vs_droSec1.bwa.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SucSec_and_SynthSim_vs_droSec1.bwaUniq.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SucSec_and_SynthSim_vs_droSec1.ngm.lift2dm6.dm6_w100000_s10000.windowCounts.bed"],
		heterozygosity_results = ["analysis_out/10A_vs_droSim1.bwa.indiv.10A.droSim1_w100000_s10000.windowedZygosity.bed","analysis_out/10A_vs_droSim1.bwa.joint.10A.droSim1_w100000_s10000.windowedZygosity.bed","analysis_out/SRR5860570_vs_droSim1.bwa.indiv.SRR5860570.droSim1_w100000_s10000.windowedZygosity.bed","analysis_out/SRR5860570_vs_droSim1.bwa.joint.SRR5860570.droSim1_w100000_s10000.windowedZygosity.bed","analysis_out/SRR6426002_vs_droSim1.bwa.indiv.SRR6426002.droSim1_w100000_s10000.windowedZygosity.bed","analysis_out/SRR6426002_vs_droSim1.bwa.joint.SRR6426002.droSim1_w100000_s10000.windowedZygosity.bed","analysis_out/SRR869587_vs_droSim1.bwa.indiv.SRR869587.droSim1_w100000_s10000.windowedZygosity.bed","analysis_out/SRR869587_vs_droSim1.bwa.joint.SRR869587.droSim1_w100000_s10000.windowedZygosity.bed","analysis_out/SucSec_vs_droSim1.bwa.indiv.SucSec.droSim1_w100000_s10000.windowedZygosity.bed","analysis_out/SucSec_vs_droSim1.bwa.joint.SucSec.droSim1_w100000_s10000.windowedZygosity.bed"],
		variant_based_results = [],
	params:
		runmem_gb=1,
		runtime="56:00:00"

	output:
		pdf_report = "PsiSeq2.pdf"
	shell:
		"sh scripts/markerDown.sh PsiSeq2_TheSequel.Rmd {output.pdf_report}"
# 		'R -e Sys.setenv"(RSTUDIO_PANDOC="/nas/longleaf/apps/rstudio/1.0.136/bin/pandoc")" -e  rmarkdown::render"("scripts/PsiSeq2_TheSequel.Rmd",output_file="PsiSeq2.pdf")"'



#different_mapping_strategies=[],
# non-suc example?



rule build_PsiSeq2_analysis_withoutRich:
	input:
		lift_tests = ["analysis_out/10A_and_SynthSec_vs_droSim1.droSim1_w100000_s10000.windowCounts.bed","analysis_out/10A_and_SynthSim_vs_droSec1.droSec1_w100000_s10000.windowCounts.bed","variant_comparisons/10A_and_SynthSec_vs_droSim1.dm6.sharedSnps.tooHeavy.bed"],
		psiSeq_lines = ["analysis_out/10A_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/10B_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/13A_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/13B_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/17A_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/17B_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR303333_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/10A_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/10B_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/13A_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/13B_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/17A_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/17B_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR303333_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed"],
		check_Flybase_genomes=['analysis_out/SRR869587_and_SynthSecFB_vs_dsimr2.bwa.dsimr2_w100000_s10000.windowCounts.bed','analysis_out/SRR6426002_and_SynthSecFB_vs_dsimr2.bwa.dsimr2_w100000_s10000.windowCounts.bed','analysis_out/SRR5860570_and_SynthSecFB_vs_dsimr2.bwa.dsimr2_w100000_s10000.windowCounts.bed','analysis_out/SRR869587_and_SynthSimFB_vs_dsecr1.bwa.dsecr1_w100000_s10000.windowCounts.bed','analysis_out/SRR6426002_and_SynthSimFB_vs_dsecr1.bwa.dsecr1_w100000_s10000.windowCounts.bed','analysis_out/SRR5860570_and_SynthSimFB_vs_dsecr1.bwa.dsecr1_w100000_s10000.windowCounts.bed'],
		control_sec_lines = ["analysis_out/SRR869587_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR869587_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR6426002_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR6426002_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR5860570_and_SynthSec_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/SRR5860570_and_SynthSim_vs_droSec1.lift2dm6.dm6_w100000_s10000.windowCounts.bed"],
		different_parent_lines=["analysis_out/10A_and_SRR869587_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/10A_and_SRR6426002_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed","analysis_out/10A_and_SRR5860570_vs_droSim1.lift2dm6.dm6_w100000_s10000.windowCounts.bed"],
		heterozygosity_results = [],
		variant_based_results = [],

	output:
		pdf_report = "PsiSeq2.noRich.pdf"
	shell:
		"sh markerDown.sh PsiSeq2_TheSequel.Rmd {output.pdf_report}"

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

#vcftools --vcf SRR5860570_vs_droSim1.filt.indiv.vcf --diff  SRR6426002_vs_droSim1.filt.indiv.vcf  --diff-site  --stdout | awk '{if($4 == "B")print $1, $2, $2+1;}' |tr ' ' '\t' |head -n 1 > test.bed



#sbatch --mem=16G -n 4 -t 2:00:00 -J bam2vcf_SRR869587 -o bam2vcf_SRR869587.slurm.out --wrap='freebayes --standard-filters --min-coverage 4 -f /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_simulans/droSim1.fa SRR869587_vs_droSim1.sort.bam | vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout | ~/modules/vcflib/bin/vcffilter -g "GT = 1|1" > SRR869587_vs_droSim1.filt.vcf'

#for ess in $(ls | grep SRR); do 
#> cat "$ess"/"$ess"_vs_droSim1.filt.vcf | vcf2bed > "$ess"/"$ess"_vs_droSim1.filt.vcf.bed; 
#> echo $ess;
#> done

#bedtools intersect -a mapped_reads/SRR5860570/SRR5860570_vs_droSim1.filt.vcf.bed -b mapped_reads/SRR6426002/SRR6426002_vs_droSim1filt.vcf.bed | bedtools intersect -a - -b mapped_reads/SRR869587/SRR869587_vs_droSim1.filt.vcf.bed | cut -f 1-7 > intersection_of_all_SRR_vs_droSim1.vcf.bed



#sbatch -n 3 -t 6:00:00 -o freebayes_everything.slurm.out "freebayes --standard-filters --min-coverage 4 -f /proj/cdjones_lab/Genomics_Data_Commons/genomes/drosophila_simulans/droSim1.fa mapped_reads/C*/*sort.bam mapped_reads/PSI*/*sort.bam mapped^C vcftools --remove-indels --vcf - --recode --recode-INFO-all --stdout  > {output.joint_vcf}"


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

