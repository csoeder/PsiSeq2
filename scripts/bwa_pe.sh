out_title=$1
r1_fastq=$2
r2_fastq=$3
ref_genome=$4

bwa aln "$ref_genome" "$r1_fastq" > "$out_title".R1.sai
bwa aln "$ref_genome" "$r2_fastq" > "$out_title".R2.sai
bwa sampe "$ref_genome" "$out_title".R1.sai "$out_title".R2.sai "$r1_fastq" "$r2_fastq"| samtools view -Shb | samtools sort -o "$out_title".sort.bam -
samtools index "$out_title".sort.bam



