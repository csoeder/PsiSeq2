out_title=$1
in_fastq=$2
ref_genome=$3

bwa aln "$ref_genome" "$in_fastq" > "$out_title".sai
bwa samse "$ref_genome" "$out_title".sai "$in_fastq" | samtools view -Shb | samtools sort -o "$out_title".sort.bam -
samtools index "$out_title".sort.bam



