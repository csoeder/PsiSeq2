vcf_in=$1
windows=$2
bed_out=$3
name=$4

cat "$vcf_in" | vcftools --indv "$name"  --vcf - --recode --recode-INFO-all --stdout | vcf2bed | grep  '1/1' | bedtools map -c 2 -o count -a "$windows" -b - | awk '{print$0"\tHOM"}' > "$bed_out";
cat "$vcf_in" | vcftools --indv "$name"  --vcf - --recode --recode-INFO-all --stdout | vcf2bed | grep  '0/1' | bedtools map -c 2 -o count -a "$windows" -b - | awk '{print$0"\tHET"}' >> "$bed_out";

