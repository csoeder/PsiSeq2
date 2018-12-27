parent_vcf=$1
offspring_vcf=$2
bed_out=$3


#vcftools --vcf "$parent_vcf" --diff "$offspring_vcf" --diff-site  --stdout | awk '{if(($4 == "B") && ($7 == $8))print $1, $2-1, $2, "NA", 1; }' | tr " " "\t" > "$bed_out".tmp
#vcftools --vcf "$parent_vcf" --diff "$offspring_vcf" --diff-site  --stdout | awk '{if(($4 == "1"))print $1, $2-1, $2, "NA", 0;}' | tr " " "\t" >> "$bed_out".tmp

#### "calls" variants even where there isn't sufficient info - got to pull the actually called varz 
#	only non-refs from the parent:
cat "$parent_vcf" | vcf2bed | grep -e "1/1" -e "0/1" -e "1/0" | cut -f 1-3 > "$parent_vcf".tmp.bed 
sleep 10
#	all callable from the offspring:
cat "$offspring_vcf" | vcf2bed | grep -e "1/1" -e "0/1" -e "1/0" -e "0/0"  | cut -f 1-3 > "$offspring_vcf".tmp.bed
sleep 10
#cat "$offspring_vcf" | vcf2bed | grep -e "1/1" -e "0/1" -e "1/0"  | cut -f 1-3 > "$offspring_vcf".tmp.bed
bedtools intersect -wa -a "$parent_vcf".tmp.bed  -b "$offspring_vcf".tmp.bed | awk '{print $0"\tNA\t1"}' > "$bed_out".tmp
sleep 10
bedtools intersect -v -wa -a "$parent_vcf".tmp.bed  -b "$offspring_vcf".tmp.bed | awk '{print $0"\tNA\t0"}' >> "$bed_out".tmp
sleep 10

bedtools sort -i "$bed_out".tmp > "$bed_out"
sleep 10

rm "$parent_vcf".tmp.bed "$offspring_vcf".tmp.bed "$bed_out".tmp
