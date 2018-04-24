parent_vcf=$1
offspring_vcf=$2
bed_out=$3


vcftools --vcf "$parent_vcf" --diff "$offspring_vcf" --diff-site  --stdout | awk '{if(($4 == "B") && ($7 == $8))print $1, $2-1, $2, "NA", 1; }' | tr " " "\t" > "$bed_out".tmp
vcftools --vcf "$parent_vcf" --diff "$offspring_vcf" --diff-site  --stdout | awk '{if(($4 == "1"))print $1, $2-1, $2, "NA", 0;}' | tr " " "\t" >> "$bed_out".tmp

bedtools sort -i "$bed_out".tmp > "$bed_out"

rm "$bed_out".tmp



