
sharedSnps1=$(echo $1 | rev | cut -f 1 -d "/" | rev)
sharedSnps2=$(echo $2 | rev | cut -f 1 -d "/" | rev)

#sharedSnps1=variant_comparisons/17A.bwaUniq.joint.sharedWith.SRR869587.bwaUniq.joint.vs_droSec1.lift2dm6.vcfNaive.sharedSnps.bed
#sharedSnps2=variant_comparisons/17A.bwaUniq.joint.sharedWith.SRR869587.bwaUniq.joint.vs_droSim1.lift2dm6.vcfNaive.sharedSnps.bed



# just ones
#cat "$sharedSnps1" | awk '{if($5>0)print;}' > test.bed
#cat "$sharedSnps2" | awk '{if($5<1)print$1"\t"$2"\t"$3"\t"$4"\t0";}' >> test.bed
#bedtools sort -i test.bed > test.sort.bed
#bedtools sort -i test.bed | uniq > test.sort.bed


#with zeros
cat variant_comparisons/"$sharedSnps1" > tmp/"$sharedSnps1"vs"$sharedSnps2".bed
cat variant_comparisons/"$sharedSnps2" | awk '{if($5<1)print$1"\t"$2"\t"$3"\t"$4"\t1";}' >> tmp/"$sharedSnps1"vs"$sharedSnps2".bed
cat variant_comparisons/"$sharedSnps2" | awk '{if($5>0)print$1"\t"$2"\t"$3"\t"$4"\t0";}' >> tmp/"$sharedSnps1"vs"$sharedSnps2".bed
cat tmp/"$sharedSnps1"vs"$sharedSnps2".bed | bedtools sort -i - | uniq | cut -f 1,2,3 | uniq -d > tmp/"$sharedSnps1"vs"$sharedSnps2".probs.bed
cat tmp/"$sharedSnps1"vs"$sharedSnps2".bed | bedtools intersect -v -a - -b tmp/"$sharedSnps1"vs"$sharedSnps2".probs.bed | bedtools sort -i - | uniq 
rm tmp/"$sharedSnps1"vs"$sharedSnps2".bed tmp/"$sharedSnps1"vs"$sharedSnps2".probs.bed








