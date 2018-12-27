cat $1 | awk '{print$1,$2,$2+1,"NA",$3}' | tr ' ' '\t' | bedtools sort -i - > $2
