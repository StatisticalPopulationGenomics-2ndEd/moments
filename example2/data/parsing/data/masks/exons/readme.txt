done with 
zcat exome_hg38.tsv.gz | cut -f3,4,5 | tail -n +2 | awk '$0="chr"$0' | bedtools merge | awk '{print $0 >> "exons_"$1".bed"}'
