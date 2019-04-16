head -n1 Rec.Chr1.txt > head
cat Rec.Chr*.txt |grep -v "^#"|cat head - > Rec.all_chr.tsv
