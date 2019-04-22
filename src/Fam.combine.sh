head -n1 Rec.FamTest.Chr1.sup.txt > head.sup
cat Rec.FamTest.Chr*.sup.txt |grep -v "^#"|cat head.sup - > Rec.Sup.FamTest.tsv
mkdir -p sup
mv Rec.FamTest.Chr*.sup.txt sup
head -n1 Rec.FamTest.Chr1.txt > head
cat Rec.FamTest.Chr*.txt |grep -v "^#"|cat head - > Rec.FamTest.tsv
