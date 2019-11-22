
function extrac_site () {
	name=$(basename $1|sed  s/.vcf.gz/.sites.vcf.gz/g)
	#nohup less $1 |cut -f 1-8 |bgzip > $name &
	nohup bcftools view -G -O v $1 |bgzip > $name &
}

#for file in `less /home/local/users/jw/Genetics_Projects/SPARK/30K/VCF/SPARK30K.TrioSamples.vcf.list`;
for file in `less ${1}`;
do
	extrac_site $file
done

