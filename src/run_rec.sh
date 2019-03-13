for chr in $(seq 22);do
	#echo $chr
	nohup python Recessive.py --chr $chr > Chr.${chr}.txt &
done
