#for chr in $(seq 22);do
	#echo $chr
	#nohup python Recessive.py --chr $chr > Chr.${chr}.txt &
	#nohup python Recessive.py --chr $chr &
#done
mkdir -p AF1e-2; cd AF1e-2;
nohup parallel -j 10 python ../Recessive.py --chr {} --af 1e-2 ::: $(seq 22) &
cd ../
mkdir -p AF1e-3; cd AF1e-3;
nohup parallel -j 10 python ../Recessive.py --chr {} --af 1e-3 ::: $(seq 22) &
