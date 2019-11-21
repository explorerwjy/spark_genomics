#mkdir -p AF1e-2; cd AF1e-2;
nohup parallel -j 15 python FilterSites.py --chr {} --af 1e-2 ::: $(seq 22) &
