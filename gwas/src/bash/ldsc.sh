#Choose which set of gene sets to analyze. Options include Multi_tissue_gene_expr, Multi_tissue_chromatin, GTEx_brain, Cahoy, ImmGen, or Corces_ATAC
#cts_name=Cahoy 
#cts_name=GTEx_brain 
cts_name=Multi_tissue_gene_expr 

#Download the LD scores
#wget -c https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/${cts_name}_1000Gv3_ldscores.tgz
#tar -xvzf ${cts_name}_1000Gv3_ldscores.tgz

#Download and format the summary statistics, or use your own.
#python munge_sumstats.py \
#	--sumstats body_BMIz.sumstats.gz \
#	--merge-alleles w_hm3.snplist \
#	--out UKBB_BMI

#Run the regression
ldsc.py \
	--h2-cts spark.sumstats.gz \
	--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
	--out SPARK_${cts_name} \
	--ref-ld-chr-cts $cts_name.ldcts \
	--w-ld-chr weights_hm3_no_hla/weights.
