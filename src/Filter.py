import argparse
import csv
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import numpy as bp
import scipy 
import pysam
import multiprocessing
import gzip

class RecessiveModel:
	def __init__(self, af = 1e-2, SBPV_cutoff=1e-3, DP_cutoff=7, AB_cutoff1=0.2, AB_cutoff2=0.8):
		self.SBPV_cutoff = SBPV_cutoff
		self.DP_cutoff = DP_cutoff
		self.AB_cutoff1 = AB_cutoff1
		self.AB_cutoff2 = AB_cutoff2
		self.LGD = set(["splice_acceptor_variant", "splice_donor_variant", "stop_gained", 
			"stop_lost", "start_lost", "frameshift_variant"])
		self.CADD_cutoff = 25
		self.REVEL_cutoff = 0.5
		self.AF_cutoff = af
		#self.C = ["syn","lgd","dmis","lgd/dmis"]
		self.C = ["syn","lgd","mis","cadd15","cadd20","cadd25","revel.5"]
		return

	# Return genotype as [A1, A2] if pass QC , return False otherwise
	def getINFO(self, info_string):
		infolist = info_string.split(';')
		infodict = {}
		for kv in infolist:
			kv = kv.split('=')
			if len(kv) == 2:
				k, v = kv
				infodict[k] = v
		return infodict

	def match_allele_csq(self, Ref, Alts, csq_head, csq_string):
		# Trim Leading Base
		Alts = Alts.split(",")
		if len(list(set([x[0] for x in Alts])))==1 and Ref[0] == list(set([x[0] for x in Alts]))[0]:
			_Ref = Ref[1:] if len(Ref[1:]) >0 else "-"
			_Alts = [Alt[1:] if len(Alt[1:]) >0 else "-" for Alt in Alts]
			#print (Ref, Alts, ";", _Ref, _Alts)
		else:
			_Alts = Alts
		res = {}
		csqs = csq_string.split(",")
		csqs = [dict(zip(csq_head, vep.split("|"))) for vep in csqs]
		for i, Alt in enumerate(Alts):
			res[Alt] = []
			for j, csq in enumerate(csqs):
				if csq["Allele"] == _Alts[i]:
					csq["Consequence"] = csq["Consequence"] .split("&")
					res[Alt].append(csq)
		return res

	def Recessive(self, Chr, GenotypeFil, VEPFil, AFFil, GenecodeFil):
		GenotypeFil = pysam.TabixFile(GenotypeFil)
		VEPFil = pysam.TabixFile(VEPFil)
		AFFil = pysam.TabixFile(AFFil)
		Genes, Transtripts = LoadGeneCode(GenecodeFil)
		CSQ_header = [X.strip().split("Format: ")[1].rstrip('>\"').split("|") for X in VEPFil.header if X.startswith("##INFO=<ID=CSQ")][0]
		Samples = GenotypeFil.header[-1].split("\t")[9:]
		OutFil1 = open("Rec.FamTest.Chr{}.txt".format(Chr), 'w')
		OutFil1.write("#Gene\t" + "\t".join(["{}.obs\t{}.haps".format(t,t) for t in self.C]) + "\n")
		OutFil2 = open("Rec.FamTest.Chr{}.sup.txt".format(Chr), 'w')
		OutFil2.write("#Gene\t" + "\t".join(["{}.NCantPhase\t{}.CantPhaseFams\t{}.N3vars".format(t,t,t) for t in self.C]) + "\n")
		OutFil3 = open("Rec.FamTest.Chr{}.sup.txt".format(Chr), 'w')
		OutFil3.write("#Gene\t" + "\t".join(["{}.NCantPhase\t{}.CantPhaseFams\t{}.N3vars".format(t,t,t) for t in self.C]) + "\n")
		#OutFil2.write("{}\t{}\t{}\n".format(Gene, "\t".join("{}\t{}\t{}".format(CantPhase[t], CantPhase_fams[t], MoreThanThree[t]) for t in self.C)))

		for i, (Gene,GTF) in enumerate(Genes.items()): # iterate through genes
			Gene_Fam_dat = {} # store genotypes for each fam, group by variant categories
			for cat in self.C:
				Gene_Fam_dat[cat] = {}
				for i, trio in enumerate(Trios):
					Gene_Fam_dat[cat][trio.FamID] = []
			start, end = int(GTF.start), int(GTF.end)
			veps, cohort, genotypes = [],[], []
			for term in VEPFil.fetch(Chr, start, end):
				veps.append(term)
			for term in AFFil.fetch(Chr, start, end):
				cohort.append(term)
			for term in GenotypeFil.fetch(Chr, start, end):
				genotypes.append(term)
			for var in zip(veps, cohort, genotypes):
				llist = var[0].split("\t")
				llist2 = var[2].split("\t")
				Chr, Pos, Ref, Alts = llist[0], llist[1], llist[3], llist[4]
				cohort_af = list(map(float, self.getINFO(var[1].split("\t")[7])["EUR_AF"].split(",")))
				fmt = llist2[8]
				Sample_genotypes = llist2[9:]
				infodict = self.getINFO(llist[7])
				Allele_CSQ_dict = self.match_allele_csq(Ref, Alts, CSQ_header, infodict["CSQ"])
				for i, Alt in enumerate(Alts.split(",")):
					var_k = "{}:{}:{}:{}".format(Chr, Pos, Ref, Alt)
					try:
						Allele_CSQ_dict[Alt][0]["gnomADg_AF_NFE"] = Allele_CSQ_dict[Alt][0]["gnomADg_AF_NFE"].split("&")[0]
						Allele_CSQ_dict[Alt][0]["gnomADe_AF_NFE"] = Allele_CSQ_dict[Alt][0]["gnomADe_AF_NFE"].split("&")[0]
						vep = Allele_CSQ_dict[Alt][0]
						gnomADg_af = 0 if (vep["gnomADg_AF_NFE"] == "" or vep["gnomADg_AF_NFE"] == ".")\
								else float(vep["gnomADg_AF_NFE"])
						gnomADe_af = 0 if (vep["gnomADe_AF_NFE"] == "" or vep["gnomADe_AF_NFE"] == ".")\
								else float(vep["gnomADe_AF_NFE"])
						af = cohort_af[i]
						#if gnomADg_af > 1e-2 or af > 1e-2 or af == 0:
						#if max(gnomADg_af, af) > 1e-2 or af == 0:
						if max(gnomADg_af, af) > self.AF_cutoff or af == 0:
							continue
						cons = Allele_CSQ_dict[Alt][0]["Consequence"]
						#print (cons)
						if len(set(cons).intersection(self.LGD))>= 1:
							#print (gnomADe_af, af, cons)
							#self.LookUpBiallic(i, "lgd", fmt, Sample_genotypes, Trios)
							Gene_Fam_dat = self.AddVar(i, var_k, "lgd", fmt, Sample_genotypes, Trios, Gene_Fam_dat)
						if "synonymous_variant" in set(cons):
							Gene_Fam_dat = self.AddVar(i, var_k, "syn", fmt, Sample_genotypes, Trios, Gene_Fam_dat)
						if "missense_variant" in set(cons):
							Gene_Fam_dat = self.AddVar(i, var_k, "mis", fmt, Sample_genotypes, Trios, Gene_Fam_dat)
						if "missense_variant" in set(cons) and float(Allele_CSQ_dict[Alt][0]["CADD_PHRED"]) > 15:
							Gene_Fam_dat = self.AddVar(i, var_k, "cadd15", fmt, Sample_genotypes, Trios, Gene_Fam_dat)
						if "missense_variant" in set(cons) and float(Allele_CSQ_dict[Alt][0]["CADD_PHRED"]) > 20:
							Gene_Fam_dat = self.AddVar(i, var_k, "cadd20", fmt, Sample_genotypes, Trios, Gene_Fam_dat)
						if "missense_variant" in set(cons) and float(Allele_CSQ_dict[Alt][0]["CADD_PHRED"]) > 25:
							Gene_Fam_dat = self.AddVar(i, var_k, "cadd25", fmt, Sample_genotypes, Trios, Gene_Fam_dat)
						try:
							if "missense_variant" in set(cons) and float(Allele_CSQ_dict[Alt][0]["REVEL_score"]) > 0.5:
								Gene_Fam_dat = self.AddVar(i, var_k, "revel.5", fmt, Sample_genotypes, Trios, Gene_Fam_dat)
						except ValueError:
							continue
					except KeyError as e:
						print(e)
						print("KeyError", Ref, Alts, Alt, Allele_CSQ_dict)
						return
					except IndexError:
						print("IndexError", Ref, Alts, llist[7], Allele_CSQ_dict)
						return
			res = self.Phasing_N_Count(Gene_Fam_dat, Trios)
			OBS = {}
			EXP = {}
			CantPhase = {}
			MoreThanThree = {}
			CantPhase_fams = {}
			for t in self.C:
				OBS[t] = (res[t][0] + res[t][1])
				EXP[t] = (res[t][2])
				CantPhase[t] = res[t][3]
				CantPhase_fams[t] = ",".join(res[t][4])
				MoreThanThree[t] = res[t][5]
			OutFil.write("{}\t{}\n".format(Gene, "\t".join("{}\t{}".format(OBS[t], EXP[t]) for t in self.C)))
			OutFil2.write("{}\t{}\n".format(Gene, "\t".join("{}\t{}\t{}".format(CantPhase[t], CantPhase_fams[t], MoreThanThree[t]) for t in self.C)))
		return
	

def GetOptions():
	parser = argparse.ArgumentParser()
	#parser.add_argument("-v", "--vcf", type=str, required=True, help="<Required> VCF file")
	#parser.add_argument("-l", "--list", type=str, required=True, help="<Required> Indv List file")
	#parser.add_argument("-o", "--out", type=str, help="<Required> Output VCF file")
	parser.add_argument("--chr", type=str, help="<Required> Output VCF file")
	parser.add_argument("--af", type=float, default=1e-2, help="<Required> Output VCF file")
	args = parser.parse_args()
	#if args.out == None:
	#	args.out = "test.out.vcf"
	return args

def main():
	args = GetOptions()
	ins = RecessiveModel(args.af)
	#List = [l.strip() for l in open(args.list)]
	#ins.ComputeSiteAF(args.vcf, List, "EUR", args.out)
	Chr = args.chr
	GenotypeFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K/VCF/TrioVCF/Genotypes/SPARK30K.TrioSamples.Chr{}.vcf.gz".format(Chr)
	VEPFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K/VCF/TrioVCF/sites/SPARK30K.TrioSamples.Chr{}.vep.vcf.gz".format(Chr)
	#VEPFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K/VCF/TrioVCF/sites/Annotated2/SPARK30K.TrioSamples.Chr{}.vep.vcf.gz".format(Chr)
	AFFil = "/home/local/users/jw/Genetics_Projects/SPARK/spark_genomics/dat/SPARK30K.TrioSamples.Chr{}.eurAF.vcf.gz".format(Chr)
	genecode = "/home/local/users/jw/vep_data/homo_sapiens/GeneCodeV29/CHRs/genecodev29.{}.gtf".format(Chr)
	ins.Recessive(Chr, GenotypeFil ,VEPFil, AFFil, genecode)

if __name__=='__main__':
	main()
