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
#from vcftidy import * 

class RecessiveModel:
	def __init__(self, af = 1e-2, SBPV_cutoff=1e-5, DP_cutoff=8, AB_cutoff1=0.25, AB_cutoff2=0.8):
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
		self.C = ["syn","lgd","mis","cadd15","cadd20","cadd25","revel.5","mvp2.85","mpc1","lgd_cadd25"]
		return

	# compute NHet, NHom, AC and AF of each site, within a population, Write to a VCF file.
	def ComputeSiteAF(self, InpVCF, Indvs, prefix, OutVCF):
		fin = gzip.open(InpVCF, 'rt')
		fout = open(OutVCF, 'w')
		for l in fin:
			if l.startswith("##"):
				fout.write(l)
				continue
			elif l.startswith("#"):
				llist = l.strip().split("\t")
				#Head = llist[]
				indvs = llist[9:]
				index = self.GetIndex(Indvs, indvs)
				fout.write("##ComputeAF={}\n".format(prefix))
				fout.write("\t".join(llist[:9])+"\n")
			else:
				llist = l.strip().split("\t")
				length = (len(llist[4].split(","))+1)
				AFs = [0] * length
				ACs = [0] * length
				AC_Het = [0] * length
				AC_Hom = [0] * length
				AN = 0 
				GTs = llist[9:]
				for idx in index:
					GT = self.GenotypeQC(llist[8], GTs[idx])
					if GT:
						A1, A2 = GT[0], GT[1]
						ACs[A1] += 1
						ACs[A2] += 1
						if A1 != A2:
							AC_Het[A2] += 1
						else:
							AC_Hom[A2] += 1
						AN += 2
				for i in range(length):
					try:
						AFs[i] = str(float(ACs[i])/AN)
					except ZeroDivisionError:
						AFs[i] = "0"
					ACs[i] = str(ACs[i])
					AC_Het[i] = str(AC_Het[i])
					AC_Hom[i] = str(AC_Hom[i])
				New = ";{}_AN={};{}_AC={};{}_AF={};{}_AC_Het={};{}_AC_Hom={}".format(prefix, AN, prefix, ",".join(ACs[1:]), prefix, ",".join(AFs[1:]), prefix, ",".join(AC_Het[1:]), prefix, ",".join(AC_Hom[1:]))
				llist[7] = llist[7] + New
				fout.write("\t".join(llist[:8])+"\n") 
		return
	def GetIndex(self, Indvs, indvs):
		res = []
		for indv in Indvs:
			try:
				res.append(indvs.index(indv))
			except:
				continue
		return res

	# Return genotype as [A1, A2] if pass QC , return False otherwise
	def GenotypeQC(self, fmt, gt_dat):
		tmp = {}
		for k,v in zip(fmt.split(":"), gt_dat.split(":")):
			tmp[k] = v
		GT = tmp["GT"].split("/")
		if tmp["GQ"] == ".":
			#return False
			return [0,0]
		elif float(tmp["GQ"]) < 60:
			return False
		if GT[0] == "." or GT[1] == ".":
			return False
		if tmp["GT"] != "0/0":
			#print(fmt, gt_dat)
			#if "." in tmp["SBPV"]:
			#	pass
			#elif float(tmp["SBPV"]) < self.SBPV_cutoff:
			#	return False
			if float(tmp["DP"]) < self.DP_cutoff:
				return False
			if float(tmp["AD"].split(",")[1])/float(tmp["DP"]) < self.AB_cutoff1: # or float(tmp["AD"].split(",")[1])/float(tmp["DP"]) > self.AB_cutoff2:
				return False
		return [int(i) for i in GT]

	def LoadPedigree(self, PedFil, Samples):
		PedFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K_07/recessive/EUR_Fams.ped"
		Fams = []
		reader = csv.reader(open(PedFil, 'rt'), delimiter="\t")
		PreFamID, tmp = None, None
		for row in reader:
			row.append(Samples.index(row[1])) # Add sample index in VCF header, to locate genotype
			FamID = row[0]
			if FamID != PreFamID:
				if tmp != None:
					Fams.append(tmp)
				tmp = Family(FamID)
				PreFamID = FamID
				tmp.Proband = Sample(row)
			else:
				if row[1] == tmp.Proband.Father:
					tmp.Father = Sample(row)
				elif row[1] == tmp.Proband.Mother:
					tmp.Mother = Sample(row)
				else:
					tmp.Siblings.append(Sample(row))
		return Fams

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

	def search_severe_consequence(self, var_k, Allele_CSQ_dict, Alt):
		severe_consequence = None
		severe_trans = None
		for i in range(len(Allele_CSQ_dict[Alt])):
			consequence = Allele_CSQ_dict[Alt][i]["Consequence"]
			Transcript = Allele_CSQ_dict[Alt][i]["Feature"]
			#if var_k == "11:7994466:T:G":
			#   print(consequence, Transcript)
			#if set(consequence).intersection(self.LGD) >= 1:
			if len(set(consequence).intersection(self.LGD))>= 1:
				return i, consequence, Transcript
			elif consequence[0] == "missense_variant":
				return i, consequence, Transcript
			elif consequence[0] == "synonymous_variant":
				#print(consequence)
				severe_consequence = consequence
				severe_trans = Transcript
		if severe_consequence == None:
			return 0, "non-coding", None
		else:
			return 0, severe_consequence, severe_trans

	def Recessive(self, Chr, GenotypeFil, VEPFil, AFFil, GenecodeFil):
		#GenotypeFil = pysam.TabixFile(GenotypeFil)
		#VEPFil = pysam.TabixFile(VEPFil)
		#AFFil = pysam.TabixFile(AFFil)
		#Genes, Transtripts = LoadGeneCode(GenecodeFil)
		#CSQ_header = [X.strip().split("Format: ")[1].rstrip('>\"').split("|") for X in VEPFil.header if X.startswith("##INFO=<ID=CSQ")][0]
		GenotypeFil = gzip.open(GenotypeFil, 'rb')
		VEPFil = gzip.open(VEPFil, 'rb')
		DIR = "/home/local/users/jw/Genetics_Projects/SPARK/30K_07/VCF/Filt/"
		OutFil1 = open(DIR+"SPARK30K.Anno.Filt.Chr{}.vcf".format(Chr), 'wb')
		OutFil2 = open(DIR+"SPARK30K.Genotypes.Filt.Chr{}.vcf".format(Chr), 'wb')
		# Skip Header Genotype
		while 1:
			l = GenotypeFil.readline()#.decode("utf-8")
			OutFil2.write(l)
			if l.decode("utf-8").startswith("#C"):
				break
		# Skip Header VEP 
		while 1:
			l = VEPFil.readline()
			OutFil1.write(l)
			if l.decode("utf-8").startswith("##INFO=<ID=CSQ"):
				CSQ_header = l.decode("utf-8").split("Format: ")[1].rstrip('>\"').split("|")
			if l.decode("utf-8").startswith("#C"):
				break

		while 1: # iterate through genes
			VEPLine = VEPFil.readline()
			GenotypesLine = GenotypeFil.readline()
			if not VEPLine:
				break
			llist = VEPLine.decode("utf-8").split("\t")
			#llist2 = GenotypesLine.decode("utf-8").split("\t")
			Chr, Pos, Ref, Alts = llist[0], llist[1], llist[3], llist[4]
			infodict = self.getINFO(llist[7])
			Allele_CSQ_dict = self.match_allele_csq(Ref, Alts, CSQ_header, infodict["CSQ"])
			KEEP = False
			for i, Alt in enumerate(Alts.split(",")):
				var_k = "{}:{}:{}:{}".format(Chr, Pos, Ref, Alt)
				Allele_CSQ_dict[Alt][0]["gnomADg_AF_NFE"] = Allele_CSQ_dict[Alt][0]["gnomADg_AF_NFE"].split("&")[0]
				Allele_CSQ_dict[Alt][0]["gnomADe_AF_NFE"] = Allele_CSQ_dict[Alt][0]["gnomADe_AF_NFE"].split("&")[0]
				vep = Allele_CSQ_dict[Alt][0]
				gnomADg_af = 0 if (vep["gnomADg_AF_NFE"] == "" or vep["gnomADg_AF_NFE"] == ".")\
						else float(vep["gnomADg_AF_NFE"])
				gnomADe_af = 0 if (vep["gnomADe_AF_NFE"] == "" or vep["gnomADe_AF_NFE"] == ".")\
						else float(vep["gnomADe_AF_NFE"])
				if gnomADg_af > self.AF_cutoff:
					continue
				idx_anno, cons, trans = self.search_severe_consequence(var_k, Allele_CSQ_dict, Alt)
				#if cons[0] == "synonymous_variant":
				#	print(cons)
				if len(set(cons).intersection(self.LGD))>= 1 or "synonymous_variant" in set(cons) or "missense_variant" in set(cons):
					KEEP = True
			if KEEP:
				OutFil1.write(VEPLine)
				OutFil2.write(GenotypesLine)
		return

	
def GetOptions():
	parser = argparse.ArgumentParser()
	#parser.add_argument("-v", "--vcf", type=str, required=True, help="<Required> VCF file")
	#parser.add_argument("-l", "--list", type=str, required=True, help="<Required> Indv List file")
	#parser.add_argument("-o", "--out", type=str, help="<Required> Output VCF file")
	parser.add_argument("--chr", type=str, help="<Required> Output VCF file")
	parser.add_argument("--af", type=float, help="<Required> Output VCF file")
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
	GenotypeFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K_07/VCF/GenotypesSplitbyChr/GATK4_20190729.chr{}.vcf.gz".format(Chr)
	VEPFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K_07/VCF/SitesSplitbyChr/annotated/GATK4_20190729.chr{}.mappability.vcf.gz".format(Chr)
	AFFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K_07/recessive/AF/GATK4_20190729.chr{}.eurAF.vcf.gz".format(Chr)
	genecode = "/home/local/users/jw/vep_data/homo_sapiens/GeneCodeV29/CHRs/genecodev29.{}.gtf".format(Chr)
	ins.Recessive(Chr, GenotypeFil ,VEPFil, AFFil, genecode)

if __name__=='__main__':
	main()
