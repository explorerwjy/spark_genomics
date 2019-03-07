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
	def __init__(self, SBPV_cutoff=1e-3, DP_cutoff=7, AB_cutoff=0.1):
		self.SBPV_cutoff = SBPV_cutoff
		self.DP_cutoff = DP_cutoff
		self.AB_cutoff = AB_cutoff
		self.LGD = set(["splice_acceptor_variant", "splice_donor_variant", "stop_gained", 
			       "stop_lost", "start_lost", "frameshift_variant"])
		self.CADD_cutoff = 25
		self.AF_cutoff = 1e-2
		self.C = ["syn","lgd","dmis","lgd/dmis","dmis,lgd"]
		return

	# compute NHet, NHom, AC and AF of each site, within a population, Write to a VCF file.
	def ComputeSiteAF(self, InpVCF, Indvs, prefix, OutVCF):
		fin = gzip.open(InpVCF, 'rt')
		fout = open(OutVCF, 'wt')
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
			return False
		elif float(tmp["GQ"]) < 30:
			return False
		if GT[0] == "." or GT[1] == ".":
			return False
		if tmp["GT"] != "0/0":
			#print(fmt, gt_dat)
			if "." in tmp["SBPV"]:
				pass
			elif float(tmp["SBPV"]) < self.SBPV_cutoff:
				return False
			if float(tmp["DP"]) < self.DP_cutoff:
				return False
			if float(tmp["AD"].split(",")[1])/float(tmp["DP"]) < self.AB_cutoff:
				return False
		return [int(i) for i in GT]

	def LoadPedigree(self, PedFil):
		PedFil = "/home/local/users/jw/Genetics_Projects/SPARK/spark_genomics/dat/EUR_Trios.ped"
		Trios = []
		reader = csv.reader(open(PedFil, 'rt'), delimiter="\t")
		counter = 0
		for row in reader:
			counter += 1
			if counter == 1:
				tmp = Family(row[0])
				tmp.Proband = Sample(row)
			if counter == 2:
				tmp.Father = Sample(row)
			if counter == 3:
				tmp.Mother = Sample(row)
				counter = 0
		return Trios
	def getINFO(info_string):
		infolist = info_string.split(';')
		infodict = {}
		for kv in infolist:
			kv = kv.split('=')
			if len(kv) == 2:
				k, v = kv
				infodict[k] = v
		return infodict

	def match_allele_csq(Ref, Alts, csq_head, csq_string):
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

	def Recessive(self, GenotypeFil, VEPFil, AFFil):
		GenotypeFil = pysam.TabixFile("/home/local/users/jw/Genetics_Projects/SPARK/30K/VCF/TrioVCF/Genotypes/SPARK30K.TrioSamples.Chr21.vcf.gz")
		VEPFil = pysam.TabixFile("/home/local/users/jw/Genetics_Projects/SPARK/30K/VCF/TrioVCF/sites/SPARK30K.TrioSamples.Chr21.vep.vcf.gz")
		AFFil = pysam.TabixFile("/home/local/users/jw/Genetics_Projects/SPARK/spark_genomics/dat/SPARK30K.TrioSamples.Chr21.eurAF.vcf.gz")
		genecodechr21 = "/home/local/users/jw/vep_data/homo_sapiens/GeneCodeV29/CHRs/genecodev29.21.gtf"
		Genes, Transtripts = LoadGeneCode(genecodechr21)
		CSQ_header = [X.strip().split("Format: ")[1].rstrip('>\"').split("|") for X in VEPFil.header if X.startswith("##INFO=<ID=CSQ")][0]
		#print(CSQ_header)
		Samples = GenotypeFil.header[-1].split("\t")[9:]
		infodict = getINFO(llist[7])
		Allele_CSQ_dict = match_allele_csq(Ref, Alts, CSQ_header, infodict["CSQ"])
		"""
		for trio in Trios:
			Nhaps, haps = 0,0
			for var in zip(veps, cohort, genotypes):
				llist = var[0].split("\t")
				Chr, Pos, Ref, Alts = llist[0], llist[1], llist[3], llist[4]
				Sample_genotypes = genotypes[9:]
				infodict = getINFO(llist[7])
				Allele_CSQ_dict = match_allele_csq(Ref, Alts, CSQ_header, infodict["CSQ"])
		"""
		for i,trio in enumerate(Trios):
			Trios[i].Haps = {}
			for cat in ["syn","lgd","dmis","lgd/dmis","dmis,lgd"]
				Trios[i].pro_haps[cat] = [0,0]
				Trios[i].fa_haps[cat] = [0,0]
				Trios[i].mo_haps[cat] = [0,0]
		for var in zip(veps, cohort, genotypes):
			llist = var[0].split("\t")
			Chr, Pos, Ref, Alts = llist[0], llist[1], llist[3], llist[4]
			cohort_af = list(map(float, getINFO(var[1].split("\t")[7])["EUR_AF"].split(",")))
			Sample_genotypes = genotypes[9:]
			infodict = getINFO(llist[7])
			Allele_CSQ_dict = match_allele_csq(Ref, Alts, CSQ_header, infodict["CSQ"])
			for i, Alt in enumerate(Alts.split(",")):
				#print(Alt, Allele_CSQ_dict[Alt][0])
				try:
					Allele_CSQ_dict[Alt][0]["gnomADg_AF_NFE"] = Allele_CSQ_dict[Alt][0]["gnomADg_AF_NFE"].split("&")[0]
					gnomADg_af = 0 if (vep["gnomADg_AF_NFE"] == "" or vep["gnomADg_AF_NFE"] == ".")\
							else float(vep["gnomADg_AF_NFE"])
					Allele_CSQ_dict[Alt][0]["gnomADe_AF_NFE"] = Allele_CSQ_dict[Alt][0]["gnomADe_AF_NFE"].split("&")[0]
					gnomADe_af = 0 if (vep["gnomADe_AF_NFE"] == "" or vep["gnomADe_AF_NFE"] == ".")\
							else float(vep["gnomADe_AF_NFE"])
					af = cohort_af[i]
					if gnomADg_af > 1e-2 or af > 1e-2 or af == 0:
						continue
					cons = Allele_CSQ_dict[Alt][0]["Consequence"]
					#print (cons)
					if len(set(cons).intersection(LGD))>= 1:
						print (gnomADe_af, af, cons)
						LookUpBiallic(i, "lgd", Sample_genotypes, Trops)
					elif "synonymous_variant" in set(cons):
						print (gnomADe_af, af, cons)
						LookUpBiallic(i, "syn", Sample_genotypes, Trops)
					elif "missense_variant" in set(cons) and float(Allele_CSQ_dict[Alt][0]["CADD_PHRED"]) > 25:
						print (gnomADe_af, af, cons, float(Allele_CSQ_dict[Alt][0]["CADD_PHRED"]))
						LookUpBiallic(i, "dmis", Sample_genotypes, Trops)
				except KeyError:
					print("KeyError", Alt, Allele_CSQ_dict)
				except IndexError:
					print("IndexError", Ref, Alts, llist[7], Allele_CSQ_dict)
		Obs = {}
		for t in self.C:
			OBs[t] = 0
		for i,Trio in enumerate(Trios):
			for t in self.C
				if Trios[i].pro_haps[t] == [1,1]:
					Obs[t] += 1
		print("{}\t{}".format(gene, "\t".join(Obs[t] or t in self.C)))
						
		return
	def LookUpBiallic(i, Vartype, var, gts, Trios):
		for i, (prob, fa, mo) in enumerate(trios):
			gt_prob, gt_fa, gt_mo = r.GenotypeQC(gts[prob[6]]), r.GenotypeQC(gts[fa[6]]), r.GenotypeQC(gts[mo[6]])
			if gt_prob == False or gt_fa == False or gt_mo == False:
				continue # gt failed QC
			elif (gt_prob[0] not in gt_fa and gt[1] not in gt_mo) or (gt_prob[1] not in gt_fa and gt[0] not in gt_mo):
				continue # mendelian error
			else:
				if gt_fa[0] == i:
					Trios[i].fa_haps[Vartpye][0] == 1
				if gt_fa[1] != i:
					Trios[i].fa_haps[Vartype][1] == 1
				if gt_mo[0] != i:
					Trios[i].mo_haps[Vartype][0] == 1
				if gt_mo[1] != i:
					Trios[i].mo_haps[Vartype][1] == 1
				if gt_prob[1] == gt_fa[1]: #paternal var
					Trios[i].pro_haps[Vartype][0] == 1
				if gt_prob[1] == gt_mo[1]: #maternal var
					Trios[i].pro_haps[Vartype][1] == 1
                
	# compute Aggregated Allele freq for each gene
	def ComputeCumulativeAF():
		return 

class Sample:
    def __init__(self, row):
        self.FamID = row["FamID"]
        self.sampleID = row["SampleID"]
        self.Father = row["Paternal"]
        self.Mother = row["Maternal"]
        self.Sex = row["Sex"]
        self.Affected = row["Affected"]
    def show(self):
        print self.FamID, self.sampleID, self.Father, self.Mother, self.Sex, self.Affected
    def display(self):
        #return "\t".join([self.FamID, self.sampleID, self.Father, self.Mother, self.Sex, self.Affected])
        return list(map(str, [self.FamID, self.sampleID, self.Father, self.Mother, self.Sex, self.Affected]))

class Family:
    def __init__(self, FamID):
        self.FamID = FamID
        self.Father = None
        self.Mother = None
        self.Proband = None
        self.Siblings = []
    def show(self):
        print "FamID:{} Proband:{} Father:{} Mother:{} Siblings:{}".format(
        self.FamID, self.Proband.sampleID, self.Father.sampleID, self.Mother.sampleID, ", ".join(
            [x.sampleID for x in self.Siblings]))

class GTFRecord:
    def __init__(self, Chr, source, Type, start, end, strand, info):
        self.Chr = Chr
        self.source = source
        self.Type = Type
        self.start = start
        self.end = end
        self.strand = strand
        self.info = info

def gtf_info_parser(info):
    res = {}
    for term in info.split(";"):
        if term == "":
            continue
        #print(">",term)
        key,v = term.split()
        v = v.strip('"')
        res[key]=v
    return res
        
def LoadGeneCode(genecodefil):
    Genes = {}
    Transcripts = {}
    hand = open(genecodefil, 'rt')
    for l in hand:
        if l.startswith("#"):
            continue
        llist = l.strip().split("\t")
        info = gtf_info_parser(llist[8])
        if llist[2] == "gene":
            Genes[info["gene_name"]] = GTFRecord(llist[0], llist[1], llist[2], llist[3], llist[4], llist[6], info)
            Transcripts[info["gene_name"]] = []
        elif llist[2] == "transcript":
            Transcripts[info["gene_name"]].append(GTFRecord(llist[0], llist[1], llist[2], llist[3], llist[4], llist[6], info))
    return Genes, Transcripts 
	# 
def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument("-v", "--vcf", type=str, required=True, help="<Required> VCF file")
	parser.add_argument("-l", "--list", type=str, required=True, help="<Required> Indv List file")
	parser.add_argument("-o", "--out", type=str, help="<Required> Output VCF file")
	args = parser.parse_args()
	if args.out == None:
		args.out = "test.out.vcf"
	return args

def main():
	args = GetOptions()
	ins = RecessiveModel()
	List = [l.strip() for l in open(args.list)]
	ins.ComputeSiteAF(args.vcf, List, "EUR", args.out)

if __name__=='__main__':
	main()
