#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# Filter.py
# Filter Family VCF and ready for recessive genotype testing. 
# Input:
#	Fam.VCF; Fam.Ped; 
#	Transcripts; Annotations;
# Output:
#	QC and filtered Fam.VCF
#========================================================================================================

import argparse
import gz
import pysam
import csv
from Recessive import *

class Filter:
	def __init__(self, args):
		self.vcf = pysam.TabixFile(args.vcf)
		self.ped = LoadPedigree(args.pedigree)
		self.ann = pysam.TabixFile(args.annotation)
		self.gtf = ReadGeneGTF(args.gtf)
		self.out = open(args.out, 'wt')

	def run2(self):
		for gene in self.gtf:
			for Chr, Start, End in gene.intervals:
				variants = self.vcf.fetch(Chr, Start, End) 
				annos = self.ann.fetch(Chr, Start, End)
				MakeAnnDict(annos)
				for variant in variants:
					record = record.split("\t")
					Chr, pos, ID, ref, Alts, Qual, Filter, INFO, Format = record[:9]
					Genotypes = record[9:]
						

def getINFO(self, info_string):
	infolist = info_string.split(';')
	infodict = {}
	for kv in infolist:
		kv = kv.split('=')
		if len(kv) == 2:
			k, v = kv
			infodict[k] = v
	return infodict

def TrimAllele(pos, ref, alt):
	while 1:
		if len(ref) == 0:
			ref = "-"
			return pos, ref, alt
		elif len(alt) == 0:
			alt = "-"
			return pos, ref, alt
		elif ref[0] == alt[0]: # left trim
			ref, alt = ref[1:], alt[1:]
			pos -= 1
		elif ref[-1] == alt[-1]: # right trim
			ref, alt = ref[:-1], alt[:-1]
		else:
			return pos, ref, alt


def match_allele_csq(res, Chr, Pos, Ref, Alts, csq_head, csq_string):
	Alts = Alts.split(",")
	csqs = csq_string.split(",")
	csqs = [dict(zip(csq_head, vep.split("|"))) for vep in csqs]
	for i, Alt in enumerate(Alts):
		res[":".join(Chr, pos, ref, alt)] = []
		pos, ref, alt = TrimAllele(Pos, Ref, Alt)
		for j, csq in enumerate(csqs):
			if csq["Allele"] == _Alts[i]:
				csq["Consequence"] = csq["Consequence"] .split("&")
				res[":".join(pos, ref, alt)].append(csq)
	return res
	

def MakeAnnDict(annos):
	res = {}
	for record in annos:
		record = record.split("\t")
		Chr, pos, ID, ref, Alts, Qual, Filter, INFO = record
		CSQ = getINFO(INFO)["CSQ"]
		res = match_allele_csq(res, Chr, pos, ref, Alts, csq_head, csq_string)
	return res

def MakeOutName(args.vcf):
	names = args.vcf.split("/")[-1].split(".")
	if names[-1] == "gz":
		names = names[:-1]
	names[-1] = "qc"
	names.append("vcf")
	return ".".join(names) 

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('-p','--pedigree', type=str, help = 'family pedigree file')
	parser.add_argument('-v','--vcf', type=str, help = 'family vcf file')
	parser.add_argument('-a','--annotation', type=str, help = 'annotation vcf file')
	parser.add_argument('--gtf', type=str, help = 'transcript gtf file')
	parser.add_argument('-o','--out', type=str, help = 'output vcf file')
	args = parser.parse_args()
	if args.out == None:
		args.out = MakeOutName(args.vcf)
	return args

def main():
	args = GetOptions()
	ins = Filter(args)
	ins.run()	

	return

if __name__=='__main__':
	main()
