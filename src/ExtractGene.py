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
	def __init__(self, af = 1e-2, SBPV_cutoff=1e-3, DP_cutoff=5, AB_cutoff1=0.1, AB_cutoff2=0.7):
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
		#self.C = ["syn","lgd","mis","cadd15","cadd20","cadd25","revel.5"]
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
		if GT == ['1', '1']:
			print(GT)
		if tmp["GQ"] == ".":
			#return False
			return [0, 0]
		elif float(tmp["GQ"]) < 60:
			return False
		if GT[0] == "." or GT[1] == ".":
			return False
		if tmp["GT"] != "0/0":
			print(tmp)
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
		#PedFil = "/home/local/users/jw/Genetics_Projects/SPARK/spark_genomics/dat/EUR_Fams.ped"
		PedFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K_07/recessive/EUR_Fams.ped"
		#PedFil = "SF0044997.ped"
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
		Fams.append(tmp)
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
			#	print(consequence, Transcript)
			#if set(consequence).intersection(self.LGD) >= 1:
			if len(set(consequence).intersection(self.LGD))>= 1:
				return i, consequence, Transcript
			elif consequence[0] == "missense_variant":
				return i, consequence, Transcript
			elif consequence[0] == "synonymous_variant":
				severe_consequence == consequence
				severe_trans = Transcript 
		if severe_consequence == None:
			return 0, "non-coding", None
		else:
			return 0, severe_consequence, severe_trans

	def Recessive(self, Chr, Gene, GenotypeFil, VEPFil, AFFil, GenecodeFil):
		Gene = Gene
		GenotypeFil = pysam.TabixFile(GenotypeFil)
		#GenotypeFil = pysam.TabixFile("test.vcf.gz")
		VEPFil = pysam.TabixFile(VEPFil)
		#VEPFil = pysam.TabixFile("test.vep.vcf.gz")
		AFFil = pysam.TabixFile(AFFil)
		Genes, Transtripts = LoadGeneCode(GenecodeFil)
		CSQ_header = [X.strip().split("Format: ")[1].rstrip('>\"').split("|") for X in VEPFil.header if X.startswith("##INFO=<ID=CSQ")][0]
		Samples = GenotypeFil.header[-1].split("\t")[9:]
		print(Samples[2631])
		OutFil = csv.writer(open("Rec.Chr{}.{}.variants.tsv".format(Chr, Gene), 'w'), delimiter="\t")
		Header = ["Chr", "Pos", "Ref", "Alt", "effect", "Consequence", "gnomADg_AF_NFE", "CADD_PHRED", "REVEL_score","MPC_score","MVP2_RankScore", "FamID", "SampleID", "Genotype", "Role", "N_kid_in_Fam"]
		OutFil.writerow(Header)
		Trios = self.LoadPedigree("a", Samples)

		GTF  = Genes[Gene]
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
					#print(llist2, gnomADg_af, af)
					#if max(gnomADg_af, af) > self.AF_cutoff or af == 0:
					if max(gnomADg_af, af) > self.AF_cutoff:
						continue
					#cons = Allele_CSQ_dict[Alt][0]["Consequence"]
					idx_anno, cons, trans = self.search_severe_consequence(var_k, Allele_CSQ_dict, Alt)
					if var_k == "11:7994466:T:G":
						print(cons, trans)
					#print(llist2, cons)
					if len(set(cons).intersection(self.LGD))>= 1:
						Gene_Fam_dat = self.AddVar(i, var_k, idx_anno, "lgd", fmt, Sample_genotypes, Trios, Gene_Fam_dat, Allele_CSQ_dict)
					if "synonymous_variant" in set(cons):
						Gene_Fam_dat = self.AddVar(i, var_k, idx_anno, "syn", fmt, Sample_genotypes, Trios, Gene_Fam_dat, Allele_CSQ_dict)
					if "missense_variant" in set(cons):
						Gene_Fam_dat = self.AddVar(i, var_k, idx_anno, "mis", fmt, Sample_genotypes, Trios, Gene_Fam_dat, Allele_CSQ_dict)
					if ("missense_variant" in set(cons) and float(Allele_CSQ_dict[Alt][0]["CADD_PHRED"]) > 25) or (len(set(cons).intersection(self.LGD))>= 1):
						Gene_Fam_dat = self.AddVar(i, var_k, idx_anno, "lgd_cadd25", fmt, Sample_genotypes, Trios, Gene_Fam_dat, Allele_CSQ_dict)
				except KeyError as e:
					print(e)
					print("KeyError", Ref, Alts, Alt, Allele_CSQ_dict)
					return
				except IndexError:
					print("IndexError", Ref, Alts, llist[7], Allele_CSQ_dict)
					return
		res = self.Phasing_N_Count(Gene_Fam_dat, Trios, OutFil)
		return

	def AddVar(self, i, var_k, idx_anno, Vartype, fmt, gts, Trios, Gene_Fam_dat, Allele_CSQ_dict):
		print(var_k, Vartype)
		var_info = []
		var_coord = var_k.split(":")
		var_info.extend(var_coord)
		alt = var_coord[3]
		var_info.append(Vartype)
		var_info.append(Allele_CSQ_dict[alt][idx_anno]["Consequence"][0])
		var_info.append(Allele_CSQ_dict[alt][idx_anno]["gnomADg_AF_NFE"])
		var_info.append(Allele_CSQ_dict[alt][idx_anno]["CADD_PHRED"])
		var_info.append(Allele_CSQ_dict[alt][idx_anno]["REVEL_score"])
		var_info.append(Allele_CSQ_dict[alt][idx_anno]["MPC_score"])
		var_info.append(Allele_CSQ_dict[alt][idx_anno]["MVP2_rankscore"])
		for j, trio in enumerate(Trios):
			#print(trio.Proband.ID)
			prob, fa, mo, sibs = trio.Proband, trio.Father, trio.Mother, trio.Siblings
			gt_prob, gt_fa, gt_mo = self.GenotypeQC(fmt, gts[prob.index]), self.GenotypeQC(fmt, gts[fa.index]), self.GenotypeQC(fmt, gts[mo.index])
			GT_prob, GT_fa, GT_mo = gts[prob.index], gts[fa.index], gts[mo.index]
			gt_sibs = [self.GenotypeQC(fmt, gts[x.index]) for x in sibs]
			GT_sibs = [gts[x.index] for x in sibs]
			#print (gt_prob, gt_fa, gt_mo)
			#if gt_prob == [1,1]:
			
			if trio.Proband.sampleID == "SP0015755":
				print(var_k, prob.index, gts[prob.index], gts[fa.index], gts[mo.index])
			
			if gt_prob != [0,0] and gt_prob != False:
				#print(gt_prob, gt_fa, gt_mo)
				pass

			if gt_prob == False or gt_fa == False or gt_mo == False: 
				#print (gt_prob, gt_fa, gt_mo)
				continue # Failed QC
			elif ( (gt_prob[1] not in [0, i+1]) or (gt_fa[1] not in [0, i+1]) or (gt_mo[1] not in [0, i+1]) ) or (gt_prob[1] == 0 and gt_fa[1] == 0 and gt_mo[1] == 0):
				#print (gt_prob, gt_fa, gt_mo)
				continue # Not this allele 
			sib_fail_qc = False
			for gt in gt_sibs:
				if gt == False:
					sib_fail_qc = True
			if sib_fail_qc:
				continue
			elif (gt_prob[0] not in gt_fa or gt_prob[1] not in gt_mo) and (gt_prob[1] not in gt_fa or gt_prob[0] not in gt_mo):
				#print ("ME", gt_prob, gt_fa, gt_mo)
				continue # Mendelian Error
			else:
				gt_prob, gt_fa, gt_mo = self.gt_recode(gt_prob), self.gt_recode(gt_fa), self.gt_recode(gt_mo)
				print(gt_prob, gt_fa, gt_mo)
				gt_sibs = [self.gt_recode(gt) for gt in gt_sibs]
				#Gene_Fam_dat[trio.FamID][Vartype].append([var_k, gt_prob, gt_fa, gt_mo])
				Gene_Fam_dat[Vartype][trio.FamID].append([var_k, var_info, gt_prob, gt_fa, gt_mo, gt_sibs, GT_prob, GT_fa, GT_mo, GT_sibs])
				#print(gt_prob, gt_fa, gt_mo)
		return Gene_Fam_dat
	
	def gt_recode(self, gt):
		if gt[0] != 0 :
			gt[0] = 1
		if gt[1] != 0 :
			gt[1] = 1
		return gt
	def Modify(self, GT):
		GT = GT.split(":")
		return GT[0] + ":" + GT[4]

	def Phasing_N_Count(self, Gene_Fam_dat, Trios, OutFil):
		res = {}
		for t in self.C:
			N_hom = 0
			N_chet = 0
			N_hom_chet = 0
			N_haps = 0
			N_cant_phase = 0
			cant_phase_fam = []
			N_more_than_three = 0
			for i, trio in enumerate(Trios):
				Nindv = len(trio.Siblings) + 1
				variants_in_fam = Gene_Fam_dat[t][trio.FamID] #list of variants in this gene in this fam
				#for item in variants_in_fam
				if len(variants_in_fam) == 1: #only 1 variant
					var_k, var_info, gt_pro, gt_fa, gt_mo, gt_sibs, GT_prob, GT_fa, GT_mo, GT_sibs = variants_in_fam[0]
					GT_prob = self.Modify(GT_prob)
					GT_fa = self.Modify(GT_fa)
					GT_mo = self.Modify(GT_mo)
					GT_sibs = [self.Modify(X) for X in GT_sibs]
					N_haps += sum(gt_fa + gt_mo)
					OutParents = False
					for i, (gt, GT, SPID) in enumerate(zip([gt_pro] + gt_sibs, [GT_prob] + GT_sibs, [trio.Proband.sampleID] + [x.sampleID for x in trio.Siblings])):
						if gt == [1,1]:
							OutParents = True
							N_hom += 1
							if i == 0:
								OutFil.writerow(var_info + [trio.FamID, SPID, GT, "Proband", Nindv])
							else:
								OutFil.writerow(var_info + [trio.FamID, SPID, GT, "Sibling", Nindv])
					if OutParents:
						OutFil.writerow(var_info + [trio.FamID, trio.Father.sampleID, GT_fa, "Father", Nindv])
						OutFil.writerow(var_info + [trio.FamID, trio.Mother.sampleID, GT_mo, "Mother", Nindv])
				elif len(variants_in_fam) == 2: # 2 variants 
					v1, var_info1, gt_p1, gt_f1, gt_m1, gt_sibs1, GT_prob1, GT_fa1, GT_mo1, GT_sibs1 = variants_in_fam[0]
					v2, var_info2, gt_p2, gt_f2, gt_m2, gt_sibs2, GT_prob2, GT_fa2, GT_mo2, GT_sibs2 = variants_in_fam[1]
					gts1 = zip([gt_p1] + gt_sibs1, [GT_prob1] + GT_sibs1, [trio.Proband.sampleID] + [x.sampleID for x in trio.Siblings])
					gts2 = zip([gt_p2] + gt_sibs2, [GT_prob2] + GT_sibs2, [trio.Proband.sampleID] + [x.sampleID for x in trio.Siblings])
					if (gt_f1 == [0,0] and gt_m1 == [0,1] and gt_f2 == [0,1] and gt_m2 == [0,0]) or (gt_f1 == [0,0] and gt_m1 == [0,1] and gt_f2 == [0,1] and gt_m2 == [0,0]):
						# 0/0 0/1 -> 0/1 
						# 0/1 0/0 -> 0/1
						N_haps += 2
						OutParents = False
						for i, (gt1, gt2) in enumerate(zip(gts1 ,gts2)):
							gt1, GT1, SPID1 = gt1
							gt2, GT2, SPID2 = gt2
							if gt1 == [0,1] and gt2 == [0,1]:
								N_chet += 1
								OutParents = True
								GT1 = self.Modify(GT1)
								GT2 = self.Modify(GT2)
								if i == 0:
									OutFil.writerow(var_info1 + [trio.FamID, SPID1, GT1, "Proband", Nindv])
									OutFil.writerow(var_info2 + [trio.FamID, SPID2, GT2, "Proband", Nindv])
								else:
									OutFil.writerow(var_info1 + [trio.FamID, SPID1, GT1, "Sibling", Nindv])
									OutFil.writerow(var_info2 + [trio.FamID, SPID2, GT2, "Sibling", Nindv])
						if OutParents:
							GT_fa1 = self.Modify(GT_fa1)
							GT_fa2 = self.Modify(GT_fa2)
							GT_mo1 = self.Modify(GT_mo1)
							GT_mo2 = self.Modify(GT_mo2)
							OutFil.writerow(var_info1 + [trio.FamID, trio.Father.sampleID, GT_fa1, "Father", Nindv])
							OutFil.writerow(var_info2 + [trio.FamID, trio.Father.sampleID, GT_fa2, "Father", Nindv])
							OutFil.writerow(var_info1 + [trio.FamID, trio.Mother.sampleID, GT_mo1, "Mother", Nindv])
							OutFil.writerow(var_info2 + [trio.FamID, trio.Mother.sampleID, GT_mo2, "Mother", Nindv])
					elif (gt_f1 == [0,1] and gt_m1 == [0,1] and gt_p1 == [0,1]) or (gt_f2 == [0,1] and gt_m2 == [0,1] and gt_p2 == [0,1]):
						# Unable to phase
						N_cant_phase += 1
						N_haps += 4
						cant_phase_fam.append(trio.FamID)
						OutParents = False
						#for gt1, gt2 in zip(gts1, gts2):
						for i, (gt1, gt2) in enumerate(zip(gts1 ,gts2)):
							gt1, GT1, SPID1 = gt1
							gt2, GT2, SPID2 = gt2
							if gt1 == [0,1] and gt2 == [0,1]:
								N_chet += 1
								OutParents = True
								GT1 = self.Modify(GT1)
								GT2 = self.Modify(GT2)
								if i == 0:
									OutFil.writerow(var_info1 + [trio.FamID, SPID1, GT1, "Proband", Nindv])
									OutFil.writerow(var_info2 + [trio.FamID, SPID2, GT2, "Proband", Nindv])
								else:
									OutFil.writerow(var_info1 + [trio.FamID, SPID1, GT1, "Sibling", Nindv])
									OutFil.writerow(var_info2 + [trio.FamID, SPID2, GT2, "Sibling", Nindv])
						if OutParents:
							GT_fa1 = self.Modify(GT_fa1)
							GT_fa2 = self.Modify(GT_fa2)
							GT_mo1 = self.Modify(GT_mo1)
							GT_mo2 = self.Modify(GT_mo2)
							OutFil.writerow(var_info1 + [trio.FamID, trio.Father.sampleID, GT_fa, "Father", Nindv])
							OutFil.writerow(var_info2 + [trio.FamID, trio.Father.sampleID, GT_fa, "Father", Nindv])
							OutFil.writerow(var_info1 + [trio.FamID, trio.Father.sampleID, GT_mo, "Mother", Nindv])
							OutFil.writerow(var_info2 + [trio.FamID, trio.Father.sampleID, GT_mo, "Mother", Nindv])

				elif len(variants_in_fam) >= 2: # more than 2 variants
					N_more_than_three += 1
			res[t] = (N_hom, N_chet, N_haps, N_cant_phase, cant_phase_fam, N_more_than_three)
		return res

	def LookUpBiallic(self, i, Vartype, fmt, gts, Trios):
		for j, trio in enumerate(Trios):
			prob, fa, mo = trio.Proband, trio.Father, trio.Mother
			#print(fmt, gts[prob.index])
			gt_prob, gt_fa, gt_mo = self.GenotypeQC(fmt, gts[prob.index]), self.GenotypeQC(fmt, gts[fa.index]), self.GenotypeQC(fmt, gts[mo.index])
			if gt_prob == False or gt_fa == False or gt_mo == False:
				continue # gt failed QC
			elif (gt_prob[0] not in gt_fa or gt_prob[1] not in gt_mo) or (gt_prob[1] not in gt_fa or gt_prob[0] not in gt_mo):
				continue # mendelian error
			else:
				# Phasing
				if gt_prob[1] == i+1 and gt_prob[0] == i+1: # Hom
					Trios[j].pro_haps[Vartype] = [1,1]
					#print("12", Trios[j].pro_haps[Vartype], gt_prob, gt_fa, gt_mo)
					if gt_fa[0] == i+1 and gt_fa[1] == i+1:
						Trios[j].fa_haps[Vartype] = [1,1]
					else:
						Trios[j].fa_haps[Vartype][0] = 1
					if gt_mo[0] == i+1 and gt_mo[1] == i+1:
						Trios[j].mo_haps[Vartype] = [1,1]
					else:
						Trios[j].mo_haps[Vartype][0] = 1

				elif gt_prob[1] == gt_fa[1] and gt_mo[1] == 0 and gt_prob[1] == i+1 : #paternal transmitted
					if gt_fa[0] == i+1: # transmitted from hom parternal
						Trios[j].fa_haps[Vartype] = [1,1]
					else: # transmitted from het paternal
						Trios[j].fa_haps[Vartype][0] = 1
					Trios[j].pro_haps[Vartype][0] = 1
					
				elif gt_prob[1] == gt_mo[1] and gt_fa[1] == 0 and gt_prob[1] == i+1: #maternal transmitted
					if gt_mo[0] == i+1: # transmitted from hom maternal
						Trios[j].mo_haps[Vartype] = [1,1]
					else: # transmitted from het maternal
						Trios[j].mo_haps[Vartype][0] = i+1
					Trios[j].pro_haps[Vartype][1] = 1

				elif gt_prob[1] == 0: # proband 0/0
					if gt_fa[1] == 1: # father has one hap
						if Trios[j].fa_haps[Vartype][0] == 1:
							Trios[j].fa_haps[Vartype][1] = 0
						elif Trios[j].fa_haps[Vartype][0] == 0:
							Trios[j].fa_haps[Vartype][1] = 1
					if gt_mo[1] == 1: # mother has one hap
						if Trios[j].mo_haps[Vartype][0] == 1:
							Trios[j].mo_haps[Vartype][1] = 0
						elif Trios[j].mo_haps[Vartype][0] == 0:
							Trios[j].mo_haps[Vartype][1] = 1

	def LookUpBiallicLGD_DMIS(self, i, Vartype, fmt, gts, Trios):
		for j, trio in enumerate(Trios):
			prob, fa, mo = trio.Proband, trio.Father, trio.Mother
			#print(fmt, gts[prob.index])
			gt_prob, gt_fa, gt_mo = self.GenotypeQC(fmt, gts[prob.index]), self.GenotypeQC(fmt, gts[fa.index]), self.GenotypeQC(fmt, gts[mo.index])
			if gt_prob == False or gt_fa == False or gt_mo == False:
				continue # gt failed QC
			elif (gt_prob[0] not in gt_fa and gt_prob[1] not in gt_mo) or (gt_prob[1] not in gt_fa and gt_prob[0] not in gt_mo):
				continue # mendelian error
			else:
				# Phasing
				if gt_prob[1] == i+1 and gt_prob[0] == i+1: # Hom
					Trios[j].pro_haps[Vartype] = [1,1]
					#print("12", Trios[j].pro_haps[Vartype], gt_prob, gt_fa, gt_mo)
					if gt_fa[0] == i+1 and gt_fa[1] == i+1:
						Trios[j].fa_haps[Vartype] = [1,1]
					else:
						Trios[j].fa_haps[Vartype][0] = 1
					if gt_mo[0] == i+1 and gt_mo[1] == i+1:
						Trios[j].mo_haps[Vartype] = [1,1]
					else:
						Trios[j].mo_haps[Vartype][0] = 1

				elif gt_prob[1] == gt_fa[1] and gt_mo[1] == 0 and gt_prob[1] == i+1 : #paternal transmitted
					if gt_fa[0] == i+1: # transmitted from hom parternal
						Trios[j].fa_haps[Vartype] = [1,1]
					else: # transmitted from het paternal
						Trios[j].fa_haps[Vartype][0] = 1
					Trios[j].pro_haps[Vartype][0] = 1
					
				#elif gt_fa[1] == 0 and gt_prob[1] == 0:
				#	Trios[j].fa_haps[Vartype][1] = 1


				elif gt_prob[1] == gt_mo[1] and gt_fa[1] == 0 and gt_prob[1] == i+1: #maternal transmitted
					if gt_mo[0] == i+1: # transmitted from hom maternal
						Trios[j].mo_haps[Vartype] = [1,1]
					else: # transmitted from het maternal
						Trios[j].mo_haps[Vartype][0] = i+1
					Trios[j].pro_haps[Vartype][1] = 1

				elif gt_prob[1] == 0: # proband 0/0
					if gt_fa[1] == 1: # father has one hap
						if Trios[j].fa_haps[Vartype][0] == 1:
							Trios[j].fa_haps[Vartype][1] = 0
						elif Trios[j].fa_haps[Vartype][0] == 0:
							Trios[j].fa_haps[Vartype][1] = 1
					if gt_mo[1] == 1: # mother has one hap
						if Trios[j].mo_haps[Vartype][0] == 1:
							Trios[j].mo_haps[Vartype][1] = 0
						elif Trios[j].mo_haps[Vartype][0] == 0:
							Trios[j].mo_haps[Vartype][1] = 1

class Sample:
	def __init__(self, row):
		#self.FamID = row["FamID"]
		#self.sampleID = row["SampleID"]
		#self.Father = row["Paternal"]
		#self.Mother = row["Maternal"]
		#self.Sex = row["Sex"]
		#self.Affected = row["Affected"]
		self.FamID = row[0]
		self.sampleID = row[1]
		self.Father = row[2]
		self.Mother = row[3]
		self.Sex = row[4]
		self.Affected = row[5]
		self.index = row[-1]
	def show(self):
		print(self.FamID, self.sampleID, self.Father, self.Mother, self.Sex, self.Affected)
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
		self.pro_haps = {} 
		self.fa_haps = {}
		self.mo_haps = {}
	def show(self):
		print("FamID:{} Proband:{} Father:{} Mother:{} Siblings:{}".format(
			self.FamID, self.Proband.sampleID, self.Father.sampleID, self.Mother.sampleID, ", ".join(
				[x.sampleID for x in self.Siblings])))

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
			if info["gene_name"] not in Genes:
				Genes[info["gene_name"]] = GTFRecord(llist[0], llist[1], llist[2], llist[3], llist[4], llist[6], info)
				Transcripts[info["gene_name"]] = []
			Transcripts[info["gene_name"]].append(GTFRecord(llist[0], llist[1], llist[2], llist[3], llist[4], llist[6], info))
	return Genes, Transcripts 
# 
def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument("--chr", type=str, required=True, help="<Required> Chromesome")
	parser.add_argument("--gene", type=str, required=True, help="<Required> Gene")
	parser.add_argument("--af", type=float, help="<Required> Allele Freq")
	args = parser.parse_args()
	#if args.out == None:
	#	args.out = "test.out.vcf"
	return args

def main():
	args = GetOptions()
	ins = RecessiveModel(args.af)
	Chr = args.chr
	Gene = args.gene
	#GenotypeFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K/VCF/TrioVCF/Genotypes/SPARK30K.TrioSamples.Chr{}.vcf.gz".format(Chr)
	#VEPFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K/VCF/TrioVCF/sites/SPARK30K.TrioSamples.Chr{}.vep.mappability.vcf.gz".format(Chr)
	#VEPFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K/VCF/TrioVCF/sites/Annotated2/SPARK30K.TrioSamples.Chr{}.vep.vcf.gz".format(Chr)
	#AFFil = "/home/local/users/jw/Genetics_Projects/SPARK/spark_genomics/dat/SPARK30K.TrioSamples.Chr{}.eurAF.vcf.gz".format(Chr)
	GenotypeFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K_07/VCF/GenotypesSplitbyChr/GATK4_20190729.chr{}.vcf.gz".format(Chr)
	VEPFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K_07/VCF/SitesSplitbyChr/annotated/GATK4_20190729.chr{}.mappability.vcf.gz".format(Chr)
	AFFil = "/home/local/users/jw/Genetics_Projects/SPARK/30K_07/recessive/AF/GATK4_20190729.chr{}.eurAF.vcf.gz".format(Chr)
	genecode = "/home/local/users/jw/vep_data/homo_sapiens/GeneCodeV29/CHRs/genecodev29.{}.gtf".format(Chr)
	ins.Recessive(Chr, Gene, GenotypeFil ,VEPFil, AFFil, genecode)

if __name__=='__main__':
	main()
