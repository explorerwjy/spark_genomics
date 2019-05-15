import sklearn as sk
from sklearn import decomposition
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import gzip as gz
import scipy
from scipy import stats
import math
import random
from scipy.optimize import minimize
from scipy.special import factorial
import scipy.stats as ss
import scipy.optimize as so
import csv
from itertools import groupby, cycle
from operator import itemgetter


class AncestryPCA:
    def __init__(self, referenceVCF, referencePanel, VCF):
        self.referenceVCF = referenceVCF
        self.referencePanel = referencePanel
        self.VCF = VCF
        self.Colors = {'SAS':"green" , 'EAS':"blue", 'AMR':"orange", 'AFR':"purple", 'EUR':"red", 'AJ':"pink", 'DOMI':"yellow",'CASE':"grey"}
    def LoadPanel(self):
        df = pd.read_csv(self.referencePanel, delimiter="\t")
        Populations = df["super_pop"].unique()
        self.Pop = {}
        for pop in Populations:
            self.Pop[pop] = df[df["super_pop"]==pop]["sample"].values
    def LoadVCF(self, vcf_file):
        hand = gz.open(vcf_file, 'rt')
        d = {}
        for l in hand:
            if l.startswith("##"):
                continue
            elif l.startswith("#"):
                indvs = l.strip().split("\t")[9:]
                for indv in indvs:
                    d[indv] = []
            else:
                if l.strip().split("\t")[0] in ["Y", "chrY", "chrX", "X"]:
                    continue
                GTs = [self.convertGT(GT) for GT in l.strip().split("\t")[9:]]
                for indv, GT in zip(indvs, GTs):
                    d[indv].append(GT)
        df = pd.DataFrame(data=d)
        return df.transpose()
    def convertGT(self, GT):
        GT = GT.split(":")[0]
        if "." in GT:
            return -9
        elif "/" in GT:
            GT = map(int, GT.split("/"))
            return sum(GT)
        elif "|" in GT:
            GT = map(int, GT.split("|"))
            return sum(GT)
        else:
            return -9
    def pca(self):
        self.LoadPanel()
        df = self.LoadVCF(self.referenceVCF)
        a, b = df.shape
        print ("Reference Panel has {} individuals {} SNPs".format(a,b))
        self.model = decomposition.PCA()
        self.model.fit(df)

    
def simGT(p, N):
    exp = 0
    for i in range(N):
        if random.uniform(0,1) <= p:
            exp += 1
    return exp

def onerun(NumofTotalHap, NumofRareHap, NumofChildsDict):
    dat = np.concatenate( (np.zeros(NumofTotalHap-NumofRareHap), np.ones(NumofRareHap) ), axis=0)
    np.random.shuffle(dat)
    NumofChilds = sorted(NumofChildsDict.items(), key=lambda x: x[0])
    exp = 0
    start = 0
    for i, count in enumerate(NumofChilds):
        for j in range(start, start + count[1]*4, 4):
            f1, f2, m1, m2 = dat[j: j+4]
            if f1+f2 == 1 and m1 + m2 == 1:
                exp += simGT(0.25, count[0])#0.25 * count[0]
            elif (f1+f2 == 1 and m1 + m2 == 2) or (f1+f2 == 2 and m1 + m2 == 1):
                exp += simGT(0.5, count[0])
            elif f1+f2 == 2 and m1 + m2 == 2:
                exp += simGT(1, count[0])
        start = j
    return exp

def onerun2(NumofTotalHap, NumofRareHap, NumofChildsDict, frac=0.1):
    dat = np.concatenate( (np.zeros(NumofTotalHap-NumofRareHap), np.ones(NumofRareHap) ), axis=0)
    np.random.shuffle(dat)
    NumofChilds = sorted(NumofChildsDict.items(), key=lambda x: x[0])
    exp = 0
    start = 0
    for i, count in enumerate(NumofChilds):
        for j in range(start, start + count[1]*4, 4):
            f1, f2, m1, m2 = dat[j: j+4]
            if f1+f2 == 1 and m1 + m2 == 1:
                exp += simGT(0.25, count[0])#0.25 * count[0]
            elif (f1+f2 == 1 and m1 + m2 == 2) or (f1+f2 == 2 and m1 + m2 == 1):
                exp += simGT(0.5, count[0])
            elif f1+f2 == 2 and m1 + m2 == 2:
                exp += simGT(1, count[0])
        start = j
    return exp

def Permutation(NumofTotalHap, NumofRareHap, NumofChildsDict, Nperm=20000):
    EXP = []
    for i in range(Nperm):
        exp = onerun(NumofTotalHap, NumofRareHap, NumofChildsDict)
        EXP.append(exp)
    EXP = np.array(EXP)
    return EXP

def PlotNFit(dat, af, mu = None, fit=False):
    bins=np.arange(min(dat), max(dat)+1)
    counts, bins = np.histogram(dat,bins=bins,density=1)
    count = dict(zip(bins,counts))
    x = np.arange(min(dat), max(dat))
    y = np.array([count[i] for i in x])
    plt.plot(x, y, 'bo', ms=8, label='permute pmf', color="red")
    plt.vlines(x, 0, y, colors='b', lw=5, alpha=0.5, color="red")
    if fit:
        mu = max(0, FitPoisson(dat))
    else:
        mu = np.mean(dat)
    x = np.arange(max(0, stats.poisson.ppf(0.0, mu)), max(1, stats.poisson.ppf(0.9999, mu)))
    plt.plot(x, stats.poisson.pmf(x, mu), 'bo', ms=8, label='poisson pmf')
    plt.vlines(x, 0, stats.poisson.pmf(x, mu), colors='b', lw=5, alpha=0.5)
    mean, var = np.mean(dat), np.var(dat)
    plt.title("AF={}, mean:{}, var:{}, lambda={}".format(af, round(mean,3), round(var,3), round(float(mu),3)))
    plt.legend(loc='best', frameon=False)
    plt.show()
    if fit:
        return mu
    else:
        return None

def NegBinom_Likelihood(P, x, neg=1):
    n=np.round(P[0]) #by definition, it should be an integer
    p=P[1]
    loc=np.round(P[2])
    return neg*(np.log(ss.nbinom.pmf(x, n, p, loc))).sum()
def FitNegBinom(dat):
    result=[]
    for i in range(20, 160): #in fact (80, 120) should probably be enough
        _=so.fmin(NegBinom_Likelihood, [i, 0.5, 0], args=(dat,-1), full_output=True, disp=False)
        result.append((_[1], _[0]))
    P2=sorted(result, key=lambda x: x[0])[0][1]
    return np.round(P2[0]), P2[1]

def plotscatter(E, O, title, xlim=None, ylim=None):
    #fig = plt.figure(figsize=(5,5), dpi=200)
    #fig.clear()

    plt.figure(dpi=120)
    plt.scatter(E, O, s=5)
    xlim = max(E) if xlim == None else xlim
    ylim = max(O) if ylim == None else ylim
    _max = max([xlim, ylim])
    plt.text(xlim*0.8, ylim*0.8, "Ratio:{0:.2f}".format(sum(O)/sum(E)))
    plt.plot((0, _max), (0, _max), color="black")
    plt.xlabel("Expected")
    plt.ylabel("Observed")
    plt.title(title)
    plt.xlim((0, xlim))
    plt.ylim((0, ylim))
    plt.show()
    plt.savefig("../Slides/{}.png".format("".join(title.split())), dpi=200, format="png")
    plt.clf()

def QQplot(pvalues, title="QQ plot"):
    pvalues.sort(reverse=True)
    Qvalues = []
    for x in pvalues:
        try:
            Qvalues.append(min(10, -math.log(x,10)))
        except:
            print(x)
    top = int(Qvalues[-1]) + 1
    NumTest = len(Qvalues)
    Qvalues = [0] * (19000-NumTest) + Qvalues
    Qexp = []
    for i in range(len(Qvalues)):
        Qexp.append(float(i+1)/NumTest)
    Qexp.sort(reverse=True)
    Qexp = [-1*math.log(x,10) for x in Qexp]
    plt.subplot()
    plt.scatter(Qexp, Qvalues, alpha=0.5)
    plt.plot([0, top], [0, top], ls="-")
    plt.title(title)
    plt.xlabel('Exp Q')
    plt.ylabel('Obs Q')
    plt.show()

def QQplot2(pvalues, title="QQ plot", threshold=5e-8):
    pvalues.sort(reverse=True)
    Qvalues = []
    for x in pvalues:
        try:
            Qvalues.append(min(10, -math.log(x,10)))
        except:
            print(x)
    top = int(Qvalues[-1]) + 1
    NumTest = len(Qvalues)
    Qexp = []
    for i in range(NumTest):
        Qexp.append(float(i+1)/NumTest)
    Qexp.sort(reverse=True)
    Qexp = [-math.log(x,10) for x in Qexp]
    #plt.subplot(dpi=120)
    plt.figure(figsize=(4, 4), dpi=100)
    plt.scatter(Qexp, Qvalues, alpha=0.5, s=2)
    plt.plot([0, top], [0, top], ls="--", color='black')
    if threshold != False:
        plt.axhline(y=-math.log(threshold, 10), linestyle="--", color="black")
    plt.title(title)
    plt.xlabel('Exp Q')
    plt.ylabel('Obs Q')
    plt.show()
    #return Qvalues, Qexp

def volcano(pvalues, effects, title="volcano plot"):
    Qvalues = [-1*math.log(x,10) for x in pvalues]
    effects = [math.log(x,2) for x in effects]
    plt.scatter(effects, Qvalues)
    plt.title(title)
    plt.show()
    
def GetPvalues(exp, obs):
    ratio, pvalue = [], []
    for e,o in zip(exp, obs):
        r = o/e if e != 0 else 'nan'
        p = stats.poisson.sf(obs-1,exp)
        ratio.append(r)
        pvalue.append(p)
    return ratio, pvalue

def _gen_data(df):
    """
    iterate over the files and yield chr, start, pvalue
    """
    for row in df.iterrows():
        #print(row)
        #print(row[1])
        yield row[1]["CHR"], row[1]["BP"], row[1]["P"]
        #yield row[1]["CHR"], row[1]["POS"], row[1]["P-value"]
            #yield toks[columns[0]], int(toks[columns[1]]), float(toks[columns[2]])
def cmp(a, b):
    if a > b:
        return 1
    elif a == b:
        return 0
    elif a < b:
        return -1
def chr_loc_cmp(alocs, blocs):
    return cmp(alocs[0], blocs[0]) or cmp(alocs[1], blocs[1])
def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K:
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K
def manhattan(df, image_path="./manhattan.png", no_log=False, colors="rgbk", title="manhattan", lines=False, ymax=10):
    xs = []
    ys = []
    cs = []
    colors = cycle(colors)
    xs_by_chr = {}
    last_x = 0
    #data = sorted(_gen_data(df), cmp=chr_loc_cmp)
    data = sorted(_gen_data(df), key=cmp_to_key(chr_loc_cmp))
    for seqid, rlist in groupby(data, key=itemgetter(0)):
        #color = colors.next()
        color = next(colors)
        rlist = list(rlist)
        region_xs = [last_x + r[1] for r in rlist]
        xs.extend(region_xs)
        ys.extend([r[2] for r in rlist])
        cs.extend([color] * len(rlist))
        xs_by_chr[seqid] = (region_xs[0] + region_xs[-1]) / 2
        # keep track so that chrs don't overlap.
        last_x = xs[-1]
    #xs_by_chr = [(k, xs_by_chr[k]) for k in sorted(xs_by_chr.keys(), cmp=chr_cmp)]
    xs_by_chr = [(k, xs_by_chr[k]) for k in sorted(xs_by_chr.keys(), key=cmp_to_key(cmp))]
    xs = np.array(xs)
    ys = np.array(ys) if no_log else -np.log10(ys)
    plt.close()
    f = plt.figure(figsize=(10,4), dpi=120)
    ax = f.add_axes((0.1, 0.09, 0.88, 0.85))
    if title is not None:
        plt.title(title)
    ax.set_ylabel('-log10(p-value)')
    if lines:
        ax.vlines(xs, 0, ys, colors=cs, alpha=0.5)
    else:
        ax.scatter(xs, ys, s=2, c=cs, alpha=0.8, edgecolors='none')
    # plot 0.05 line after multiple testing.
    #ax.axhline(y=-np.log10(0.05 / len(data)), color='0.5', linewidth=2)
    ax.axhline(y=-np.log10(5e-8), color='0.5', linewidth=2)
    plt.axis('tight')
    plt.xlim(0, xs[-1])
    plt.ylim(ymin=0)
    if ymax is not None: plt.ylim(ymax=ymax)
    plt.xticks([c[1] for c in xs_by_chr], [c[0] for c in xs_by_chr], rotation=0, size=8.5)
    #print >>sys.stderr, "saving to: %s" % image_path
    #plt.savefig(image_path)
    plt.show()
def get_filehandles(args):
    return (open(a) if a != "-" else sys.stdin for a in args)

def processPRS_PTDT(S):
    inp = open("/Users/jiayao/Work/spark/dat/30K/GWAS/spark_prs/plink.{}.profile".format(S), 'rt')
    out = csv.writer(open("/Users/jiayao/Work/spark/dat/30K/GWAS/spark_prs/plink.{}.profile.tsv".format(S), 'wt'), delimiter="\t")
    for l in inp:
        out.writerow(l.split())
    PRS = pd.read_csv("/Users/jiayao/Work/spark/dat/30K/GWAS/spark_prs/plink.{}.profile.tsv".format(S), delimiter="\t")
    ID2Score = dict(zip(PRS["IID"].values, PRS["SCORE"].values))
    Case = PRS[PRS["PHENO"]==2]["SCORE"].values
    NonCase = PRS[PRS["PHENO"]==1]["SCORE"].values
    #plt.hist(Case, color="red", alpha=0.5, bins=50, normed=1)
    #plt.hist(NonCase, color="blue", alpha=0.5, bins=50, normed=1)
    #plt.show()
    mean_case = np.mean(Case)
    mean_control = np.mean(NonCase)
    #print(mean_case, mean_control)
    FamDat = pd.read_csv("/Users/jiayao/Work/spark/dat/30K/GWAS/spark_prs/GenoHQ.fam.tsv", delimiter="\t", header=None)
    FamDat.columns = ["FamID", "SampleID", "FatherID", "MotherID", "Gender", "Pheno"]
    MidPrs = []
    ProbPrs = []
    for row in FamDat.iterrows():
        row = row[1]
        try:
            if (row["Pheno"] == 2) and (row["FatherID"]!='0') and (row["MotherID"]!='0'):
                prob_prs = ID2Score[row["SampleID"]]
                FaPrs = ID2Score[row["FatherID"]]
                MoPrs = ID2Score[row["MotherID"]]
                mid_prs = (FaPrs+MoPrs)/2
                MidPrs.append(mid_prs)
                ProbPrs.append(prob_prs)
        except:
            continue
    DEVs = []
    SD_PRS_MP = np.std(MidPrs)
    for prob_prs, mid_prs in zip(ProbPrs, MidPrs):
        pTDT_dev = (prob_prs - mid_prs)/SD_PRS_MP
        DEVs.append(pTDT_dev)
    prob_dev = np.mean(DEVs)
    prob_std = np.std(DEVs) / math.sqrt(len(DEVs))
    t, prob_p = stats.ttest_1samp(DEVs, 0)
    MidPrs = []
    SibPrs = []
    for row in FamDat.iterrows():
        row = row[1]
        try:
            if (row["Pheno"] == 1) and (row["FatherID"]!='0') and (row["MotherID"]!='0'):
                sib_prs = ID2Score[row["SampleID"]]
                FaPrs = ID2Score[row["FatherID"]]
                MoPrs = ID2Score[row["MotherID"]]
                mid_prs = (FaPrs+MoPrs)/2
                MidPrs.append(mid_prs)
                SibPrs.append(sib_prs)
        except:
            continue
            print(row)
    SD_PRS_MP = np.std(MidPrs)
    #print(SD_PRS_MP)
    DEVs = []
    for sib_prs, mid_prs in zip(SibPrs, MidPrs):
        pTDT_dev = (sib_prs - mid_prs)/SD_PRS_MP
        DEVs.append(pTDT_dev)
    #plt.hist(DEVs, bins=100)
    #plt.show()
    sib_dev = np.mean(DEVs)
    sib_std = np.std(DEVs) / math.sqrt(len(DEVs))
    t, sib_p = stats.ttest_1samp(DEVs, 0)
    return (prob_dev, prob_std, prob_p), (sib_dev, sib_std, sib_p)

def processPRS_PTDT_stratified(S, group1, group2):
    PRS = pd.read_csv("/Users/jiayao/Work/spark/dat/30K/GWAS/spark_prs/plink.{}.profile.tsv".format(S), delimiter="\t")
    ID2Score = dict(zip(PRS["IID"].values, PRS["SCORE"].values))
    FamDat = pd.read_csv("/Users/jiayao/Work/spark/dat/30K/GWAS/spark_prs/GenoHQ.fam.tsv", delimiter="\t", header=None)
    FamDat.columns = ["FamID", "SampleID", "FatherID", "MotherID", "Gender", "Pheno"]
    MidPrs1 = []
    ProbPrs1 = []
    MidPrs2 = []
    ProbPrs2 = []
    for row in FamDat.iterrows():
        row = row[1]
        try:
            prob_prs = ID2Score[row["SampleID"]]
            FaPrs = ID2Score[row["FatherID"]]
            MoPrs = ID2Score[row["MotherID"]]
            mid_prs = (FaPrs+MoPrs)/2
            if row["SampleID"] in group1:
                MidPrs1.append(mid_prs)
                ProbPrs1.append(prob_prs)
            elif row["SampleID"] in group2:
                MidPrs2.append(mid_prs)
                ProbPrs2.append(prob_prs)
        except:
            continue
    g1_devs = []
    SD_PRS_MP1 = np.std(MidPrs1)
    for prob_prs, mid_prs in zip(ProbPrs1, MidPrs1):
        pTDT_dev = (prob_prs - mid_prs)/SD_PRS_MP1
        g1_devs.append(pTDT_dev)
    g1_mean_dev = np.mean(g1_devs)
    g1_std = np.std(g1_devs) / math.sqrt(len(g1_devs))
    t1, p1 = stats.ttest_1samp(g1_devs, 0)

    g2_devs = []
    SD_PRS_MP2 = np.std(MidPrs2)
    for prob_prs, mid_prs in zip(ProbPrs2, MidPrs2):
        pTDT_dev = (prob_prs - mid_prs)/SD_PRS_MP2
        g2_devs.append(pTDT_dev)
    g2_mean_dev = np.mean(g2_devs)
    g2_std = np.std(g2_devs) / math.sqrt(len(g2_devs))
    t2, p2 = stats.ttest_1samp(g2_devs, 0)
    t, p = stats.ttest_ind(g1_devs, g2_devs)
    return (g1_mean_dev, g1_std, p1), (g2_mean_dev, g2_std, p2), p
