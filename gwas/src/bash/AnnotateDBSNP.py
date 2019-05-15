import pysam
import sys
import csv

dbsnp="/home/local/users/jw/resources/references/GRCh38/dbsnp_146.hg38.vcf.gz"
tbx = pysam.TabixFile(dbsnp)

reader = csv.reader(open(sys.argv[1], 'rt'), delimiter="\t")
writer = csv.writer(open("ipsych.common.scores", 'wt'), delimiter="\t")
idx_chr = 1 #header.index("CHR")
idx_pos = 2 #header.index("POS")
idx_ref = 3 #header.index("A1")
idx_alt = 4 #header.index("A2")
idx_z = 9
Count1 = 0
Count2 = 0
for row in reader:
	CHR = row[idx_chr]
	POS = int(row[idx_pos])
	REF = row[idx_ref]
	ALT = row[idx_alt]
	if len(REF) > len(ALT):
		offset = len(REF) - len(ALT)
		_pos = POS + offset
	else:
		_pos = POS
	#print(CHR, POS, _pos, REF, ALT)
	for record in tbx.fetch(CHR, POS-1, POS+1):
		record = record.split("\t")
		#print(CHR, POS,_pos, REF, ALT, record)
		_alts = record[4].split(",")
		for alt in _alts:
			#if (record[0] == CHR and int(record[1]) == _pos and record[3] == REF and alt == ALT) or (record[0] == CHR and int(record[1]) == _pos and alt == REF and record[3] == ALT):
			#	snpid = record[2]
			#	continue
			#print(snpid)
			if record[0] == CHR and int(record[1]) == _pos:
				if record[3] == REF and alt == ALT:
					z_score = float(row[idx_z])
					Count1 += 1
				elif record[3] == ALT and alt == REF:
					z_score = -float(row[idx_z])
					Count2 += 1
	new_row = [row[0], record[3], z_score]
	writer.writerow(new_row)

print(Count1, Count2)
