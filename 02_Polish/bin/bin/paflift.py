#!/usr/bin/env python

import re, sys
from Bio import SeqIO

def fa2dict(path):
	seqs = {}
	handle = gzip.open(path, "rt") if path.endswith(('gz', 'gzip', 'GZ', 'GZIP')) else open(path, 'r')
	for record in SeqIO.parse(handle, "fasta"):
		seqs[record.id] = str(record.seq).upper()
	return seqs

def get_cigar(alns):
	for aln in alns:
		if aln.startswith("cg:Z:"):
			return aln.strip("cg:Z:")
	return None

def read_paf(path):
	pafs = {}
	with open(path) as IN:
		for line in IN:
			lines = line.strip().split()
			cigar = get_cigar(lines)
			assert cigar != None
			qname, qlen, qs, qe, std, rname, rlen, rs, re, mch, aln, qscore = lines[:12]
			qlen, qs, qe, rlen, rs, re, mch, aln, qscore = map(int, (qlen, qs, qe, rlen, rs, re, mch, aln, qscore))
			if qscore <= 0 or std != '+':
				continue
			if qname not in pafs:
				pafs[qname] = {}
			if rname not in pafs[qname]:
				pafs[qname][rname] = []
			pafs[qname][rname].append((qlen, qs, qe, rlen, rs, re, mch, aln, qscore, cigar))
	for qname in pafs:
		for rname in pafs[qname]:
			pafs[qname][rname] = sorted(pafs[qname][rname], key = lambda x: -x[6])
	return pafs

def parse_name(name):
	names = name.split("_")
	return names[0], int(names[1]), int(names[2]), int(names[4]), int(names[5])


if len(sys.argv) != 4:
	print ("python %s asm.fa asm.paf ref.fa > result.fbed" % sys.argv[0])
	sys.exit(1)

seqs = fa2dict(sys.argv[1])
pafs = read_paf(sys.argv[2])
ref = fa2dict(sys.argv[3])


for qname in pafs:
	for rname in pafs[qname]:
		qsca, qsta, qend, qtsta, qtend = parse_name(qname)
		rsca, rsta, rend, rtsta, rtend = parse_name(rname)
		if qsca != rsca or not rtsta <= qtsta <= qtend <= rtend:
			continue
		s, e = qsta - 100, qend + 100
		for qlen, qsta, qend, rlen, rsta, rend, mch, aln, qscore, cigar in pafs[qname][rname]:
			qsta += qtsta
			qend += qtsta
			if not (qsta <= s <= qend or qsta <= e <= qend):
				continue
			g = re.findall(r'(\d+)([SMIDH])', cigar)
			for n, c in g:
				n = int(n)
				if c =='M':
					for i in range(n):
						if s <= qsta + i <= e:
							print("%s\t%d\t%s\t%d" % (qname, qsta + i, rname, rsta + i))
					rsta += n
					qsta += n
				elif c == 'S':
					qsta += n
				elif c == 'I':
					for i in range(n):
						if s <= qsta + i <= e:
							print("%s\t%d\t%s\t%d" % (qname, qsta + i, rname, rsta))
					qsta += n
				elif c == 'D':
					for i in range(n):
						if s <= qsta <= e:
							print("%s\t%d\t%s\t%d" % (qname, qsta, rname, rsta + i))
					rsta += n
				else:
					print ('ERROR:', n, c)
					sys.exit(1)


