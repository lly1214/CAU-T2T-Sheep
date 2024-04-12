#!/usr/bin/env python

import sys, os 
from Bio import SeqIO

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(SCRIPT_PATH + '/../lib/')
from kit import *

def fa2dict(path):
    seqs = {}
    handle = gzip.open(path, "rt") if is_gzfile(path) else open(path, 'r')
    for record in SeqIO.parse(handle, "fasta"):
        seqs[record.id] = str(record.seq).upper()
    return seqs

def save_seq(name, seq, outfile):
	with open(outfile, 'w') as tmp:
		tmp.write(">%s\n%s\n" % (name, seq))
	return outfile

if len(sys.argv) != 4:
	print ("python %s map_cmd ref query" % sys.argv[0], file=sys.stderr)
	sys.exit(1)


map_cmd = sys.argv[1]
ref = fa2dict(sys.argv[2])
query = fa2dict(sys.argv[3])
outfile = sys.argv[3] + '.nextp2w'

pg = None
alns = []
for name, seq in ref.items():
	assert not name.endswith("np12")
	reffile = save_seq(name, seq, outfile + '.ref.fa')
	if name not in query:
		name = name + '_np12'
	seq = query[name]
	queryfile = save_seq(name, seq, outfile + '.query.fa')
	cmd = map_cmd % (reffile, queryfile)
	is_suc, out, err = run_cmd(cmd)
	if not is_suc:
		print("Failed run %s, stdout:%s stderr:%s" % (cmd, out, err), file=sys.stderr)
		sys.exit(1)

	for line in out.split("\n"):
		if line.startswith("@SQ"):
			print(line)
		elif line.startswith("@PG"):
			pg = line
		else:
			alns.append(line)
	os.remove(reffile)
	os.remove(queryfile)

print(pg)
for aln in alns:
	print(aln)
