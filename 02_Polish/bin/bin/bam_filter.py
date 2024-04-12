#!/usr/bin/env python
import sys
import pysam
import argparse
import gzip
from Bio import SeqIO

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
	pass

def get_mlen(read):
	nm = read.get_tag('NM')
	mlen = blen = 0 
	cigartuples = read.cigartuples
	for op, lens in cigartuples:
		if op in [0, 7, 8]: 
			mlen += lens
			blen += lens
		elif op in [1, 2]: 
			blen += lens
			nm -= lens
	return mlen - nm

def filt_map_reads(ref_len, map_reads, out_bam, min_map_depth, flank, win=50, depth_thread=30):
	depth = [0] * (int(ref_len/win) + 1)
	map_reads.sort(key = lambda x: (x[1], x[0].mapping_quality, x[0].is_secondary), reverse=True)
	filt_map_reads = []
	for read, _ in map_reads:
		s, e = int(read.reference_start/win) + 1, int(read.reference_end/win)
		min_depth = sys.maxsize
		for i in range(s, e + 1):
			if depth[i] < min_depth:
				min_depth = depth[i]
		if min_depth > depth_thread:
			continue
		for i in range(s, e + 1):
			depth[i] += 1
		filt_map_reads.append(read)
	check_offset = int(flank/win)
	for d in depth[check_offset : len(depth) - check_offset]:
		if d <= min_map_depth:
			return False
	filt_map_reads.sort(key = lambda x: (x.reference_start, x.reference_end))
	for read in filt_map_reads:
		out_bam.write(read)
	return True

def fa2dict(path):
    seqs = {}
    handle = gzip.open(path, "rt") if path.endswith(('gz', 'gzip', 'GZ', 'GZIP')) else open(path, 'r')
    for record in SeqIO.parse(handle, "fasta"):
        seqs[record.id] = str(record.seq).upper()
    return seqs

def is_head(qs, qe, ql, rs, re, rl, max_clip_len):
	return rs < max_clip_len and ql - qe < max_clip_len

def is_tail(qs, qe, ql, rs, re, rl, max_clip_len):
	return qs < max_clip_len and rl - re < max_clip_len
		
def main(args):
	ref = None
	map_reads = []
	genome = fa2dict(args.genome) if args.genome else None
	bam_reader = pysam.AlignmentFile(args.bam)
	if args.genome:
		out_fa = open(args.out_prefix + '.fa', 'w')
	out_bam = pysam.AlignmentFile(args.out_prefix + '.bam', "wb", header=bam_reader.header)
	for read in bam_reader:
		if ref and read.reference_name != ref:
			if filt_map_reads(bam_reader.get_reference_length(ref), map_reads, out_bam, args.min_map_depth, args.flank) and genome:
				print(">%s\n%s" % (ref, genome[ref]), file=out_fa)
			ref = read.reference_name
			map_reads = []
		elif not ref:
			ref = read.reference_name
		rlen = read.infer_read_length()
		if read.is_unmapped or read.is_supplementary or rlen <= args.min_read_len or read.reference_end - read.reference_start <= args.min_map_len:
			continue
		ref_len = bam_reader.get_reference_length(read.reference_name)
		if read.query_alignment_end - read.query_alignment_start > rlen * (1 - args.max_clip_fra) or \
			(read.mapping_quality > args.min_map_qual and (is_head(read.query_alignment_start, read.query_alignment_end, rlen, read.reference_start, \
			read.reference_end, ref_len, args.max_clip_len) or is_tail(read.query_alignment_start, read.query_alignment_end, rlen, \
				read.reference_start, read.reference_end, ref_len, args.max_clip_len))):
			map_reads.append((read, get_mlen(read)))
	if map_reads:
		if filt_map_reads(bam_reader.get_reference_length(ref), map_reads, out_bam, args.min_map_depth, args.flank) and genome:
			print(">%s\n%s" % (ref, genome[ref]), file=out_fa)
	if args.genome:
		out_fa.close()
	out_bam.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
				formatter_class = HelpFormatter,
				description = '''
bam_region_depth:
	filter bam records using different thresholds

exmples:
	%(prog)s read.map2ref.sorted.bam 

'''
	)
	parser.add_argument("bam", type = str, help ="sorted .bam file, index file is not required.")
	parser.add_argument("out_prefix", type = str, help ="out results to PREFIX.fa, PREFIX.bam")
	parser.add_argument("--genome", metavar = 'FILE', type = str, help ="genome file in FASTA format.")
	parser.add_argument('--sr', action='store_true', default=False,
		help='short reads bam, long reads by default. use preset options: --min_read_len 100 --min_map_len 100 \
			--max_clip_len 10 --max_clip_fra 0.05')
	parser.add_argument('--min_map_depth', metavar = 'INT', type=int, default=-1,
		help='filter sequences with mapping depth <= INT')
	parser.add_argument('--min_read_len', metavar = 'INT', type=int, default=10000,
			help='filter reads with length <= INT.')
	parser.add_argument('--min_map_len', metavar = 'INT', type=int, default=5000,
			help='filter alignments with alignment length <= INT.')
	parser.add_argument('--min_map_qual', metavar = 'INT', type=int, default=1,
			help='filter alignments with mapping quality <= INT.')
	parser.add_argument('--max_clip_len', metavar = 'INT', type=int, default=100,
			help='filter alignments with unaligned length >= INT.')
	parser.add_argument('--max_clip_fra', metavar = 'FLOAT', type=float, default=0.1,
			help='filter alignments with unaligned length >= FLOAT x $read_length .')
	parser.add_argument('--flank', metavar = 'INT', type=int, default=10000,
		help='the flanking length of extracted sequences')
	args = parser.parse_args()
	if args.sr:
		if args.min_read_len == 10000:
			args.min_read_len = 100
		if args.min_map_len == 5000:
			args.min_map_len = 100
		if args.max_clip_len == 100:
			args.max_clip_len = 10
		if args.max_clip_fra == 0.1:
			args.max_clip_fra = 0.05
	main(args)
