#!/usr/bin/env python
import sys
import pysam
import argparse
import gzip
from Bio import SeqIO

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
	pass

class Region(object):
	def __init__(self, st, en, mapq, is_secondary):
		self.st = st
		self.en = en
		self.mapq = mapq
		self.is_secondary = is_secondary

def out_map_regions(reflength, regions, mapq, ont_flank):

	def fill_regions(regions, s, e):
		if len(regions) == 0 or s > regions[-1][1]:
			regions.append([s, e])
		elif e > regions[-1][1]:
			regions[-1][1] = e

	def get_unmap_regions(regions, flank):
		unmap_regions = []
		for i in range(len(regions)):
			s, e = regions[i]
			if i == 0:
				if s > 100 :
					unmap_regions.append([0, s - 1])
			else:
				unmap_regions.append([regions[i-1][1], s - 1])
			if i == len(regions) - 1:
				if e < reflength - 100:
					unmap_regions.append([e, reflength - 1])

		merge_unmap_regions = []
		ns = ne = 0
		for s, e in unmap_regions:
			if ns == ne:
				ns, ne = s, e
			elif s - ne < flank:
				ne = e
			else:
				merge_unmap_regions.append((ns, ne))
				ns, ne = s, e
		if ns != ne:
			merge_unmap_regions.append((ns, ne))
		return merge_unmap_regions

	def contain_regions(regions, s, e):
		for rs, re in regions:
			if rs <= s < e <= re:
				return True
		return False

	unmap_regions = []
	regions.sort(key = lambda x: (x.st, x.en))
	regions_no_filt = []
	regions_filt_mapq_secondary = []

	for r in regions:
		if r.mapq > mapq and not r.is_secondary:
			fill_regions(regions_filt_mapq_secondary, r.st, r.en)
		fill_regions(regions_no_filt, r.st, r.en)
	
	for s, e in get_unmap_regions(regions_filt_mapq_secondary, flank=ont_flank):
		if contain_regions(regions_no_filt, s, e):
			unmap_regions.append([s, e, False])
		else:
			unmap_regions.append([s, e, True])
	return unmap_regions

def main(args):
	ref = None
	unmap_regions = {}
	map_regions = []
	bam_reader = pysam.AlignmentFile(args.bam)
	for read in bam_reader:
		if ref and read.reference_name != ref:
			unmap_regions[ref] = out_map_regions(bam_reader.get_reference_length(ref), map_regions, args.min_map_qual, max(args.flank))
			ref = read.reference_name
			map_regions = []
		elif not ref:
			ref = read.reference_name
		rlen = read.infer_read_length()
		if read.is_unmapped or read.is_supplementary or rlen <= args.min_read_len or \
			(read.reference_end - read.reference_start) <= args.min_map_len or \
			rlen - (read.query_alignment_end - read.query_alignment_start) >= args.max_clip_len:
			continue
		map_regions.append(Region(read.reference_start, read.reference_end, \
			read.mapping_quality, read.is_secondary))
	if map_regions:
		unmap_regions[ref] = out_map_regions(bam_reader.get_reference_length(ref), map_regions, args.min_map_qual)

	with open (args.out_prefix + '.nomap.bed', 'w') as OUT:
		for ref in unmap_regions:
			for s, e, l in unmap_regions[ref]:
				print ("%s\t%d\t%d\t%s" % (ref, s, e, l), file=OUT)

	for flank in args.flank:
		with open (args.out_prefix + '.nomap.%d.fa' % flank, 'w') as OUT:
			handle = gzip.open(args.genome, "rt") if args.genome.endswith(('gz', 'gzip', 'GZ', 'GZIP')) else open(args.genome, 'r')
			for record in SeqIO.parse(handle, "fasta"):
				ref = record.id
				seq = str(record.seq)
				seq_len = len(seq)
				for s, e, l in unmap_regions[ref]:
					if l:
						true_s = s - flank if s > flank else 0
						true_e = e + flank if e + flank < seq_len else seq_len
						print (">%s_%d_%d_true_%d_%d\n%s" % (ref, s, e, true_s, true_e, seq[true_s:true_e + 1]), file=OUT)

	with open (args.out_prefix + '.nomap.false.fa', 'w') as OUT:
		handle = gzip.open(args.genome, "rt") if args.genome.endswith(('gz', 'gzip', 'GZ', 'GZIP')) else open(args.genome, 'r')
		for record in SeqIO.parse(handle, "fasta"):
			ref = record.id
			seq = str(record.seq)
			seq_len = len(seq)
			out = False
			for s, e, l in unmap_regions[ref]:
				if not l:
					out = True
					break
			if out:
				print (">%s\n%s" % (ref, seq), file=OUT)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
				formatter_class = HelpFormatter,
				description = '''
bam_region_depth:
	compute mapping depth using different thresholds

exmples:
	%(prog)s read.map2ref.sorted.bam genome.fasta out_prefix

'''
	)
	parser.add_argument("bam", type = str, help ="sorted .bam file, index file is not required.")
	parser.add_argument("genome", type = str, help ="genome assembly file in [GZIP] FASTA format.")
	parser.add_argument("out_prefix", type = str, help ="out result to PREFIX.nomap.flank.fa, PREFIX.nomap.bed")
	parser.add_argument('--min_read_len', metavar = 'INT', type=int, default=1000,
			help='filter reads with length <= INT.')
	parser.add_argument('--min_map_len', metavar = 'INT', type=int, default=500,
			help='filter alignments with alignment length <= INT.')
	parser.add_argument('--min_map_qual', metavar = 'INT', type=int, default=1,
			help='filter alignments with mapping quality <= INT.')
	parser.add_argument('--max_clip_len', metavar = 'INT', type=int, default=100,
			help='filter alignments with unaligned length >= INT.')
	parser.add_argument('--flank', metavar = 'INT', type=int, default=[10000, 300000], action='append', nargs='+',
			help='The flanking length of extracted sequences')
	args = parser.parse_args()
	main(args)
