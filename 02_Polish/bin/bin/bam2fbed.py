import sys
import gzip
import pysam

def remove_np(name):
	return name[:-5] if name.endswith("_np12") else name

def out_result(tname, result):
	for i in range(len(result)):
		qp, qb = result[i]
		if qb != 'N' and qb != '-':
			print ("%s\t%c\t%d\t%d" % (tname, qb, i, qp))

def read_bams(infile):
	bam = pysam.AlignmentFile(bamfile)
	result = {}
	for read in bam:
		read.qname = remove_np(read.qname)
		if read.is_unmapped or read.is_secondary or read.qname != read.reference_name:
			continue
		if read.reference_name not in result:
			result[read.reference_name] = []
		qstart, qend = read.query_alignment_start, read.query_alignment_end
		if read.is_reverse:
			qstart = read.query_length - read.query_alignment_end
			qend = read.query_length - read.query_alignment_start
		result[read.reference_name].append([qstart, qend , read.reference_start, read.reference_end, read.cigartuples, read.is_reverse, read.query_sequence, bam.get_reference_length(read.reference_name)])
	for refname in result:
		result[refname] = sorted(result[refname], key = lambda x: (x[2]-x[3]))
	return result

def filt_alns(bam_data):
	for refname in bam_data:
		result = []
		for i, aln in enumerate(bam_data[refname]):
			qs, qe, rs, re, cigar, rev, seq, rlen = aln
			is_valid = True
			for j, (qs_t, qe_t, rs_t, re_t, _, _, _, _) in enumerate(bam_data[refname]):
				if j >= i:
					break
				if qs_t < qs < qe < qe_t and rs_t < rs < re < re_t:
					is_valid = False
					break
			if is_valid:
				result.append(aln)
		bam_data[refname] = result

def filt_reverse(bam_data):

	def not_overlap_with_rev_alns(qs, qe, rs, re, rev_alns):
		for qs_t, qe_t, rs_t, re_t in rev_alns:
			if qs_t <= qs <= qe_t or qs_t <= qe <= qe_t or rs_t <= rs <= re_t or rs_t <= re <= re_t:
				return False
		return True

	for refname in bam_data:
		result = []
		rev_alns = []
		for qs, qe, rs, re, cigar, rev, seq, rlen in bam_data[refname]:
			if rev:
				rev_alns.append([qs, qe, rs, re])
			elif not_overlap_with_rev_alns(qs, qe, rs, re, rev_alns): 
				result.append([qs, qe, rs, re, cigar, rev, seq, rlen])
			elif re - rs > 30000 and abs(qs - rs) < 10000 and abs(re - qe) < 10000:
				result.append([qs, qe, rs, re, cigar, rev, seq, rlen])
			 
		bam_data[refname] = result
		
#out file format, input: asm_map_to_ref
# refid, asm_base, ref_pos, asm_pos

if len(sys.argv) != 2:
	print ("python %s asm.sort.bam > result.fbed" % sys.argv[0])
	sys.exit(1)

bamfile = sys.argv[1]
bam_data = read_bams(bamfile)
filt_alns(bam_data)
filt_reverse(bam_data)

# for i in bam_data:
# 	for l in bam_data[i]:
# 		print (i, l[:4], l[5])

for qname, alns in bam_data.items():
	result = [[0, 'N']] * (alns[0][7])
	for qs, qe, rs, re, cigar, rev, seq, rlen in alns:
		assert not rev, "%s has reversed alns" % qname
		for (c, l) in cigar:
			if c == 0:
				for i in range(l):
					tp = rs + i
					qp = qs + i
					if result[tp][1] != 'N': #duplication mapping regions
						result[tp][1] = '-'
					else:
						result[tp] = [qp, seq[qp]]
				qs += l
				rs += l
			elif c == 1:
				qs += l
			elif c == 2:
				rs += l
			elif c == 4:
				continue
			elif c == 5:
				print("minimap2 need to be run with `-Y`", file=sys.stderr)
				sys.exit(1)
			else:
				print("%s contains a unknown cigar: %d %d" % (qname, c, l), file=sys.stderr)
				sys.exit(1)
		assert qs == qe and rs == re, (qname, qs, qe, rs, re)

	out_result(qname, result)
