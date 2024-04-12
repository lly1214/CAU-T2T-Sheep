import operator
import sys
import os
import argparse
import copy
from collections import defaultdict
import pysam

def read_regions(in_low_cov_lst):
	low_cov_regions = dict()
	if os.path.exists(in_low_cov_lst):
		with read_file(in_low_cov_lst) as fin:
			for line in fin:
				arr1 = line.strip().split()
				if len(arr1) < 2:
					continue
				arr2 = arr1[1].split(';')
				regions = list()
				for x in arr2:
					s, e = x.split(',')
					regions.append((int(s), int(e)))
				regions.sort(key=operator.itemgetter(0, 1))
				low_cov_regions[arr1[0]] = regions
	return low_cov_regions  # []


def read_aligns(in_align_lst):
	if not os.path.exists(in_align_lst):
		return
	align_list = set()
	with read_file(in_align_lst) as fin:
		for line in fin:
			arr1 = line.strip().split()
			# rid, chr_id, ts, strand
			align_list.add((arr1[0], arr1[1], int(arr1[2]), arr1[3]))
	return align_list


def adjust_cigar_start(rstart, r_new_start, cigar_tuples):
	new_cigar = list()
	copy_flag = False
	q_new_start = -1
	r_coor, q_coor = rstart, 0
	for tag, length in cigar_tuples:
		if copy_flag:
			new_cigar.append((tag, length))
			continue
		if tag == 0 or tag == 7 or tag == 8:  # M, =, X
			if r_coor <= r_new_start <= r_coor + length - 1:
				q_new_start = q_coor + r_new_start - r_coor
				new_cigar.append((tag, r_coor + length - r_new_start))
				copy_flag = True
				continue
			r_coor += length
			q_coor += length
		elif tag == 1 or tag == 4:  # I, S
			q_coor += length
		elif tag == 2:  # D
			if r_coor <= r_new_start <= r_coor + length - 1:
				q_new_start = q_coor
				new_cigar.append((tag, r_coor + length - r_new_start))
				copy_flag = True
				continue
			r_coor += length
	return q_new_start, new_cigar


def adjust_cigar_end(rstart, r_new_start, cigar_tuples):
	new_cigar = list()
	q_new_start = -1                         # query
	r_coor, q_coor = rstart, 0               # ref, query
	for tag, length in cigar_tuples:
		if tag == 0 or tag == 7 or tag == 8:  # M, =, X
			if r_coor <= r_new_start <= r_coor + length - 1:
				q_new_start = q_coor + r_new_start - r_coor
				new_cigar.append((tag, r_new_start - r_coor + 1))
				break
			r_coor += length
			q_coor += length
		elif tag == 1 or tag == 4:  # I, S
			q_coor += length
		elif tag == 2:  # D
			if r_coor <= r_new_start <= r_coor + length - 1:
				q_new_start = q_coor - 1
				new_cigar.append((tag, r_new_start - r_coor + 1))
				break
			r_coor += length
		new_cigar.append((tag, length))
	return q_new_start, new_cigar


def get_ql_tl(cigartuples):
	r_coor = q_coor = 0
	for tag, length in cigartuples:
		if tag == 0 or tag == 7 or tag == 8:  # M, =, X
			r_coor += length
			q_coor += length
		elif tag == 1 or tag == 4:  # I, S
			q_coor += length
		elif tag == 2:  # D
			r_coor += length
	return r_coor, q_coor


def is_ter(ts, start_coors, fuzz=5000):
	for x in start_coors:
		if abs(ts - x) < fuzz:
			return True
	return False


def read_coors(start_coor_lst):
	start_coors = defaultdict(list)
	with read_file(start_coor_lst) as fin:
		for line in fin:
			arr1 = line.strip().split("\t")
			arr2 = arr1[1].split(',')
			for x in arr2:
				if x.strip() != '':
					start_coors[arr1[0]].append(int(x))
	for vals in start_coors.values():
		vals2 = list(set(vals))
		vals2.sort()
		vals[:] = vals2[:]
	return start_coors


def read_haplotypes(hap_read_lst_file):
	if not os.path.exists(hap_read_lst_file):
		return
	hap_read_lst = set()
	with read_file(hap_read_lst_file) as fin:
		for line in fin:
			if line.startswith('@'):
				continue
			arr1 = line.strip().split()
			hap_read_lst.add(arr1[0])
	return hap_read_lst


def main(args):
	start_coors = read_coors(args.start_coor_lst)
	end_coors = read_coors(args.end_coor_lst)

	out_dir = os.path.dirname(os.path.abspath(args.out_bam))
	os.system('mkdir -p %s' % out_dir)

	align_list = read_aligns(args.in_align_lst)
	regions = read_regions(args.in_region_lst)

	hap_read_lst = read_haplotypes(args.hap_read_lst)

	if os.path.exists(args.add_seq_faidx):
		query_ids = set()
		target_ids = set()
		f1 = pysam.AlignmentFile(args.in_bam, 'rb')
		for line in f1.fetch():
			if line.cigartuples[0][0] == 5 or line.cigartuples[-1][0] == 5:
				continue
			target_ids.add(line.reference_name)
			if line.query_sequence is not None:
				continue
			query_ids.add(line.query_name)
		f1.close()
		sel_ids = NdIo.get_selected_info_for_faidx(args.add_seq_faidx, query_ids, is_list=args.is_list)
		sel_seqs = NdIo.read_faidx(sel_ids, ids=None, get_fragment=False)
	else:
		sel_ids = None
		sel_seqs = None
		target_ids = None

	# cigar: M, 0; I, 1; D, 2; N, 3; S, 4; H, 5; P, 6; =, 7; X, 8; B, 9;
	kept_rgs = set([x for x in args.kept_rg.split(',') if len(x.strip()) > 0])
	f1 = pysam.AlignmentFile(args.in_bam, 'rb')
	f1_h = f1.header.to_dict()
	if 'RG' in f1_h:
		new_rg = [x for x in f1_h['RG'] if x['ID'] in kept_rgs]
		f1_h['RG'] = new_rg
	if target_ids is not None:
		new_sq = [x for x in f1_h['SQ'] if x['SN'] in target_ids]
		f1_h['SQ'] = new_sq
	ref_name_to_id = dict()
	for i, sq_dict in enumerate(f1_h['SQ']):
		ref_name_to_id[sq_dict['SN']] = i
	f2 = pysam.AlignmentFile(args.out_bam, "wb", header=f1_h)
	for line in f1.fetch():
		assert(isinstance(line, pysam.AlignedSegment))
		if hap_read_lst is not None and line.query_name not in hap_read_lst:
			continue
		ql = line.infer_read_length()
		tmp_qs, tmp_qe = line.query_alignment_start, line.query_alignment_end - 1
		ts, te = line.reference_start, line.reference_start + line.reference_length - 1  # []
		if line.cigartuples[0][0] == 5:
			tmp_qs += line.cigartuples[0][1]
			tmp_qe += line.cigartuples[0][1]
		if line.is_reverse:
			qe = ql - 1 - tmp_qs
			qs = ql - 1 - tmp_qe  # []
		else:
			qs = tmp_qs
			qe = tmp_qe  # []
		is_start_ter = is_ter(ts, start_coors[line.reference_name])
		is_end_ter = is_ter(te, end_coors[line.reference_name])
		if not is_start_ter and not is_end_ter and (qe - qs + 1) / ql < 0.9:
			continue
		if line.is_reverse:
			left_q_flank = ql - qe - 1
			right_q_flank = qs
		else:
			left_q_flank = qs
			right_q_flank = ql - qe - 1
		if (not is_start_ter and left_q_flank >= 500) or (not is_end_ter and right_q_flank >= 500):
			continue
		strand = "-" if line.is_reverse else "+"
		if args.change_flag:
			if line.flag & 0x100:  # secondary
				line.flag ^= 0x100
			if line.flag & 0x800:  # supplementary
				line.flag ^= 0x800
			line.mapping_quality = 60
			line.set_tag('tp', 'P')
		if align_list is not None and (line.query_name, line.reference_name, line.reference_start, strand) not in align_list:
			continue
		if len(regions) == 0:
			f2.write(line)
			continue
		if line.reference_name not in regions:
			continue
		if line.query_sequence is None and sel_ids is not None:
			# if line.query_name == '383373':
			# 	import pdb; pdb.set_trace()
			if line.cigartuples[0][0] == 5 or line.cigartuples[-1][0] == 5:
				continue
			if line.is_reverse:
				line.query_sequence = NdIo.reverse_seq(sel_seqs[line.query_name])
			else:
				line.query_sequence = sel_seqs[line.query_name]
		orig_cigar, orig_qseq, orig_qqual = line.cigartuples, line.query_sequence, line.query_qualities
		for s, e in regions[line.reference_name]:
			if ts <= s <= te:
				q_new_start, line.cigartuples = adjust_cigar_start(ts, s, orig_cigar)
				if q_new_start != -1:
					line.reference_start = s
					line.query_sequence = orig_qseq[q_new_start:]
					if orig_qqual is not None:
						line.query_qualities = orig_qqual[q_new_start:]
					# line.query_alignment_start = q_new_start
					line.reference_id = ref_name_to_id[line.reference_name]
					f2.write(line)
			elif ts <= e <= te:
				q_new_end, line.cigartuples = adjust_cigar_end(ts, e, orig_cigar)
				if q_new_end != -1:
					line.query_sequence = orig_qseq[:q_new_end + 1]
					if orig_qqual is not None:
						line.query_qualities = orig_qqual[:q_new_end + 1]
					# line.query_alignment_end = q_new_end + 1
					line.reference_id = ref_name_to_id[line.reference_name]
					f2.write(line)
				# if line.query_name == '2899399':
					# print(q_new_end, ts, e, strand, len(line.query_sequence), len(orig_qseq), *get_ql_tl(line.cigartuples), *get_ql_tl(orig_cigar))
					#print(line.cigartuples)
			elif s <= ts <= e:
				line.cigartuples, line.query_sequence, line.query_qualities = orig_cigar, orig_qseq, orig_qqual
				line.reference_id = ref_name_to_id[line.reference_name]
				f2.write(line)
				break
	f1.close()
	f2.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('in_bam', help='input bam file')
	parser.add_argument('out_bam', help='output bam file')
	parser.add_argument('-al', dest='in_align_lst', type=str, metavar='FILE', default='null.lst',
		help='input alignment list file, default=null.lst')
	parser.add_argument('-cf', dest='change_flag', action='store_true', default=False,
		help='change flags and quality or not, default=False')
	parser.add_argument('-rl', dest='in_region_lst', type=str, metavar='FILE', default='null.lst',
		help='input region list file, default=null.lst')
	parser.add_argument('-sl', dest='start_coor_lst', type=str, metavar='FILE', default='',
		help='contig start coor list, default=""')
	parser.add_argument('-el', dest='end_coor_lst', type=str, metavar='FILE', default='',
		help='contig end coor list, default=""')
	parser.add_argument('-as', dest='add_seq_faidx', type=str, metavar='FILE', default='',
		help='add sequence from faidx file, default=""')
	parser.add_argument('-il', dest='is_list', action='store_true', default=False,
		help='faidx file is list or not, default=False')
	parser.add_argument('-hl', dest='hap_read_lst', type=str, metavar='FILE', default='',
		help='filter reads by haplotype read list, default=""')
	parser.add_argument('-rg', dest='kept_rg', type=str, metavar='STR', default='M1,P1,unknown',
		help='kept read groups, default="M1,P1,unknown"')
	args = parser.parse_args()
	main(args)
