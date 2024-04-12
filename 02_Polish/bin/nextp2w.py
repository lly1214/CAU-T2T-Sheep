#!/usr/bin/env python

import sys
import os
import gzip
import argparse
from Bio import SeqIO

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(SCRIPT_PATH + '/lib/')
from kit import *

log = ''
VERSION = 'v0.1.0'

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
	pass

def is_file(path):
	path = os.path.abspath(path)
	path_dir = os.path.dirname(path)
	if os.path.exists(path):
		if is_gzfile(path):
			return (path, [path])
		files = []
		with open (path, 'r') as IN:
			for line in IN:
				line = line.strip()
				if line.startswith('#') or not line:
					continue
				line_path = line if line.startswith('/') else path_dir + '/' + line
				if os.path.exists(line_path):
					files.append(line_path)
				elif files:
					raise FileNotFoundError(line)
				else:
					files.append(path)
					break
		return (path, files)
	else:
		raise FileNotFoundError(path)

def str_args(args):
	ret = ''
	for key, value in vars(args).items():
		if isinstance(value, tuple):
			ret += '%s: %s\n' % (key, value[0])
		else:
			ret += '%s: %s\n' % (key, str(value))
	return ret

def yak_cmd(args, out_file, ksize=21):
	cmd = SCRIPT_PATH + '/bin/yak count -k %d -o %s -t %d ' % (ksize, out_file, args.thread)
	if args.discard_singleton:
		cmd += ' -b31 '

	sr_files = args.sr[1]
	for i in range(len(sr_files) - 1):
		assert sr_files[i].split('.')[-1] == sr_files[i + 1].split('.')[-1], "short reads do not have the same type"
	cat = 'zcat' if is_gzfile(sr_files[0]) else 'cat'
	cmd += ' <(%s %s)' % (cat, ' '.join(sr_files))
	cmd += ' <(%s %s)\n' % (cat, ' '.join(sr_files)) # need twice for discard singletons
	return cmd

def qv_cmd(yak_file, genome, outfile, thread):
	cmd = SCRIPT_PATH + '/bin/yak qv -p -t %d %s %s > %s' % (thread, yak_file, genome, outfile)
	return cmd

def map_cmd(genome, infile, outfile, thread, map_opt):
	infile = " ".join(infile) if isinstance(infile, list) else infile
	return	SCRIPT_PATH + '/bin/minimap2 -ax %s -t %d %s %s|' % \
			(map_opt, thread, genome, infile) + \
		SCRIPT_PATH + '/bin/samtools sort -t %d -o %s\n' % \
			(thread, outfile)

def region_depth(bam, genome, out_prefix, **kwargs):
	return "python3 %s/bin/bam_region_depth.py %s %s %s %s" % \
		(SCRIPT_PATH, " ".join("--%s %s" % (k, v) for k, v in kwargs.items()), bam, genome, out_prefix)

def filt_bam(bam, out_prefix, **kwargs):
	return "python3 %s/bin/bam_filter.py %s %s %s" % \
		(SCRIPT_PATH, " ".join("--%s %s" % (k, v) for k, v in kwargs.items()), bam, out_prefix)

def falift(genome, infile, outfile, thread):
	bam = (outfile[:-5] if outfile.endswith('.fbed') else outfile) + '.bam'
	map_cmd = SCRIPT_PATH + '/bin/minimap2 -ax asm20 -Y -t %d %%s %%s 2>/dev/null' % thread
	cmd = "python3 %s/bin/map_one_by_one.py \"%s\" %s %s|%s/bin/samtools sort -t %d -o %s -\n" % \
		(SCRIPT_PATH, map_cmd, genome, infile, SCRIPT_PATH, thread, bam)
	cmd += "python3 %s/bin/bam2fbed.py %s > %s" % \
		(SCRIPT_PATH, bam, outfile)
	return cmd

def polish_hifi(yak_files, bam, genome, outfile, **kwargs):
	cmd = ''
	thread = kwargs.get("thread", 1)
	if os.path.exists(outfile):
		cmd += "rm -rf %s\n" % outfile
	if not os.path.exists(bam + '.bai'):
		cmd += SCRIPT_PATH + '/bin/samtools index -@ %d %s\n' % \
			(thread, bam)
	return cmd + "%s/bin/nextPolish2 --out %s %s %s %s %s\n" % \
		(SCRIPT_PATH, outfile, " ".join("--%s %s" % (k, v) for k, v in kwargs.items()), bam, genome, " ".join(yak_files))

def polish_ont_sr(args, genome, out_prefx):
	cmd = ''
	# polish with ONT
	bam = out_prefx + '.ont.bam'
	cmd = map_cmd(genome, args.ont[1], bam, args.thread, "map-ont")
	run(cmd, bam + '.done')

	out = out_prefx + '.ont.filt'
	cmd = filt_bam(bam, out)
	run(cmd, out + '.bam.done')

	bam = out + '.bam'
	bam_fofn = bam + '.fofn'
	out_genome = out_prefx + '.polishtemp1.fa'
	cmd = "\nls %s > %s\n" % (bam, bam_fofn)
	cmd += SCRIPT_PATH + '/bin/samtools index -@ %d %s\n' % (args.thread, bam)
	cmd += "python3 %s/bin/nextpolish2.py -sp -r ont -g %s -l %s -p %d -o %s\n" % \
		(SCRIPT_PATH, genome, bam_fofn, args.thread, out_genome)
	run(cmd, out_genome + '.done')

	# polish with short reads
	genome = out_genome
	cmd = '''
{path}/bin/minimap2 -ax sr -t {thread} {input} {reads} |{path}/bin/samtools sort -m 2g --threads {thread} - -o {out_prefx}.sgs.sort.bam
{bam_filter}
{path}/bin/samtools index -@ {thread} {out_prefx}.sgs.sort.filt.bam;
python3 {path}/bin/nextpolish1.py -g {input} -t 1 -p {thread} -s {out_prefx}.sgs.sort.filt.bam > {out_prefx}.polishtemp2.fa;
{path}/bin/minimap2 -ax sr -t {thread} {out_prefx}.polishtemp2.fa {reads} |{path}/bin/samtools sort -m 2g --threads {thread} - -o {out_prefx}.sgs.sort.bam
{bam_filter}
{path}/bin/samtools index -@ {thread} {out_prefx}.sgs.sort.filt.bam;
python3 {path}/bin/nextpolish1.py -g {out_prefx}.polishtemp2.fa -t 2 -p {thread} -s {out_prefx}.sgs.sort.filt.bam > {out_prefx}.nextpolish.fa;
	'''.format(path=SCRIPT_PATH, thread = args.thread, input = genome,
		reads = " ".join(args.sr[1]), out_prefx=out_prefx,
		bam_filter = filt_bam("%s.sgs.sort.bam" % out_prefx, "%s.sgs.sort.filt" % out_prefx, sr='', flank=300000))
	return cmd

def run(cmd, done_file):
	if os.path.exists(done_file):
		log.info("Skipped cmd: %s" % cmd)
		return True
	is_suc, out, err = run_cmd(cmd)
	if is_suc:
		open(done_file, 'w').close()
		log.info("Finished cmd: %s" % cmd)
		return out
	else:
		log.critical("Failed run %s, stdout:%s stderr:%s" % (cmd, out, err))

def fa2dict(path):
	seqs = {}
	handle = gzip.open(path, "rt") if is_gzfile(path) else open(path, 'r')
	for record in SeqIO.parse(handle, "fasta"):
		seqs[record.id] = str(record.seq).upper()
	return seqs

def qv2dict(infile):
	qvs = {}
	with open(infile) as IN:
		for line in IN:
			lines = line.strip().split()
			if lines[0] == "SQ":
				qvs[lines[1]] = float(lines[5])
	return qvs

def fbed2fa(path, outfile):
	OUT = open(outfile, 'w')
	with open(path) as IN:
		ref = None
		seq = ''
		for line in IN:
			lines = line.strip().split()
			if lines[0] != ref and seq:
				OUT.write(">%s\n%s\n" % (ref, seq))
				seq = ''
			ref = lines[0]
			seq += lines[1]
		if seq:
			OUT.write(">%s\n%s\n" % (ref, seq))
	OUT.close()

def generate_ont_gap_ref(hifi_gap_ref, hifi_gap_filt_ref, ont_gap_ref, out_file):
	
	def get_refid_prefix(refid):
		return "_".join(refid.split('_')[:2])

	hifi_gap_ref_dict = fa2dict(hifi_gap_ref)
	hifi_gap_filt_ref_dict = fa2dict(hifi_gap_filt_ref)
	ont_valid_refids = [get_refid_prefix(refid) for refid in hifi_gap_ref_dict if refid not in hifi_gap_filt_ref_dict]

	OUT = open(out_file, 'w')
	for refid, refseq in fa2dict(ont_gap_ref).items():
		if get_refid_prefix(refid) in ont_valid_refids:
			print (">%s\n%s" % (refid, refseq), file=OUT)

def read_nomap_bed(path):
	ret = {}
	with open(path) as IN:
		for line in IN:
			ref, st, en, lable = line.strip().split()
			if ref not in ret:
				ret[ref] = []
			ret[ref].append([int(st), int(en), lable == "True"])
	return ret

def extract_false_lable(bed):
	ret = {}
	for k, vs in bed.items():
		for s, e, l in vs:
			if not l:
				if k not in ret:
					ret[k] = []
				ret[k].append((s, e))
	return ret

def read_fbed_within_regions(path, regions, flank=100000):
	
	def is_in_regions(ref, pos):
		if ref not in regions:
			return False
		for s, e in regions[ref]:
			if s  - flank <= pos <= e + flank:
				return True
		return False

	ret = {}
	with open(path) as IN:
		for line in IN:
			ref, base, pos = line.strip().split()
			pos = int(pos)
			if is_in_regions(ref, pos):
				if ref not in ret:
					ret[ref] = {ref: []}
				ret[ref][ref].append((base, pos))
	return ret #input is ordered

def read_fbed(path):
	
	def parse_ref(ref):
		ref, st, en, lable, t_st, t_en = ref.split('_')
		return ref, int(t_st)

	ret = {}
	with open(path) as IN:
		for line in IN:
			lines = line.strip().split()
			ref, base, pos = lines[:3]
			pos = int(pos)
			ref, offset = parse_ref(ref)
			if ref not in ret:
				ret[ref] = {}
			if lines[0] not in ret[ref]:
				ret[ref][lines[0]] = []
			if len(lines) == 3:
				ret[ref][lines[0]].append((base, pos + offset))
			else:
				ret[ref][lines[0]].append((base, pos + offset, int(lines[3])))
	return ret

def merge_polish_genomes(hifi_polish_filt_ref, gap_bed, hifi_polish_nofilt_ref, hifi_polish_gap_fbed, hifi_polish_gap_fa, 
		ont_polish_gap_fbed, ont_polish_gap_fa, out_file):
	
	def iter_fbed(infile):
		seqs = []
		previos_ref = None
		with open(infile) as IN:
			for line in IN:
				lines = line.strip().split()
				ref, base, pos = lines[:3]
				pos = int(pos)
				if previos_ref and ref != previos_ref:
					yield (seqs, previos_ref)
					seqs = []
				previos_ref = ref
				seqs.append((base, pos))
		if seqs:
			yield (seqs, previos_ref)

	def in_ont_gap(ref, st, en):
		if ref not in ont_polish_gap_fbed:
			return False
		for name in ont_polish_gap_fbed[ref]:
			names = name.split('_')
			if names[0] == ref and int(names[1]) == st and int(names[2]) == en:
				return True
		return False

	def find_fbed_by_gap(seqs, st, en):
		if len(seqs) == 1:
			return list(seqs.items())[0]
		for name, fbed in seqs.items():
			names = name.split('_')
			if names[0] == ref and int(names[1]) == st and int(names[2]) == en:
				return name, fbed
		
		print("seqs doesn't contain regions: %d %d" % (st, en))
		sys.exit(1)

	def out_seqs(ref, seqs, OUT):
		OUT.write(">%s\n" % ref)
		for base in seqs:
			if len(base) == 2:
				base = base[0]
			OUT.write(base)
		OUT.write("\n")

	def find_index(seqs, pos):
		for index, seq in enumerate(seqs):
			if seq[1] == pos:
				return index
			elif seq[1] > pos:
				if index > 0 and seqs[index - 1][1] < pos:
					return index -1
				else:
					break
		return None

		print("seqs doesn't contain pos: %d" % pos)
		sys.exit(1)

	def matched_left_kmers(seqs1, seqs2, raw_pos, min_len=10):
		pos = raw_pos
		while True:
			seqs1_index = find_index(seqs1, pos)
			seqs2_index = find_index(seqs2, pos)
			if seqs1_index != None and seqs2_index != None:
				if raw_pos != pos:
					print("left shift pos from %d to %d " % (raw_pos, pos))
				break
			pos -= 1
			if pos < 0:
				return (-1, -1)

		while seqs1_index >= min_len and seqs2_index >= min_len:
			for i in range(min_len):
				seqs1_index -= i
				seqs2_index -= i
				base1, pos1 = seqs1[seqs1_index][:2]
				base2, pos2 = seqs2[seqs2_index][:2]
				if base1 == base2 and pos1 == pos2:
					continue
				else:
					if pos1 > pos2:
						seqs1_index -= 1
					elif pos1 < pos2:
						seqs2_index -= 1
					else:
						seqs1_index -= 1
						seqs2_index -= 1
					break
			else:
				print(2, raw_pos, pos, seqs1_index, seqs1[seqs1_index], seqs2_index, seqs2[seqs2_index])
				return (seqs1_index, seqs2_index)
		print("matched_left_kmers: seqs1 doesn't contain the same %d-kmers with seqs2, seqs1_index:%d seqs2_index:%d" % (min_len, seqs1_index, seqs2_index))
		return (-1, -1)

	def matched_right_kmers(seqs1, seqs2, raw_pos, min_len=10):
		pos = raw_pos
		while True:
			seqs1_index = find_index(seqs1, pos)
			seqs2_index = find_index(seqs2, pos)
			if seqs1_index != None and seqs2_index != None:
				if raw_pos != pos:
					print("right shift pos from %d to %d " % (raw_pos, pos))
				break
			pos += 1
			if pos >= seqs1[-1][1] or pos >= seqs2[-1][1]:
				return (-1, -1)

		while seqs1_index >= min_len and seqs2_index >= min_len:
			for i in range(min_len):
				seqs1_index += i
				seqs2_index += i
				base1, pos1 = seqs1[seqs1_index][:2]
				base2, pos2 = seqs2[seqs2_index][:2]
				if base1 == base2 and pos1 == pos2:
					continue
				else:
					if pos1 > pos2:
						seqs2_index += 1
					elif pos1 < pos2:
						seqs1_index += 1
					else:
						seqs1_index += 1
						seqs2_index += 1
					break
			else:
				return (seqs1_index, seqs2_index)
		print("matched_right_kmers: seqs1 doesn't contain the same %d-kmers with seqs2 seqs1_index:%d seqs2_index:%d" % (min_len, seqs1_index, seqs2_index))
		return (-1, -1)

	OUT = open(out_file, 'w')
	for seqs, ref in iter_fbed(hifi_polish_filt_ref):
		if ref in gap_bed:
			last_pos = 0
			new_seqs = ''
			for st, en, la in gap_bed[ref]:
				if la:
					if in_ont_gap(ref, st, en): #polished by ONT
						print("in_ont_gap", ref, st, en)
						gap_name, replace_fbed = find_fbed_by_gap(ont_polish_gap_fbed[ref], st, en)
						start1, start2 = (0, 0) if st == 0 else matched_left_kmers(seqs, replace_fbed, st)
						end1, end2 = (len(seqs) - 1, len(replace_fbed) - 1) if en == seqs[-1][1] else matched_right_kmers(seqs, replace_fbed, en)
						if start1 == -1 or end1 == -1:
							print("failed polish regions using ONT: %s %d %d" % (ref, st, en))
							continue

						start2 = replace_fbed[start2][2] # ont_polish_gap_fbed maybe contain some un-aligned regions 
						end2 = replace_fbed[end2][2]
						# ref2 = replace_fbed[end2][2] #ont_polish_gap_fa has a different ref name 
						if start1 != 0:
							assert start1 > last_pos, "in_ont_gap: gaps are overlapped: %d vs %d" % (last_pos, start1)
							new_seqs += "".join([base[0] for base in seqs[last_pos:start1]])
						print("in_ont_gap", ref, st, en, start1, start2, end1, end2)
						new_seqs += ont_polish_gap_fa[gap_name][start2:end2]
						last_pos = end1
					else:#polished by HiFi
						print("in_hifi_gap", ref, st, en)
						gap_name, replace_fbed = find_fbed_by_gap(hifi_polish_gap_fbed[ref], st, en)
						start1, start2 = (0, 0) if st == 0 else matched_left_kmers(seqs, replace_fbed, st)
						end1, end2 = (len(seqs) - 1, len(replace_fbed) - 1) if en == seqs[-1][1] else matched_right_kmers(seqs, replace_fbed, en)
						if start1 == -1 or end1 == -1:
							print("failed polish regions using HiFi: %s %d %d" % (ref, st, en))
							continue
						if start1 != 0:
							assert start1 > last_pos, "in_hifi_gap: gaps are overlapped: %d vs %d" % (last_pos, start1)
							new_seqs += "".join([base[0] for base in seqs[last_pos:start1]])
						print("in_hifi_gap", ref, st, en, start1, start2, end1, end2)
						new_seqs += hifi_polish_gap_fa[gap_name][start2:end2] # hifi_polish_gap_fbed maybe contain some un-aligned regions 
						last_pos = end1
				else:#polished by hifi_polish_nofilt_ref
					print("hifi_polish_nofilt_ref", ref, st, en)
					gap_name, replace_fbed = find_fbed_by_gap(hifi_polish_nofilt_ref[ref], st, en)
					start1, start2 = (0, 0) if st == 0 else matched_left_kmers(seqs, replace_fbed, st)
					end1, end2 = (len(seqs) - 1, len(replace_fbed) - 1) if en == seqs[-1][1] else matched_right_kmers(seqs, replace_fbed, en)
					if start1 != 0:
						assert start1 > last_pos, "hifi_polish_nofilt_ref: gaps are overlapped: %d vs %d" % (last_pos, start1)
						new_seqs += "".join([base[0] for base in seqs[last_pos:start1]])
					print("hifi_polish_nofilt_ref", ref, st, en, start1, start2, end1, end2)
					new_seqs += "".join([base[0] for base in replace_fbed[start2:end2]])
					last_pos = end1
			new_seqs +=  "".join([base[0] for base in seqs[last_pos:]])
			out_seqs(ref, new_seqs, OUT)
		else:
			out_seqs(ref, seqs, OUT)
	OUT.close()

def main(args):
	global log 
	args.log = 'pid' + str(os.getpid()) + '.' + args.log.strip('pidXXX.')
	log = plog(args.log)
	log.info('%s start...' % os.path.basename(sys.argv[0]))
	log.info('version:%s logfile:%s' % (VERSION, args.log))
	log.info('options: \n' + str_args(args))

	#generate yak files
	short_yak_file = args.workdir + '/01.sr.short.yak'
	cmd = yak_cmd(args, short_yak_file, 21)
	run(cmd, short_yak_file + '.done')
	long_yak_file = args.workdir + '/01.sr.long.yak'
	cmd = yak_cmd(args, long_yak_file, 31)
	run(cmd, long_yak_file + '.done')
	yak_files = [short_yak_file, long_yak_file]

	# run mapping HiFi-to-genome
	hifi_map_genome = args.workdir + '/02.hifi_map_to_genome.sort.bam'
	cmd = map_cmd(args.genome, args.hifi[1], hifi_map_genome, args.thread, "map-pb");
	run(cmd, hifi_map_genome + '.done')

	#run `bam_region_depth` on `hifi_map_genome` generated by last step
	infile = hifi_map_genome
	outfile = infile
	cmd = region_depth(infile, args.genome, outfile)
	run(cmd, outfile + '.nomap.done')
	nomap_bed = outfile + '.nomap.bed'
	nomap_bed_ref = outfile + '.nomap.false.fa'
	hifi_gap_ref = outfile + '.nomap.10000.fa'
	ont_gap_ref = outfile + '.nomap.300000.fa'

	#run mapping HiFi-to-`hifi_gap_ref`
	hifi_map_gap = args.workdir + '/03.hifi_map_to_gap.sort.bam'
	cmd = map_cmd(hifi_gap_ref, args.hifi[1], hifi_map_gap, args.thread, "map-pb");
	run(cmd, hifi_map_gap + '.done')

	# run `bam_filter` on `hifi_map_gap` generated by last step
	infile = hifi_map_gap
	hifi_map_gap = infile + '.filt'
	cmd = filt_bam(infile, hifi_map_gap, min_map_depth=3, genome=hifi_gap_ref)
	run(cmd, hifi_map_gap + '.done')

	#polishing gap regions using HiFi data with default options
	hifi_gap_filt_ref = hifi_map_gap + '.fa'
	hifi_gap_filt_bam = hifi_map_gap + '.bam'
	hifi_gap_filt_ref_polish = hifi_gap_filt_ref + ".hifipolish.filt.fa"
	cmd = polish_hifi(yak_files, hifi_gap_filt_bam, hifi_gap_filt_ref, hifi_gap_filt_ref_polish, thread=args.thread)
	run(cmd, hifi_gap_filt_ref_polish + '.done')

	hifi_polish_gap_fbed = hifi_gap_filt_ref_polish + '.fbed'
	cmd = falift(hifi_gap_filt_ref, hifi_gap_filt_ref_polish, hifi_polish_gap_fbed, args.thread)
	run(cmd, hifi_polish_gap_fbed + '.done')

	#generate gaps that do not span by HiFi
	out_file = args.workdir + '/04.hifi_map_to_genome.sort.bam.ont.ref.fa'
	out_file_done = out_file + '.done'
	if not os.path.exists(out_file_done):
		generate_ont_gap_ref(hifi_gap_ref, hifi_gap_filt_ref, ont_gap_ref, out_file)
		open(out_file_done, 'w').close()
	ont_gap_ref = out_file

	#polish `ont_gap_ref` using ONT reads and short reads
	cmd = polish_ont_sr(args, ont_gap_ref, ont_gap_ref)
	ont_gap_ref_polish = ont_gap_ref + '.nextpolish.fa'
	run(cmd, ont_gap_ref_polish + '.done')

	ont_gap_ref_fbed = ont_gap_ref_polish + '.fbed'
	cmd = falift(ont_gap_ref, ont_gap_ref_polish, ont_gap_ref_fbed, args.thread)
	run(cmd, ont_gap_ref_fbed + '.done')

	#polishing with `HiFi-to-genome` bam with default options
	hifi_map_genome_polish_filt = args.workdir + '/05.%s.hifipolish.filt.fbed' % os.path.basename(args.genome)
	cmd = polish_hifi(yak_files, hifi_map_genome, args.genome, hifi_map_genome_polish_filt, out_pos='', thread=args.thread)
	run(cmd, hifi_map_genome_polish_filt + '.done')

	#polishing regions with `False` lable in `nomap_bed`
	hifi_map_genome_polish_nofilt = args.workdir + '/06.%s.hifipolish.nofilt.fbed' % os.path.basename(args.genome)
	cmd = polish_hifi(yak_files, hifi_map_genome, nomap_bed_ref, hifi_map_genome_polish_nofilt, out_pos='', thread=args.thread, use_secondary='', min_map_qual=-1)
	run(cmd, hifi_map_genome_polish_nofilt + '.done')

	#merge polished sequences from different sources
	gap_bed = read_nomap_bed(nomap_bed)
	hifi_polish_nofilt_ref = read_fbed_within_regions(hifi_map_genome_polish_nofilt, regions=extract_false_lable(gap_bed))
	hifi_polish_gap_fbed = read_fbed(hifi_polish_gap_fbed)
	hifi_polish_gap_fa = fa2dict(hifi_gap_filt_ref_polish)
	ont_polish_gap_fbed = read_fbed(ont_gap_ref_fbed)
	ont_polish_gap_fa = dict((k[:-5] if k.endswith("_np12") else k, v) for k, v in fa2dict(ont_gap_ref_polish).items())
	hifi_polish_filt_ref = hifi_map_genome_polish_filt

	merge_polish_genome = args.workdir + '/07.%s.polish.merge.fasta' % os.path.basename(args.genome)
	merge_polish_genome_done_file = merge_polish_genome + '.done'
	if not os.path.exists(merge_polish_genome_done_file):
		merge_polish_genomes(hifi_polish_filt_ref, gap_bed, hifi_polish_nofilt_ref, hifi_polish_gap_fbed, hifi_polish_gap_fa, ont_polish_gap_fbed, ont_polish_gap_fa, merge_polish_genome)
		open(merge_polish_genome_done_file, 'w').close()

	#re-polish with HiFi data and output as the result
	hifi_map_bam = merge_polish_genome + '.sort.bam'
	cmd = map_cmd(merge_polish_genome, args.hifi[1], hifi_map_bam, args.thread, "map-hifi");
	run(cmd, hifi_map_bam + '.done')

	hifi_map_bam_polish = merge_polish_genome + '.polish.fasta'
	cmd = polish_hifi(yak_files, hifi_map_bam, merge_polish_genome, hifi_map_bam_polish, thread=args.thread)
	run(cmd, hifi_map_bam_polish + '.done')

	#caculate QV of raw & polished genomes 
	hifi_map_genome_polish_filt_fa = hifi_map_genome_polish_filt + '.fa'
	fbed2fa(hifi_map_genome_polish_filt, hifi_map_genome_polish_filt_fa)
	hifi_map_genome_polish_nofilt_fa = hifi_map_genome_polish_nofilt + '.fa'
	fbed2fa(hifi_map_genome_polish_nofilt, hifi_map_genome_polish_nofilt_fa)
	raw_genome = args.genome
	merge_polish_genome = merge_polish_genome
	hifi_map_bam_polish = hifi_map_bam_polish

	hifi_map_genome_polish_filt_fa_qv = hifi_map_genome_polish_filt_fa + '.qv'
	cmd = qv_cmd(short_yak_file, hifi_map_genome_polish_filt_fa, hifi_map_genome_polish_filt_fa_qv, args.thread)
	run(cmd, hifi_map_genome_polish_filt_fa_qv + '.done')

	hifi_map_genome_polish_nofilt_fa_qv = hifi_map_genome_polish_nofilt_fa + '.qv'
	cmd = qv_cmd(short_yak_file, hifi_map_genome_polish_nofilt_fa, hifi_map_genome_polish_nofilt_fa_qv, args.thread)
	run(cmd, hifi_map_genome_polish_nofilt_fa_qv + '.done')

	raw_genome_qv = args.workdir + '/08.%s.qv' % os.path.basename(raw_genome)
	cmd = qv_cmd(short_yak_file, raw_genome, raw_genome_qv, args.thread)
	run(cmd, raw_genome_qv + '.done')

	merge_polish_genome_qv = merge_polish_genome + '.qv'
	cmd = qv_cmd(short_yak_file, merge_polish_genome, merge_polish_genome_qv, args.thread)
	run(cmd, merge_polish_genome_qv + '.done')

	hifi_map_bam_polish_qv = hifi_map_bam_polish + '.qv'
	cmd = qv_cmd(short_yak_file, hifi_map_bam_polish, hifi_map_bam_polish_qv, args.thread)
	run(cmd, hifi_map_bam_polish_qv + '.done')

	#extract each seq with the largest QV from different sources as the finally result
	finally_polish_genome = args.workdir + '/08.%s.polish.fasta' % os.path.basename(args.genome)
	best_qvs = qv2dict(raw_genome_qv)
	out_genomes = fa2dict(raw_genome)
	for fa_file, qv_file in zip([hifi_map_bam_polish, merge_polish_genome, hifi_map_genome_polish_filt_fa, hifi_map_genome_polish_nofilt_fa], 
		[hifi_map_bam_polish_qv, merge_polish_genome_qv, hifi_map_genome_polish_filt_fa_qv, hifi_map_genome_polish_nofilt_fa_qv]):
		seqs = fa2dict(fa_file)
		for ref, qv in qv2dict(qv_file).items():
			if best_qvs[ref] < qv:
				best_qvs[ref] = qv
				out_genomes[ref] = seqs[ref]
	with open(finally_polish_genome, 'w') as OUT:
		for k, v in out_genomes.items():
			OUT.write(">%s\n%s\n" % (k, v))
	finally_polish_genome_qv = finally_polish_genome + '.qv'
	cmd = qv_cmd(short_yak_file, finally_polish_genome, finally_polish_genome_qv, args.thread)
	run(cmd, finally_polish_genome_qv + '.done')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		formatter_class = HelpFormatter,
		description = '''
NextPolish2_Wrap:
	A wrapper for NextPolish2, used to polish the genome assembly by HiFi and ONT data

exmples: 
	%(prog)s --sr sr.fofn --ont ont.fofn --hifi hifi.fofn --genome genome.fa --thread 30

'''
	)
	VERSION = getver(VERSION, SCRIPT_PATH)
	parser.add_argument('-v', '--version', action='version', version=VERSION)
	parser.add_argument ('--log', metavar = 'FILE',type = str, default = 'pidXXX.log.info',
		help = 'log file')
	parser.add_argument('--sr', metavar = 'FILE', type=is_file, required=True, default=argparse.SUPPRESS,
		help='short reads file in [GZIP] FASTA/FASTQ/FOFN format.')
	parser.add_argument('--ont', metavar = 'FILE', type=is_file, required=True, default=argparse.SUPPRESS,
		help='ONT reads file in [GZIP] FASTA/FASTQ/FOFN format.')
	parser.add_argument('--hifi', metavar = 'FILE', type=is_file, required=True, default=argparse.SUPPRESS,
		help='HiFi reads file in [GZIP] FASTA/FASTQ/FOFN format.')
	parser.add_argument('--genome', metavar = 'FILE', type=is_file, required=True, default=argparse.SUPPRESS,
		help='genome file in [GZIP] FASTA format.')
	parser.add_argument('--discard_singleton', action='store_true', default=False, 
		help='discard singletons when counting kmers from short reads')
	parser.add_argument('--workdir', metavar = 'DIR', type=str, default='INPUT_GENOME_FILE_NAME.nextp2w.work',
		help='work directory.')
	parser.add_argument('--thread', metavar = 'INT', type=int, default=10,
		help='number of threads.')
	args = parser.parse_args()
	if args.workdir == 'INPUT_GENOME_FILE_NAME.nextp2w.work':
		args.workdir = os.path.abspath(os.path.basename(args.genome[0]) + '.nextp2w.work')
	if not os.path.exists(args.workdir):
		os.mkdir(args.workdir)
	args.genome = args.genome[0]
	main(args)
