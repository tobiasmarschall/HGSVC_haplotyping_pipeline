#!/usr/bin/env python3
import sys
import datetime
from collections import defaultdict, deque
import vcf
from whatshap.args import HelpfulArgumentParser as ArgumentParser
#from whatshap.vcf import VcfReader
import pyfaidx

vcf_header = '''##fileformat=VCFv4.1
##fileDate={date}
##source=make-alt-explicit.py
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}'''

def info_to_str(info):
	keys = sorted(info.keys())
	l = []
	for key in keys:
		value = info[key]
		if value is None:
			l.append(key)
		elif type(value) == list:
			value_str = ','.join( (str(x) if x is not None else '.') for x in value)
			l.append('{}={}'.format(key, value_str))
		else:
			l.append('{}={}'.format(key, value))
	return ';'.join(l)

def nonedot(x):
	if x is None:
		return '.'
	elif type(x) == list:
		if len(x) == 0:
			return '.'
		else:
			return ','.join(x)
	else:
		return str(x)


def sv_start_end(record):
	'''Return start and end positions in pythonic coordinates'''
	assert 'SVLEN' in record.INFO
	assert 'SVTYPE' in record.INFO
	if record.INFO['SVTYPE'] == 'DEL':
		return record.POS, record.POS + abs(record.INFO['SVLEN'])
	elif record.INFO['SVTYPE'] == 'INS':
		return record.POS, record.POS + 1
	else:
		assert False, 'Unexpected SVTYPE: {}'.format(record.INFO['SVTYPE'])


def sv_distance(record1, record2):
	if record1.CHROM != record2.CHROM:
		return float('inf')
	else:
		if record1.POS > record2.POS:
			record1, record2 = record2, record1
		start1, end1 = sv_start_end(record1)
		start2, end2 = sv_start_end(record2)
		return max(0, start2 - end1)


def vcf_distance_iterator(vcf_reader):
	records = deque()
	written = deque()
	for record in vcf_reader:
		records.append(record)
		written.append(False)
		if len(records) == 3:
			assert not written[1]
			yield records[1], min(sv_distance(records[0],records[1]), sv_distance(records[1],records[2]))
			written[1] = True
			records.popleft()
			written.popleft()
	assert len(records) <= 2
	if len(records) == 2:
		d = sv_distance(records[0],records[1])
		for w, r in zip(written, records):
			if not w:
				yield r, d
	elif len(records) == 1:
		yield records[0], float('inf')


def main():
	parser = ArgumentParser(prog='make-alt-explicit.py', description=__doc__)
	parser.add_argument('ref', metavar='REF', help='Reference genome (FASTA)')
	parser.add_argument('vcf', metavar='VCF', help='VCF of merged PacBio variants')
	parser.add_argument('--min-distance', default=0, type=int, help='Minimum distance to next variant.')
	parser.add_argument('--max-length', default=None, type=int, help='Maxmimum length.')
	args = parser.parse_args()
	vcf_reader = vcf.Reader(filename=args.vcf)

	with pyfaidx.Fasta(args.ref, as_raw=True, sequence_always_upper=True) as fasta:
		reference = {}
		for chromosome in list(fasta.records):
			reference[chromosome] = str(fasta[chromosome])
			if len(reference[chromosome]) > 1000000:
				print('Loaded chromosome', chromosome, file=sys.stderr)
	
	#f = open('chr22.fa', 'w')
	#print('>chr22', file=f)
	#print(reference['chr22'], file=f)

	assert len(vcf_reader.samples) == 1
	print(vcf_header.format(date=datetime.datetime.now().strftime('%Y%m%d'), sample=vcf_reader.samples[0]))

	insertions = 0
	deletions = 0
	for record, distance in vcf_distance_iterator(vcf_reader):
		print('Processing', record.CHROM, record.POS, record.REF, record.ALT, 'length={}, distance={}'.format(record.INFO['SVLEN'],distance), file=sys.stderr)
		if distance < args.min_distance:
			print('  -- skipping: too close to next SV', file=sys.stderr)
			continue
		if args.max_length and (abs(record.INFO['SVLEN']) > args.max_length):
			print('  -- skipping: too long', file=sys.stderr)
			continue
		assert len(record.ALT) == 1
		assert len(record.REF) == 1
		assert record.REF[0] == reference[record.CHROM][record.POS - 1]
		assert 'SEQ' in record.INFO
		assert 'SVLEN' in record.INFO
		assert len(record.INFO['SEQ']) == record.INFO['SVLEN'], '{}:{}: len(SEQ)={}, SVLEN={}, SEQ={}'.format(record.CHROM,record.POS,len(record.INFO['SEQ']), record.INFO['SVLEN'], record.INFO['SEQ'])
		if str(record.ALT[0]) == '<DEL>':
			print(
				nonedot(record.CHROM),
				nonedot(record.POS),
				nonedot(record.ID),
				reference[record.CHROM][record.POS-1:record.POS+abs(record.INFO['SVLEN'])].upper(),
				nonedot(record.REF),
				nonedot(record.QUAL),
				nonedot(record.FILTER),
				info_to_str(record.INFO),
				'GT',
				record.samples[0].data.GT,
			sep='\t')
		elif str(record.ALT[0]) == '<INS>':
			#assert len(record.INFO['SEQ']) == record.INFO['SVLEN'], '{}:{}: len(SEQ)={}, SVLEN={}'.format(record.CHROM,record.POS,len(record.INFO['SEQ']), record.INFO['SVLEN'])
			#assert record.REF == record.INFO['SEQ'][0].upper(), '{}:{}: REF={}, SEQ={}'.format(record.CHROM,record.POS,record.REF,record.INFO['SEQ'])
			print(
				nonedot(record.CHROM),
				nonedot(record.POS),
				nonedot(record.ID),
				nonedot(record.REF),
				record.REF + record.INFO['SEQ'].upper(),
				nonedot(record.QUAL),
				nonedot(record.FILTER),
				info_to_str(record.INFO),
				'GT',
				record.samples[0].data.GT,
			sep='\t')
		else:
			assert False, 'Unexpected ALT field: {}'.format(record.ALT[0])
		
	#samples_freebayes = set(vcf_reader_freebayes.samples)
	#samples_10x = set()
	#for vcf_reader in vcf_readers_10x:
		#assert len(set(vcf_reader.samples).intersection(samples_10x)) == 0
		#samples_10x.update(vcf_reader.samples)
	#assert samples_freebayes == samples_10x
	#samples = list(sorted(samples_freebayes))
	#print('Samples:', ', '.join(samples), file=sys.stderr)

	## list of dicts mapping chromosome names to VariantTables for each input VCF
	#variants_10x = []
	#for vcf_reader, filename in zip(vcf_readers_10x, args.vcf10x):
		#m = {}
		#for variant_table in vcf_reader:
			#m[variant_table.chromosome] = variant_table
		#variants_10x.append(m)
		#print('Successfully read', filename, file=sys.stderr)

	#print(vcf_header.format(date=datetime.datetime.now().strftime('%Y%m%d'), sample='\t'.join(samples)))

	#concordance_table = defaultdict(int)
	#variants_to_keep = 0
	#total_variants = 0
	#gt2string = { -1: './.', 0: '0/0', 1: '0/1', 2: '1/1' }
	#for variant_table in vcf_reader_freebayes:
		#chromosome = variant_table.chromosome
		#print('Processing chromosome', chromosome, file=sys.stderr)
		#variant_tables_10x = []
		#for m in variants_10x:
			#try:
				#vt = m[chromosome]
			#except KeyError:
				#print('10X VCF lacking chromosome',chromosome, file=sys.stderr)
				#return 1
			#var2index = { variant:i for i, variant in enumerate(vt.variants) }
			#variant_tables_10x.append( (vt, var2index) )
		#genotypes_freebayes = []
		#for sample in samples:
			#genotypes_freebayes.append(variant_table.genotypes_of(sample))
		#samples2index = { sample:i for i, sample in enumerate(samples) }
		#genotypes_10x = [None] * len(samples)
		#for vt, var2index in variant_tables_10x:
			#for sample in vt.samples:
				#g = vt.genotypes_of(sample)
				#i = samples2index[sample]
				#genotypes_10x[i] = (g, var2index)
		
		#for variant_index, variant in enumerate(variant_table.variants):
			##print('Considering variant', variant, file=sys.stderr)
			#total_variants += 1
			#keep_variant = True
			#gt_strings = []
			#for sample_index, sample in enumerate(samples):
				#gt_freebayes = genotypes_freebayes[sample_index][variant_index]
				#if gt_freebayes == -1:
					#keep_variant = False
				#gt_list_10x, var2index10x = genotypes_10x[sample_index]
				#if variant in var2index10x:
					#gt_10x = gt_list_10x[var2index10x[variant]]
				#else:
					#gt_10x = -2
				#if (gt_freebayes == 0) and (gt_10x < 0):
					#gt_strings.append('0/0')
				#elif gt_freebayes == gt_10x:
					#gt_strings.append(gt2string[gt_freebayes])
				#else:
					#keep_variant = False
				#concordance_table[(gt_freebayes,gt_10x)] += 1
				##print('  sample {}: freebayes: {}, 10X: {}'.format(sample, gt_freebayes, gt_10x), file=sys.stderr)
			#if keep_variant:
				#print(chromosome, variant.position + 1, '.', variant.reference_allele, variant.alternative_allele, '.', '.', '.', 'GT', '\t'.join(gt_strings), sep='\t')
				#variants_to_keep += 1
	#print('Genotype concordance table (-1: genotype "./." or ".", -2: not in file at all)', file=sys.stderr)
	#keys = list(sorted(concordance_table.keys()))
	#for key in keys:
		#print(key, concordance_table[key], file=sys.stderr)
	#print('Keeping {} of {} variants'.format(variants_to_keep, total_variants), file=sys.stderr)
	

if __name__ == '__main__':
	main()
