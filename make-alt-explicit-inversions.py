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
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="length of the variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality: probability that the genotype is wrong (phred based).">
##FORMAT=<ID=L,Number=.,Type=Integer,Description="genotype likelihoods for absent, heterozygous, homozygous.">
##FORMAT=<ID=IE,Number=.,Type=Integer,Description="reversed read evidence in two comma-separated values in case of inversions: REF support, ALT support">
##FORMAT=<ID=SE,Number=.,Type=Integer,Description="Split-read evidence in two comma-separated values in case of inversions: REF support, ALT support">
##FORMAT=<ID=DE,Number=.,Type=Integer,Description="number of supporting and neutral reads in two comma-separated values in case of duplications: REF support, ALT support">
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

def main():
	parser = ArgumentParser(prog='make-alt-explicit-inversions.py', description=__doc__)
	parser.add_argument('ref', metavar='REF', help='Reference genome (FASTA)')
	parser.add_argument('vcf', metavar='VCF', help='VCF of merged PacBio variants')
	#parser.add_argument('--min-distance', default=0, type=int, help='Minimum distance to next variant.')
	#parser.add_argument('--max-length', default=None, type=int, help='Maxmimum length.')
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

	for record in vcf_reader:
		print('Processing', record.CHROM, record.POS, record.REF, record.ALT, 'length={}'.format(record.INFO['SVLEN']), file=sys.stderr)
		assert len(record.ALT) == 1
		assert len(record.REF) == 1
		#assert record.REF[0] == reference[record.CHROM][record.POS - 1]
		assert 'SVLEN' in record.INFO
		if str(record.ALT[0]) == '<INV>':
			ref = reference[record.CHROM][record.POS-1:record.POS+abs(record.INFO['SVLEN'])].upper()
			alt = ref[0] + reference[record.CHROM][record.POS+abs(record.INFO['SVLEN'])-1:record.POS-1:-1].upper()
			print(
				nonedot(record.CHROM),
				nonedot(record.POS),
				nonedot(record.ID),
				nonedot(ref),
				nonedot(alt),
				nonedot(record.QUAL),
				nonedot(record.FILTER),
				info_to_str(record.INFO),
				'GT',
				record.samples[0].data.GT,
			sep='\t')
		else:
			assert False, 'Unexpected ALT field: {}'.format(record.ALT[0])

if __name__ == '__main__':
	main()
