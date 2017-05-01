#!/usr/bin/env python3

import sys
from collections import defaultdict
from cyvcf2 import VCF
from whatshap.args import HelpfulArgumentParser as ArgumentParser



def hamming(s0, s1):
	assert len(s0) == len(s1)
	return sum( ((c0!=c1) and (c0!='?') and (c1!='?'))  for c0, c1 in zip(s0, s1) )

def main():
	parser = ArgumentParser(prog='het-rates.py', description=__doc__)
	#parser.add_argument('--only-genotype', default=None,
		#help='Restrict analysis to given genotype from {"0/0", "0/1", "1/1", "."}')
	parser.add_argument('--names', metavar='NAMES', default=None, help='Comma-separated list of data set names')
	#parser.add_argument('sample', metavar='SAMPLE', help='sample to be read')
	parser.add_argument('region', metavar='REGION', help='chr:start-end')
	parser.add_argument('vcf', nargs='+', metavar='VCF', help='VCF file')
	args = parser.parse_args()

	if args.names:
		names = args.names.split(',')
	else:
		names = ['file{}'.format(i) for i in range(len(args.vcf))]
	assert len(names) == len(args.vcf)

	variant_intersection = set()
	# check which variants are present
	for i, filename in enumerate(args.vcf):
		vcf = VCF(filename, gts012=True)
		s = set()
		for variant in vcf(args.region):
			if not variant.is_snp:
				continue
			s.add( (variant.start, variant.REF, variant.ALT[0]) )
		print('Found', len(s), 'SNVs in', filename, file=sys.stderr)
		if i == 0:
			variant_intersection = s
		else:
			variant_intersection.intersection_update(s)
	print('Found', len(variant_intersection), 'SNVs in intersection', file=sys.stderr)

	chromosome, region_str = args.region.split(':')
	start, end = [int(x) for x in region_str.split('-')]
	length = end - start + 1
	
	#bin_width = 100000
	#chromosome = None
	#vcf = VCF(args.vcf, gts012=True)
	#print(vcf.samples)
	#sample2idx = dict((s,i) for i,s in enumerate(vcf.samples)) 
	#sample_idx = sample2idx[args.sample]
	#vcf.set_samples(args.sample)
	hap = []
	for i, filename in enumerate(args.vcf):
		h = [defaultdict(list), defaultdict(list)]
		vcf = VCF(filename, gts012=True)
		for variant in vcf(args.region):
			if (variant.start, variant.REF, variant.ALT[0]) not in variant_intersection:
				continue
			for sample, (hap0_allele, hap1_allele, phased) in zip(vcf.samples,variant.genotypes):
				#print(sample, hap0_allele, hap1_allele, phased)
				if (hap0_allele + hap1_allele == 1) and (not phased):
					h[0][sample].append('?')
					h[1][sample].append('?')
				else:
					h[0][sample].append(hap0_allele)
					h[1][sample].append(hap1_allele)
		hap.append(h)

	print('HET RATES', file=sys.stderr)
	for name_id, name in enumerate(names):
		for sample in vcf.samples:
			hets = sum( ((a0==0 and a1==1) or (a0==1 and a1==0)) for a0, a1 in zip(hap[name_id][0][sample],hap[name_id][1][sample]))
			print(name, sample, 'het_rate:', hets/length*100000, file=sys.stderr)

	print(file=sys.stderr)
	print('HAPLOTYPES', file=sys.stderr)

	all_haps = [(name_id,sample, i) for name_id in range(len(names)) for sample in vcf.samples for i in [0,1]]
	for name_id, sample, i in all_haps:
		print(names[name_id], sample, 'hap{}'.format(i), ''.join(str(x) for x in hap[name_id][i][sample]), file=sys.stderr)
		print(names[name_id], sample, 'hap{}'.format(i), ''.join(str(x) for x in hap[name_id][i][sample]), file=sys.stderr)
	print(file=sys.stderr)
	
	for i in range(len(all_haps)):
		name_id0, sample0, hap0 = all_haps[i]
		for j in range(i+1, len(all_haps)):
			name_id1, sample1, hap1 = all_haps[j]
			hamming_distance = hamming(hap[name_id0][hap0][sample0],hap[name_id1][hap1][sample1])
			n = len(hap[name_id0][hap0][sample0])
			print(names[name_id0],sample0, 'hap{}'.format(hap0), '<-->', names[name_id1], sample1, 'hap{}'.format(hap1), hamming_distance, '{0:.2f}'.format(hamming_distance/n*100.0))


if __name__ == '__main__':
	main()
