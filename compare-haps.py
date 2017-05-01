#!/usr/bin/env python3

from collections import defaultdict
from cyvcf2 import VCF
from whatshap.args import HelpfulArgumentParser as ArgumentParser



def hamming(s0, s1):
	assert len(s0) == len(s1)
	return sum( c0!=c1 for c0, c1 in zip(s0, s1) )

def main():
	parser = ArgumentParser(prog='het-rates.py', description=__doc__)
	#parser.add_argument('--only-genotype', default=None,
		#help='Restrict analysis to given genotype from {"0/0", "0/1", "1/1", "."}')
	#parser.add_argument('sample', metavar='SAMPLE', help='sample to be read')
	parser.add_argument('vcf', metavar='VCF', help='VCF file')
	parser.add_argument('region', metavar='REGION', help='chr:start-end')
	args = parser.parse_args()

	#bin_width = 100000
	#chromosome = None
	vcf = VCF(args.vcf, gts012=True)
	#print(vcf.samples)
	#sample2idx = dict((s,i) for i,s in enumerate(vcf.samples)) 
	#sample_idx = sample2idx[args.sample]
	#vcf.set_samples(args.sample)
	hap = [defaultdict(list), defaultdict(list)]
	for variant in vcf(args.region):
		if not variant.is_snp:
			continue
		for sample, (hap0_allele, hap1_allele, phased) in zip(vcf.samples,variant.genotypes):
			print(sample, hap0_allele, hap1_allele, phased)
			if (hap0_allele + hap1_allele == 1) and (not phased):
				hap[0][sample].append('?')
				hap[1][sample].append('?')
			else:
				hap[0][sample].append(hap0_allele)
				hap[1][sample].append(hap1_allele)


	all_haps = [(sample, i) for sample in vcf.samples for i in [0,1]]
	for sample, i in all_haps:
		print(sample, 'hap{}'.format(i), ''.join(str(x) for x in hap[i][sample]))
		print(sample, 'hap{}'.format(i), ''.join(str(x) for x in hap[i][sample]))
	print()
	
	for i in range(len(all_haps)):
		sample0, hap0 = all_haps[i]
		for j in range(i+1, len(all_haps)):
			sample1, hap1 = all_haps[j]
			h = hamming(hap[hap0][sample0],hap[hap1][sample1])
			print(sample0, hap0, '<-->', sample1, hap1, h, h/len(hap[hap0][sample0]))

# R stuff:
#d = read.table('t', col.names=c('chr','start','end','hets','homs'), colClasses=c("character","integer","integer","integer","integer"))
#pdf('hets.pdf')
#hist(d$hets, breaks=200, xlim=c(0,400), xlab="HETs per 100kb")
#dev.off()

if __name__ == '__main__':
	main()
