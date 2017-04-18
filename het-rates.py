#!/usr/bin/env python3

from cyvcf2 import VCF
from whatshap.args import HelpfulArgumentParser as ArgumentParser


def main():
	parser = ArgumentParser(prog='het-rates.py', description=__doc__)
	#parser.add_argument('--only-genotype', default=None,
		#help='Restrict analysis to given genotype from {"0/0", "0/1", "1/1", "."}')
	parser.add_argument('sample', metavar='SAMPLE', help='sample to be read')
	parser.add_argument('vcf', metavar='VCF', help='VCF file')
	args = parser.parse_args()

	bin_width = 100000
	chromosome = None
	vcf = VCF(args.vcf, gts012=True)
	#print(vcf.samples)
	sample2idx = dict((s,i) for i,s in enumerate(vcf.samples)) 
	sample_idx = sample2idx[args.sample]
	#vcf.set_samples(args.sample)
	for variant in vcf:
		if not variant.is_snp:
			continue
		#print(variant)
		#print(variant.REF, variant.ALT)
		#print(variant.gt_types)
		#print(variant.var_type)
		#print(variant.format('GT'))
		if variant.CHROM != chromosome:
			if chromosome is not None:
				print(chromosome, bin_start+1, bin_start+bin_width, hets, homs)
			bin_start = 0
			hets = 0
			homs = 0
			chromosome = variant.CHROM
		while variant.start >= bin_start + bin_width:
			print(chromosome, bin_start+1, bin_start+bin_width, hets, homs)
			bin_start += bin_width
			hets = 0
			homs = 0
		if variant.gt_types[sample_idx] == 1:
			hets += 1
		elif variant.gt_types[sample_idx] == 2:
			homs += 1


# R stuff:
#d = read.table('t', col.names=c('chr','start','end','hets','homs'), colClasses=c("character","integer","integer","integer","integer"))
#pdf('hets.pdf')
#hist(d$hets, breaks=200, xlim=c(0,400), xlab="HETs per 100kb")
#dev.off()

if __name__ == '__main__':
	main()
