#!/usr/bin/env python3
import sys
import datetime
from collections import defaultdict
from whatshap.args import HelpfulArgumentParser as ArgumentParser
from whatshap.vcf import VcfReader

vcf_header = '''##fileformat=VCFv4.1
##fileDate={date}
##source=consensus-genotypes-10X-freebayes.py
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
'''

def main():
	parser = ArgumentParser(prog='consensus-genotypes-10X-freebayes.py', description=__doc__)
	parser.add_argument('vcffreebayes', metavar='FREEBAYES_VCF', help='VCF file with variants genotyped by freebayes')
	parser.add_argument('vcf10x', nargs='+', metavar='10X_VCF', help='VCF file with variants genotyped by 10X genomics')
	args = parser.parse_args()
	vcf_reader_freebayes = VcfReader(args.vcffreebayes)
	vcf_readers_10x = [VcfReader(f) for f in args.vcf10x]

	samples_freebayes = set(vcf_reader_freebayes.samples)
	samples_10x = set()
	for vcf_reader in vcf_readers_10x:
		assert len(set(vcf_reader.samples).intersection(samples_10x)) == 0
		samples_10x.update(vcf_reader.samples)
	assert samples_freebayes == samples_10x
	samples = list(sorted(samples_freebayes))
	print('Samples:', ', '.join(samples), file=sys.stderr)

	# list of dicts mapping chromosome names to VariantTables for each input VCF
	variants_10x = []
	for vcf_reader, filename in zip(vcf_readers_10x, args.vcf10x):
		m = {}
		for variant_table in vcf_reader:
			m[variant_table.chromosome] = variant_table
		variants_10x.append(m)
		print('Successfully read', filename, file=sys.stderr)

	print(vcf_header.format(date=datetime.datetime.now().strftime('%Y%m%d'), sample='\t'.join(samples)))

	concordance_table = defaultdict(int)
	variants_to_keep = 0
	total_variants = 0
	gt2string = { -1: './.', 0: '0/0', 1: '0/1', 2: '1/1' }
	for variant_table in vcf_reader_freebayes:
		chromosome = variant_table.chromosome
		print('Processing chromosome', chromosome)
		variant_tables_10x = []
		for m in variants_10x:
			try:
				vt = m[chromosome]
			except KeyError:
				print('10X VCF lacking chromosome',chromosome, file=sys.stderr)
				return 1
			var2index = { variant:i for i, variant in enumerate(vt.variants) }
			variant_tables_10x.append( (vt, var2index) )
		genotypes_freebayes = []
		for sample in samples:
			genotypes_freebayes.append(variant_table.genotypes_of(sample))
		samples2index = { sample:i for i, sample in enumerate(samples) }
		genotypes_10x = [None] * len(samples)
		for vt, var2index in variant_tables_10x:
			for sample in vt.samples:
				g = vt.genotypes_of(sample)
				i = samples2index[sample]
				genotypes_10x[i] = (g, var2index)
		
		for variant_index, variant in enumerate(variant_table.variants):
			#print('Considering variant', variant, file=sys.stderr)
			total_variants += 1
			keep_variant = True
			gt_strings = []
			for sample_index, sample in enumerate(samples):
				gt_freebayes = genotypes_freebayes[sample_index][variant_index]
				if gt_freebayes == -1:
					keep_variant = False
				gt_list_10x, var2index10x = genotypes_10x[sample_index]
				if variant in var2index10x:
					gt_10x = gt_list_10x[var2index10x[variant]]
				else:
					gt_10x = -2
				if (gt_freebayes == 0) and (gt_10x < 0):
					gt_strings.append('0/0')
				elif gt_freebayes == gt_10x:
					gt_strings.append(gt2string[gt_freebayes])
				else:
					keep_variant = False
				concordance_table[(gt_freebayes,gt_10x)] += 1
				#print('  sample {}: freebayes: {}, 10X: {}'.format(sample, gt_freebayes, gt_10x), file=sys.stderr)
			if keep_variant:
				print(chromosome, variant.position + 1, '.', variant.reference_allele, variant.alternative_allele, '.', '.', '.', 'GT', '\t'.join(gt_strings), sep='\t')
				variants_to_keep += 1
	keys = list(sorted(concordance_table.keys()))
	for key in keys:
		print(key, concordance_table[key])
	print('Keeping {} of {} variants'.format(variants_to_keep, total_variants), file=sys.stderr)
	

if __name__ == '__main__':
	main()
