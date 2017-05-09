#!/usr/bin/env python3
import sys
import datetime
from collections import defaultdict, namedtuple
from whatshap.args import HelpfulArgumentParser as ArgumentParser
import vcf
from array import array
#from whatshap.vcf import VcfReader


# ====================================================================================================
class VcfNotSortedError(Exception):
	pass


class VcfVariant:
	"""A variant in a VCF file (not to be confused with core.Variant)"""

	__slots__ = ('position', 'reference_allele', 'alternative_allele')

	def __init__(self, position, reference_allele, alternative_allele):
		"""
		position -- 0-based start coordinate
		reference_allele -- string
		alternative_allele -- string

		Multi-ALT sites are not modelled.
		"""
		self.position = position
		self.reference_allele = reference_allele
		self.alternative_allele = alternative_allele

	def __repr__(self):
		return "VcfVariant({}, {!r}, {!r})".format(self.position,
			self.reference_allele, self.alternative_allele)

	def __hash__(self):
		return hash((self.position, self.reference_allele, self.alternative_allele))

	def __eq__(self, other):
		return (self.position == other.position) and \
		       (self.reference_allele == other.reference_allele) and \
		       (self.alternative_allele == other.alternative_allele)

	def __lt__(self, other):
		return (self.position, self.reference_allele, self.alternative_allele) < (other.position, other.reference_allele, other.alternative_allele)

	def is_snv(self):
		return (self.reference_allele != self.alternative_allele) and (
			len(self.reference_allele) == len(self.alternative_allele) == 1)

class GenotypeLikelihoods:
	__slots__ = ('log_prob_g0', 'log_prob_g1', 'log_prob_g2')

	def __init__(self, log_prob_g0, log_prob_g1, log_prob_g2):
		"""Likelihoods of the three genotypes 0, 1, 2 to be given
		as log10 of the original probability."""
		self.log_prob_g0 = log_prob_g0
		self.log_prob_g1 = log_prob_g1
		self.log_prob_g2 = log_prob_g2

	def __repr__(self):
		return "GenotypeLikelihoods({}, {}, {})".format(self.log_prob_g0, self.log_prob_g1, self.log_prob_g2)

	def log10_probs(self):
		return ( self.log_prob_g0, self.log_prob_g1, self.log_prob_g2 )

	def log10_prob_of(self, genotype):
		return self.log10_probs()[genotype]

	def as_phred(self, regularizer=None):
		if regularizer is None:
			# shift log likelihoods such that the largest one is zero
			m = max(self.log_prob_g0, self.log_prob_g1, self.log_prob_g2)
			return PhredGenotypeLikelihoods(
				round((self.log_prob_g0-m) * -10),
				round((self.log_prob_g1-m) * -10),
				round((self.log_prob_g2-m) * -10)
			)
		else:
			p = [ 10**x for x in (self.log_prob_g0, self.log_prob_g1, self.log_prob_g2) ]
			s = sum(p)
			p = [ x/s + regularizer for x in p ]
			m = max(p)
			return PhredGenotypeLikelihoods( *(round(-10*math.log10(x/m)) for x in p) )


class VariantTable:
	"""
	For a single chromosome, store variants and their genotypes.
	Each row of this table contains a variant, each column
	contains the genotypes of a single sample.

	chromosome -- chromosome name
	samples -- list of sample names
	"""
	def __init__(self, chromosome, samples):
		self.chromosome = chromosome
		self.samples = samples
		self.genotypes = [ array('b', []) for _ in samples ]
		self.phases = [ [] for _ in samples ]
		self.genotype_likelihoods = [ [] for _ in samples ]
		self.allele_depths = [ [] for _ in samples ]
		self.allele_quality_sums = [ [] for _ in samples ]
		self.variants = []
		self._sample_to_index = { sample: index for index, sample in enumerate(samples) }

	def __len__(self):
		return len(self.variants)

	#def add_sample(self, name, genotypes):
		#"Add a column to the table"
		#if len(genotypes) != len(self.variants):
			#raise ValueError('Expecting as many genotypes as there are variants')
		#self._name_to_index[name] = len(self.samples)
		#self.samples.append(name)
		#self.genotypes.append(genotypes)

	def add_variant(self, variant, genotypes, phases, genotype_likelihoods, allele_depths, allele_quality_sums):
		"""
		Add a row to the table

		variant -- a VcfVariant
		genotypes -- iterable of ints that encode the genotypes of the samples:
			-1 represents an unknown genotype
			0 represents 0/0 (homozygous reference)
			1 represents 0/1 or 1/0 (heterozygous)
			2 represents 1/1 (homozygous alternative)
		phases -- iterable of VariantCallPhase objects
		genotype_likelihoods -- iterable of GenotypeLikelihoods objects
		"""
		if len(genotypes) != len(self.genotypes):
			raise ValueError('Expecting as many genotypes as there are samples')
		if len(phases) != len(self.phases):
			raise ValueError('Expecting as many phases as there are samples')
		self.variants.append(variant)
		for i, genotype in enumerate(genotypes):
			self.genotypes[i].append(genotype)
		for i, phase in enumerate(phases):
			self.phases[i].append(phase)
		for i, gl in enumerate(genotype_likelihoods):
			self.genotype_likelihoods[i].append(gl)
		for i, ad in enumerate(allele_depths):
			self.allele_depths[i].append(ad)
		for i, aqs in enumerate(allele_quality_sums):
			self.allele_quality_sums[i].append(aqs)

	def genotypes_of(self, sample):
		"""Retrieve genotypes by sample name"""
		return self.genotypes[self._sample_to_index[sample]]

	def genotype_likelihoods_of(self, sample):
		"""Retrieve genotype likelihoods by sample name"""
		return self.genotype_likelihoods[self._sample_to_index[sample]]

	def allele_depths_of(self, sample):
		return self.allele_depths[self._sample_to_index[sample]]

	def allele_quality_sums_of(self, sample):
		return self.allele_quality_sums[self._sample_to_index[sample]]

	def phases_of(self, sample):
		"""Retrieve phases by sample name"""
		return self.phases[self._sample_to_index[sample]]

	def num_of_blocks_of(self, sample):
		""" Retrieve the number of blocks of the sample"""
		return len(set([i.block_id for i in self.phases[self._sample_to_index[sample]] if i is not None]))

	def id_of(self, sample):
		"""Return a unique int id of a sample given by name"""
		return self._sample_to_index[sample]


class MixedPhasingError(Exception):
	pass


VariantCallPhase = namedtuple('VariantCallPhase', ['block_id', 'phase', 'quality'])
VariantCallPhase.__doc__ = \
	"""
	block_id is a numeric id of the phased block and
	phase is either 0 or 1 (indicating whether the REF allele is on haplotype 0 or 1).
	"""


class VcfReader:
	"""
	Read a VCF file chromosome by chromosome.
	"""
	def __init__(self, path, indels=False, phases=False, genotype_likelihoods=False):
		"""
		path -- Path to VCF file
		indels -- Whether to include also insertions and deletions in the list of
			variants.
		"""
		# TODO Always include deletions since they can 'overlap' other variants
		self._indels = indels
		self._vcf_reader = vcf.Reader(filename=path)
		self._phases = phases
		self._genotype_likelihoods = genotype_likelihoods
		self.samples = self._vcf_reader.samples  # intentionally public
		#logger.debug("Found %d sample(s) in the VCF file.", len(self.samples))

	def _group_by_chromosome(self):
		"""
		Yield (chromosome, records) tuples, where records is a list of the
		VCF records on that chromosome.
		"""
		records = []
		prev_chromosome = None
		for record in self._vcf_reader:
			if record.CHROM != prev_chromosome:
				if prev_chromosome is not None:
					yield (prev_chromosome, records)
				prev_chromosome = record.CHROM
				records = []
			records.append(record)
		if records:
			yield (prev_chromosome, records)

	def __iter__(self):
		"""
		Yield VariantTable objects for each chromosome.

		Multi-ALT sites are skipped.
		"""
		for chromosome, records in self._group_by_chromosome():
			yield self._process_single_chromosome(chromosome, records)

	@staticmethod
	def _extract_HP_phase(call):
		HP = getattr(call.data, 'HP', None)
		if HP is None:
			return None
		assert len(HP) == 2
		fields = [[int(x) for x in s.split('-')] for s in HP]
		assert fields[0][0] == fields[1][0]
		block_id = fields[0][0]
		phase1, phase2 = fields[0][1]-1, fields[1][1]-1
		assert ((phase1, phase2) == (0, 1)) or ((phase1, phase2) == (1, 0))
		return VariantCallPhase(block_id=block_id, phase=phase1, quality=getattr(call.data, 'PQ', None))

	@staticmethod
	def _extract_GT_PS_phase(call):
		if not call.is_het:
			return None
		if not call.phased:
			return None
		block_id = getattr(call.data, 'PS', 0)
		if block_id is None:
			block_id = 0
		assert call.data.GT in ['0|1','1|0']
		phase = int(call.data.GT[0])
		return VariantCallPhase(block_id=block_id, phase=phase, quality=getattr(call.data, 'PQ', None))

	def _process_single_chromosome(self, chromosome, records):
		phase_detected = None
		n_snvs = 0
		n_other = 0
		n_multi = 0
		table = VariantTable(chromosome, self.samples)
		prev_position = None
		for record in records:
			if len(record.ALT) > 1:
				# Multi-ALT sites are not supported, yet
				n_multi += 1
				continue

			pos, ref, alt = record.start, str(record.REF), str(record.ALT[0])
			if len(ref) == len(alt) == 1:
				n_snvs += 1
			else:
				n_other += 1
				if not self._indels:
					continue

			if (prev_position is not None) and (prev_position > pos):
				raise VcfNotSortedError('VCF not ordered: {}:{} appears before {}:{}'.format(chromosome, prev_position+1, chromosome, pos+1))

			if prev_position == pos:
				#logger.warning('Skipping duplicated position %s on chromosome %r', pos+1, chromosome)
				continue
			prev_position = pos

			# Read phasing information (allow GT/PS or HP phase information, but not both),
			# if requested
			if self._phases:
				phases = []
				for call in record.samples:
					phase = None
					for extract_phase, phase_name in [(self._extract_HP_phase, 'HP'), (self._extract_GT_PS_phase, 'GT_PS')]:
						p = extract_phase(call)
						if p is not None:
							if phase_detected is None:
								phase_detected = phase_name
							elif phase_detected != phase_name:
								raise MixedPhasingError('Mixed phasing information in input VCF (e.g. mixing PS and HP fields)')
							phase = p
					phases.append(phase)
			else:
				phases = [ None ] * len(record.samples)

			# Read genotype likelihoods, if requested
			if self._genotype_likelihoods:
				genotype_likelihoods = []
				for call in record.samples:
					GL = getattr(call.data, 'GL', None)
					PL = getattr(call.data, 'PL', None)
					# Prefer GLs (floats) over PLs (ints) if both should be present
					if GL is not None:
						assert len(GL) == 3
						genotype_likelihoods.append(GenotypeLikelihoods(*GL))
					elif PL is not None:
						assert len(PL) == 3
						genotype_likelihoods.append(GenotypeLikelihoods( *(pl/-10 for pl in PL) ))
					else:
						genotype_likelihoods.append(None)
			else:
				genotype_likelihoods = [ None ] * len(record.samples)

			allele_depths = []
			allele_quality_sums = []
			for call in record.samples:
				RO = getattr(call.data, 'RO', None)
				AO = getattr(call.data, 'AO', None)
				QR = getattr(call.data, 'QR', None)
				QA = getattr(call.data, 'QA', None)
				if (RO is None) or (AO is None):
					allele_depths.append(None)
				else:
					allele_depths.append((RO,AO))
				if (QR is None) or (QA is None):
					allele_quality_sums.append(None)
				else:
					allele_quality_sums.append((QR,QA))

			# PyVCF pecularity: gt_alleles is a list of the alleles in the
			# GT field, but as strings.
			# For example, when GT is 0/1, gt_alleles is ['0', '1'].
			# And when GT is 2|1, gt_alleles is ['2', '1'].
			GT_TO_INT = { 0: 0, 1: 1, 2: 2, None: -1 }
			genotypes = array('b', (GT_TO_INT[call.gt_type] for call in record.samples))
			variant = VcfVariant(position=pos, reference_allele=ref, alternative_allele=alt)
			table.add_variant(variant, genotypes, phases, genotype_likelihoods, allele_depths, allele_quality_sums)

		#logger.debug("Parsed %s SNVs and %s non-SNVs. Also skipped %s multi-ALTs.", n_snvs,
			#n_other, n_multi)

		# TODO remove overlapping variants
		return table



# ====================================================================================================

vcf_header = '''##fileformat=VCFv4.1
##fileDate={date}
##source=consensus-genotypes-freebayes-gatk.py
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}'''

# FreeBayes formats
##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">

def main():
	parser = ArgumentParser(prog='consensus-genotypes-10X-freebayes.py', description=__doc__)
	parser.add_argument('--include-ad', action='store_true', default=False, help='Add allele depth information')
	parser.add_argument('vcffreebayes', metavar='FREEBAYES_VCF', help='VCF file with variants genotyped by freebayes')
	parser.add_argument('vcfgatk', metavar='GATK_VCF', help='VCF file with variants genotyped by GATK')
	args = parser.parse_args()
	vcf_reader_freebayes = VcfReader(args.vcffreebayes)
	vcf_reader_gatk = VcfReader(args.vcfgatk)
	samples_freebayes = set(vcf_reader_freebayes.samples)
	samples_gatk = set(vcf_reader_gatk.samples)
	assert samples_freebayes == samples_gatk
	samples = list(sorted(samples_freebayes))
	print('Samples:', ', '.join(samples), file=sys.stderr)

	# dict mapping chromosome names to VariantTables for GATK variants
	variants_gatk = {}
	for variant_table in vcf_reader_gatk:
		variants_gatk[variant_table.chromosome] = variant_table
	print('Successfully read', args.vcfgatk, file=sys.stderr)

	print(vcf_header.format(date=datetime.datetime.now().strftime('%Y%m%d'), sample='\t'.join(samples)))

	concordance_table = defaultdict(int)
	variants_to_keep = 0
	total_variants = 0
	gt2string = { -1: './.', 0: '0/0', 1: '0/1', 2: '1/1' }
	for variant_table in vcf_reader_freebayes:
		chromosome = variant_table.chromosome
		print('Processing chromosome', chromosome, file=sys.stderr)
		try:
			variant_table_gatk = variants_gatk[chromosome]
		except KeyError:
			print('GATK VCF lacking chromosome',chromosome, file=sys.stderr)
			return 1
		var2index_gatk = { variant:i for i, variant in enumerate(variant_table_gatk.variants) }
		genotypes_freebayes = []
		for sample in samples:
			genotypes_freebayes.append(list(zip(variant_table.genotypes_of(sample), variant_table.allele_depths_of(sample), variant_table.allele_quality_sums_of(sample))))
		samples2index = { sample:i for i, sample in enumerate(samples) }

		genotypes_gatk = [None] * len(samples)
		for sample in variant_table_gatk.samples:
			g = variant_table_gatk.genotypes_of(sample)
			ad = variant_table_gatk.allele_depths_of(sample)
			aqs = variant_table_gatk.allele_quality_sums_of(sample)
			i = samples2index[sample]
			genotypes_gatk[i] = (g, ad, aqs)

		for variant_index, variant in enumerate(variant_table.variants):
			#print('Considering variant', variant, file=sys.stderr)
			total_variants += 1
			keep_variant = True
			format_fields = []
			for sample_index, sample in enumerate(samples):
				gt_freebayes, ad_freebayes, aqs_freebayes = genotypes_freebayes[sample_index][variant_index]
				if gt_freebayes == -1:
					keep_variant = False
				gt_list_gatk, ad_list_gatk, aqs_list_gatk = genotypes_gatk[sample_index]
				if variant in var2index_gatk:
					idx_gatk = var2index_gatk[variant]
					gt_gatk = gt_list_gatk[idx_gatk]
				else:
					idx_gatk = None
					gt_gatk = -2
				if (gt_freebayes == 0) and (gt_gatk < 0):
					format_fields.append(['0/0', ad_freebayes, None, aqs_freebayes, None])
				elif gt_freebayes == gt_gatk:
					format_fields.append([gt2string[gt_freebayes], ad_freebayes, ad_list_gatk[idx_gatk], aqs_freebayes, aqs_list_gatk[idx_gatk]])
				else:
					keep_variant = False
				concordance_table[(gt_freebayes,gt_gatk)] += 1
				#print('  sample {}: freebayes: {}, GATK: {}'.format(sample, gt_freebayes, gt_gatk), file=sys.stderr)
			if keep_variant:
				def s(t):
					if t is None:
						return '.'
					elif type(t) == str:
						return t
					else:
						return '{},{}'.format(*t)
				if args.include_ad:
					raise 'Not yet implemented'
					print(chromosome, variant.position + 1, '.', variant.reference_allele, variant.alternative_allele, '.', '.', '.', 'GT:DI:DX:QI:QX', '\t'.join( ':'.join(s(x) for x in f) for f in format_fields), sep='\t')
				else:
					print(chromosome, variant.position + 1, '.', variant.reference_allele, variant.alternative_allele, '.', '.', '.', 'GT', '\t'.join(f[0] for f in format_fields), sep='\t')
				variants_to_keep += 1
	print('Genotype concordance table (-1: genotype "./." or ".", -2: not in file at all)', file=sys.stderr)
	keys = list(sorted(concordance_table.keys()))
	for key in keys:
		print(key, concordance_table[key], file=sys.stderr)
	print('Keeping {} of {} variants'.format(variants_to_keep, total_variants), file=sys.stderr)
	

if __name__ == '__main__':
	main()
