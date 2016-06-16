#!/usr/bin/env python3
import sys
import datetime
import math
from collections import defaultdict
import vcf
from whatshap.args import HelpfulArgumentParser as ArgumentParser
from whatshap.vcf import VcfVariant

vcf_header = '''##fileformat=VCFv4.1
##fileDate={date}
##source=project-strandseq-to-consensus.py
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PQ,Number=1,Type=Integer,Description="Phred-scaled phasing quality">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}'''


def find_sample_id(needle, haystack):
	for i, x in enumerate(haystack):
		if needle == x:
			return i
	return None

def all_none(l):
	for x in l:
		if x is not None: return False
	return True


def read_synced_by_pos(vcf_readers):
	iterators = [ iter(vcf_reader) for vcf_reader in vcf_readers ]
	nexts = [ next(it) for it in iterators ]
	chroms = set( record.CHROM for record in nexts if record is not None )
	assert len(chroms) == 1
	chrom = list(chroms)[0]
	output = [ [] for _ in range(len(vcf_readers)) ]
	while not all_none(nexts):
		position = min( record.POS for record in nexts if record is not None)
		for i in range(len(vcf_readers)):
			while (nexts[i] is not None) and (nexts[i].POS == position):
				output[i].append(nexts[i])
				nexts[i] = next(iterators[i])
				assert (nexts[i] is None) or (nexts[i].CHROM == chrom)
		yield output
		output = [ [] for _ in range(len(vcf_readers)) ]


def main():
	parser = ArgumentParser(prog='project-strandseq-to-consensus.py', description=__doc__)
	parser.add_argument('--max-pq', metavar='MAX_PQ', default=50, type=int,
		help='Max value for PQ (default: %(default)s).')
	parser.add_argument('vcfconsensus', metavar='CONSENSUS_VCF', help='VCF file with consensus genotypes')
	parser.add_argument('vcfstrandseq', metavar='STRANDSEQ_VCF', help='VCF file from StrandSeq')
	args = parser.parse_args()
	
	vcf_reader_consensus = vcf.Reader(filename=args.vcfconsensus)
	vcf_reader_strandseq = vcf.Reader(filename=args.vcfstrandseq)

	samples = set(vcf_reader_consensus.samples).intersection(vcf_reader_strandseq.samples)
	assert len(samples) == 1
	sample = list(samples)[0]
	consensus_sample_id = find_sample_id(sample, vcf_reader_consensus.samples)
	strandseq_sample_id = find_sample_id(sample, vcf_reader_strandseq.samples)
	print('Working on sample', sample, file=sys.stderr)
	#self.samples = self._vcf_reader.samples
	consensus_only = 0
	strandseq_only = 0
	common = 0
	retained = 0
	gtmap = {
		'0|.':'0|1',
		'.|0':'1|0',
		'1|.':'1|0',
		'.|1':'0|1',
		'0|1':'0|1',
		'1|0':'1|0'
	}
	print(vcf_header.format(date=datetime.datetime.now().strftime('%Y%m%d'), sample=sample))
	for records_consensus, records_strandseq in read_synced_by_pos([vcf_reader_consensus, vcf_reader_strandseq]):
		if (len(records_consensus) == 0) and (len(records_strandseq) > 0):
			strandseq_only += 1
			continue
		if (len(records_consensus) > 0) and (len(records_strandseq) == 0):
			consensus_only += 1
			continue
		assert len(records_consensus) == 1
		assert len(records_strandseq) == 1
		common += 1
		assert records_consensus[0].CHROM == records_strandseq[0].CHROM
		assert records_consensus[0].POS == records_strandseq[0].POS
		assert records_consensus[0].REF == records_strandseq[0].REF
		assert len(records_consensus[0].ALT) == 1
		if len(records_strandseq[0].ALT) > 1:
			continue
		if (records_consensus[0].ALT[0] != records_strandseq[0].ALT[0]) and (records_strandseq[0].ALT[0] != 'N'):
			continue
		call_consensus = records_consensus[0].samples[consensus_sample_id]
		call_strandseq = records_strandseq[0].samples[strandseq_sample_id]
		if not call_consensus.is_het:
			continue
		if not call_strandseq.data.GT in ['0|.', '.|0', '1|.', '.|1', '0|1',  '1|0']:
			continue
		p1 = call_strandseq.data.P1
		p2 = call_strandseq.data.P2
		p_wrong = 1.0
		if p1 is not None: p_wrong *= 1.0 - p1
		if p2 is not None: p_wrong *= 1.0 - p2
		pq = -10 * math.log10(p_wrong) if p_wrong > 0 else args.max_pq
		pq = round(min(pq, args.max_pq))
		#print(records_consensus[0].POS, p1, p2, pq)
		gt = gtmap[call_strandseq.data.GT]
		r = records_consensus[0]
		print(r.CHROM, r.POS, '.', r.REF, r.ALT[0], '.', 'PASS', '.', 'GT:PQ', '{}:{}'.format(gt, pq), sep='\t')
		retained += 1
	print('Records in StrandSeq only: ', strandseq_only, file=sys.stderr)
	print('Records in consensus only: ', consensus_only, file=sys.stderr)
	print('{} common records of which {} have been retained'.format(common, retained), file=sys.stderr)


if __name__ == '__main__':
	main()
