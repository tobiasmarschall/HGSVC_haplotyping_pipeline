import re
from itertools import chain
from collections import defaultdict, namedtuple
import os
families = ['SH032', 'PR05', 'Y117']
fam2pop = {
	'SH032':'CHS',
	'Y117': 'YRI',
	'PR05':'PUR'
}
samples = {
	'SH032': ['HG00512', 'HG00513', 'HG00514'],
	'Y117': ['NA19239', 'NA19238', 'NA19240'],
	'PR05': ['HG00731', 'HG00732', 'HG00733']
}
sample2fam = {}
for fam, sample_list in samples.items():
	for sample in sample_list:
		sample2fam[sample] = fam
ref = 'ref/GRCh38_full_analysis_set_plus_decoy_hla.fa'
chromosomes = ['chr'+str(x) for x in range(1,23)] #+ ['chrX']
#chromosomes = ['chr1','chr22']

genotype_sources = ['illumina', 'consensus']
phaseinputs = ['moleculo','10X','strandseq','10X_strandseq']
phasings = ['raw/10X', 'raw/strandseq'] + ['whatshap/{}-{}-{}'.format(g,p,m) for g in genotype_sources for p in phaseinputs for m in ['single', 'trio']]

rule master:
	input:
		expand('whatshap/{genotypes}/{phaseinput}/{mode}/{family}.{chromosome}.vcf', genotypes=genotype_sources, phaseinput=phaseinputs, mode=['single','trio'], family=families, chromosome=chromosomes),
		expand('stats/{source}/{family}.tsv', source=phasings, family=families),
		expand('comparison/all-inputs-{genotypes}-{mode}.{family}.tsv', genotypes=genotype_sources, mode=['single','trio'], family=families)
		
		#expand('whatshap/{genotypes}/{phaseinput}/single/{family}.{chromosome}.vcf', genotypes=['illumina', 'consensus'], phaseinput=['moleculo','10X','strandseq','10X_strandseq'], family=families, chromosome=['chr22'])
		#expand('moleculo/chrwise/{sample}.{chromosome}.bam', sample=chain(*(samples[f] for f in families)), chromosome=chromosomes)
		#expand('freebayes/moleculo/filtered/{family}.{chromosome}.vcf', family=families, chromosome=chromosomes)
		
		#expand('pacbio/{sample}.{chromosome}.bam', sample=['HG00731', 'HG00732', 'HG00733'], chromosome=chromosomes)
		#expand('pacbio/bwa/runs/{run}.bam', run=[r for r, info in pacbio_runid_urls.items() if info.sample=='HG00731'])
		#expand('compare/multi5/{family}.tsv', chromosome=chromosomes, family=families),
		#expand('stats/{source}/{family}.tsv', family=families, source=['10X-raw', '10X-filtered', 'whatshap-10X', 'whatshap-strandseq', 'whatshap-strandseq-10X', 'strandseq'])
		
		#expand('strandseq/projected/{sample}/chr1.vcf', sample=chain(*(samples[f] for f in families))),
		#expand('tagged-bam/illumina/{sample}.tagged-20160704.bam.bai', sample=chain(*(samples[f] for f in families)))

		#expand('whatshap/{source}/{family}.{chromosome}.vcf', family=families, chromosome=chromosomes, source=['10X','strandseq']),
		#expand('compare/multi3/{family}.tsv', chromosome=chromosomes, family=families),
		#expand('release/{family}.wgs.whatshap.strandseq-10X.20160622.phased-genotypes.vcf.gz', family=families),
		##expand('compare/pedmec10X-strandseq/{sample}/{chromosome}.tsv', chromosome=chromosomes, sample=samples['SH032']),
		#expand('consensus/freebayes_10X/{family}.{chromosome}.vcf', chromosome=chromosomes, family=['Y117']),
		#expand('ftp/working/20160623_chaisson_pacbio_aligns/20160525.{sample}.PacBio.BLASR/{sample}.PUR.{chromosome}.BLASR.20160525.PacBio.10x-phased.bam.bai',sample=['HG00731', 'HG00732', 'HG00733'], chromosome=chromosomes),
		#expand('ftp/working/20160623_chaisson_pacbio_aligns/20160525.{sample}.PacBio.BLASR/{sample}.PUR.{chromosome}.BLASR.20160525.PacBio.10x-phased.bam',sample=['HG00731', 'HG00732', 'HG00733'], chromosome=chromosomes),
		#expand('whatshap/pacbio-noped/PR05.{chromosome}.vcf', chromosome=chromosomes),
		#expand('whatshap/pacbio/PR05.{chromosome}.vcf', chromosome=chromosomes)

def translate_ebi_filename(f):
	f = re.sub('HG00732.PUR.(.*).BLASR.20160525.PacBio.10x-phased.(.*)', 'HG007732.PUR.\\1.BLASR.20160525.PacBio.10x-phased.\\2', f)
	f = re.sub('HG00731.PUR.(.*).BLASR.20160525.PacBio.10x-phased.(.*)', 'HG00731.YRI.\\1.BLASR.20160525.PacBio.10x-phased.\\2', f)
	f = re.sub('HG00731.YRI.chr([56789]).BLASR.20160525.PacBio.10x-phased.(.*)', 'HG00731.chr\\1.\\2', f)
	f = re.sub('HG00733.PUR.(.*).BLASR.20160525.PacBio.10x-phased.(.*)', 'HG00733.\\1.\\2', f)
	return f

rule download_ebi_ftp:
	output: 'ftp/{file}'
	log: 'ftp/{file}.wgetlog'
	run: 
		f = translate_ebi_filename('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/{}'.format(wildcards.file))
		shell('wget --output-file={log} -O {output} {f}')
	#shell: 'touch {output}'


rule download_moleculo:
	output: 'moleculo/download/{sample}.bam'
	log: 'moleculo/download/{sample}.bam.wgetlog'
	shell: 'wget --output-file={log} -O {output} ftp://ftp.sra.ebi.ac.uk/vol1/ERA690/ERA690446/bam/{wildcards.sample}_sorted.bam'

rule moleculo_readgroup:
	input: 
		bam='moleculo/download/{sample}.bam'
	output:
		bam='moleculo/bam/{sample}.bam',
		bai='moleculo/bam/{sample}.bai'
	log: 'moleculo/bam/{sample}.log'
	shell: 'picard AddOrReplaceReadGroups CREATE_INDEX=true CREATE_MD5_FILE=true ID={wildcards.sample} PL=moleculo PU=unknown LB=unknown SM={wildcards.sample} I={input.bam} O={output.bam} > {log} 2>&1'


rule moleculo_extract_chromosome:
	input:
		bam='moleculo/bam/{sample}.bam',
		bai='moleculo/bam/{sample}.bai'
	output:
		bam='moleculo/chrwise/{sample}.{chromosome}.bam'
	log: 'moleculo/chrwise/{sample}.{chromosome}.bam.log'
	shell: '(samtools view -h {input.bam} {wildcards.chromosome} | samtools view -Sb - > {output.bam}) 2>{log}'

ruleorder: download_ebi_ftp > index_bam

rule download_ref:
	output: 'ref/GRCh38_full_analysis_set_plus_decoy_hla.fa'
	log: 'ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.wgetlog'
	shell: 'wget --directory-prefix=ref --output-file={log} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'


rule extract_strandseq:
	input: 'ftp/working/20160627_StrandSeq_vcfs/VCFfiles_StrandSeq_{sample}.zip'
	output: expand('strandseq/raw/{{sample}}/{chromosome}_phased_{{sample}}.vcf', chromosome=chromosomes)
	log: 'strandseq/raw/{sample}/unzip.log'
	shell:
		'unzip -d strandseq/raw/{wildcards.sample} {input} > {log} 2>&1; '
		'sed -i \'/^#CHROM/s/NA12878/{wildcards.sample}/g\' {output}; '
		'touch {output}'

rule freebayes_illumina:
	input:
		crams=lambda wildcards: ['ftp/data/{pop}/{sample}/high_cov_alignment/{sample}.alt_bwamem_GRCh38DH.20150715.{pop}.high_coverage.cram'.format(sample=s, pop=fam2pop[wildcards.family]) for s in samples[wildcards.family]],
		crais=lambda wildcards: ['ftp/data/{pop}/{sample}/high_cov_alignment/{sample}.alt_bwamem_GRCh38DH.20150715.{pop}.high_coverage.cram.crai'.format(sample=s, pop=fam2pop[wildcards.family]) for s in samples[wildcards.family]]
	output:
		vcf='freebayes/illumina/raw/{family}.{chromosome}.vcf'
	shell:
		'samtools merge -R {wildcards.chromosome} - {input.crams} | freebayes -f {ref} --stdin > {output.vcf}'


rule freebayes_moleculo:
	input:
		bams=lambda wildcards: ['moleculo/bam/{sample}.bam'.format(sample=s) for s in samples[wildcards.family]],
		bais=lambda wildcards: ['moleculo/bam/{sample}.bai'.format(sample=s) for s in samples[wildcards.family]],
	output:
		vcf='freebayes/moleculo/raw/{family}.{chromosome}.vcf'
	shell:
		'samtools merge -R {wildcards.chromosome} - {input.bams} | freebayes -f {ref} --stdin > {output.vcf}'


rule freebayes_pacbio_retype:
	input:
		vcf='consensus/freebayes_10X/{family}.{chromosome}.sitesonly.vcf.gz',
		vcfidx='consensus/freebayes_10X/{family}.{chromosome}.sitesonly.vcf.gz.tbi',
		bams=lambda wildcards: ['pacbio/{sample}.{chromosome}.bam'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		bais=lambda wildcards: ['pacbio/{sample}.{chromosome}.bai'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]]
	output:
		vcf='freebayes-pacbio-retype/raw/{family}.{chromosome}.vcf'
	log: 'freebayes-pacbio-retype/raw/{family}.{chromosome}.vcf.log'
	shell:
		'(time samtools merge -R {wildcards.chromosome} - {input.bams} |  freebayes -f {ref} --haplotype-basis-alleles {input.vcf} -@ {input.vcf} --stdin > {output.vcf} ) 2> {log}'

rule filter_freebayes_vcf:
	input:
		vcf='freebayes/{what}/raw/{family}.{chromosome}.vcf'
	output:
		vcf='freebayes/{what}/filtered/{family}.{chromosome}.vcf'
	shell:
		'awk \'($0 ~ /^#/) || (($6 != ".") && ($6 >= 30))\' {input.vcf} > {output.vcf}'

rule filter_10X_vcf:
	input:
		vcf='10X/{sample}/raw/{chromosome}.vcf.gz'
	output:
		vcf='10X/{sample}/filtered/{chromosome}.vcf.gz'
	shell:
		'zcat {input.vcf} | awk \'($0 ~ /^#/) || (($6 != ".") && ($6 >= 30) && ($7 == "PASS"))\' | gzip > {output.vcf}'


rule consensus_freebayes_10X:
	input:
		vcfs_10X=lambda wildcards: ['10X/{sample}/filtered/{chromosome}.vcf.gz'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		vcf_freebayes='freebayes/illumina/filtered/{family}.{chromosome}.vcf'
	output:
		vcf='consensus/freebayes_10X/{family}.{chromosome}.vcf'
	log: 'consensus/freebayes_10X/{family}.{chromosome}.log'
	shell:
		'PYTHONPATH=~/scm/whatshap ~/scm/hgsvc/consensus-genotypes-10X-freebayes.py {input.vcf_freebayes} {input.vcfs_10X} > {output.vcf} 2> {log}'

rule consensus_sitesonly:
	input:
		vcf='consensus/{what}/{family}.{chromosome}.vcf'
	output:
		vcf='consensus/{what}/{family}.{chromosome}.sitesonly.vcf.gz'
	log: 'consensus/{what}/{family}.{chromosome}.sitesonly.vcf.log'
	shell: 'picard MakeSitesOnlyVcf I={input.vcf} O={output.vcf} CREATE_INDEX=false > {log} 2>&1'

ruleorder: consensus_sitesonly > consensus_freebayes_10X

rule tabix:
	input: '{file}.vcf.gz'
	output: '{file}.vcf.gz.tbi'
	shell: 'tabix {input}'

rule extract_chromosome_10X:
	input: 
		vcf='ftp/working/20160513_10XGenomics_data_and_calls/{sample}_WGS/{sample}_WGS_phased_variants.vcf.gz'
	output:
		vcf='10X/{sample}/raw/{chromosome}.vcf.gz'
	shell:
		'zcat {input.vcf} | awk \'($0 ~ /^#/) || ($1=="{wildcards.chromosome}")\' | uniq | gzip > {output.vcf}'


rule project_strandseq:
	input:
		vcf_strandseq='strandseq/raw/{sample}/{chromosome}_phased_{sample}.vcf',
		vcf_consensus=lambda wildcards: 'consensus/freebayes_10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample], chromosome=wildcards.chromosome)
	output:
		vcf='strandseq/projected/{sample}/{chromosome}.vcf'
	log: 'strandseq/projected/{sample}/{chromosome}.log'
	shell:
		'~/scm/hgsvc/project-strandseq-to-consensus.py --max-pq 50 {input.vcf_consensus} {input.vcf_strandseq} > {output.vcf} 2> {log}'

rule extract_ped:
	input: 'ftp/sv_trio_families.ped'
	output: 'ped/{family}.ped'
	shell: 'awk \'$1 == "{wildcards.family}"\' {input} > {output}'


def get_genotype_vcf(wildcards):
	if wildcards.genotypes == 'illumina':
		return 'freebayes/illumina/filtered/{family}.{chromosome}.vcf'.format(**wildcards)
	elif wildcards.genotypes == 'consensus':
		return 'consensus/freebayes_10X/{family}.{chromosome}.vcf'.format(**wildcards)
	else:
		assert False


def get_phase_input(wildcards):
	result = []
	for p in wildcards.phaseinput.split('_'):
		if p == 'moleculo':
			result.extend(['moleculo/bam/{sample}.bam'.format(sample=s) for s in samples[wildcards.family]])
		elif p == '10X': 
			result.extend(['10X/{sample}/filtered/{chromosome}.vcf.gz'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]])
		elif p == 'strandseq':
			result.extend(['strandseq/projected/{sample}/{chromosome}.vcf'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]])
		else:
			assert False
	return result

# =========== WHATSHAP ===========
rule whatshap_single:
	input:
		vcf=get_genotype_vcf,
		phaseinput=get_phase_input,
		ref=ref
	output:
		vcf='whatshap/{genotypes}/{phaseinput}/single/{family}.{chromosome}.vcf',
		corrected_gts='whatshap/{genotypes}/{phaseinput}/single/{family}.{chromosome}.genotype-changes.tsv',
		readlist='whatshap/{genotypes}/{phaseinput}/single/{family}.{chromosome}.used-reads.tsv'
	log: 'whatshap/{genotypes}/{phaseinput}/single/{family}.{chromosome}.vcf.log'
	shell: '~/scm/whatshap.to-run/bin/whatshap phase --indels --distrust-genotypes --output-read-list {output.readlist} --changed-genotype-list {output.corrected_gts} --reference {input.ref} -o {output.vcf} {input.vcf} {input.phaseinput} > {log} 2>&1'


rule whatshap_trio:
	input:
		vcf=get_genotype_vcf,
		phaseinput=get_phase_input,
		ref=ref,
		ped='ped/{family}.ped'
	output:
		vcf='whatshap/{genotypes}/{phaseinput}/trio/{family}.{chromosome}.vcf',
		corrected_gts='whatshap/{genotypes}/{phaseinput}/trio/{family}.{chromosome}.genotype-changes.tsv',
		recomb='whatshap/{genotypes}/{phaseinput}/trio/{family}.{chromosome}.recomb',
		readlist='whatshap/{genotypes}/{phaseinput}/single/{family}.{chromosome}.used-reads.tsv'
	log: 'whatshap/{genotypes}/{phaseinput}/trio/{family}.{chromosome}.vcf.log'
	shell: '~/scm/whatshap.to-run/bin/whatshap phase --ped {input.ped} --indels --distrust-genotypes --output-read-list {output.readlist} --changed-genotype-list {output.corrected_gts} --reference {input.ref} --recombination-list {output.recomb} -o {output.vcf} {input.vcf} {input.phaseinput} > {log} 2>&1'


#rule whatshap_10X:
	#input:
		#vcfs_10X=lambda wildcards: ['10X/{sample}/filtered/{chromosome}.vcf.gz'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		#vcf_consensus='consensus/freebayes_10X/{family}.{chromosome}.vcf',
		#ped='ped/{family}.ped'
	#output:
		#vcf='whatshap/10X/{family}.{chromosome}.vcf',
		#recomb='whatshap/10X/{family}.{chromosome}.recomb',
	#log: 'whatshap/10X/{family}.{chromosome}.log'
	#shell:
		#'~/scm/whatshap.to-run/bin/whatshap phase --ped {input.ped} --tag PS --recombination-list {output.recomb} -o {output.vcf} {input.vcf_consensus} {input.vcfs_10X} > {log} 2>&1'

#rule whatshap_strandseq:
	#input:
		#vcfs_strandseq=lambda wildcards: ['strandseq/projected/{sample}/{chromosome}.vcf'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		#vcf_consensus='consensus/freebayes_10X/{family}.{chromosome}.vcf',
		#ped='ped/{family}.ped'
	#output:
		#vcf='whatshap/strandseq/{family}.{chromosome}.vcf',
		#recomb='whatshap/strandseq/{family}.{chromosome}.recomb',
	#log: 'whatshap/strandseq/{family}.{chromosome}.log'
	#shell:
		#'~/scm/whatshap.to-run/bin/whatshap phase --ped {input.ped} --tag PS --recombination-list {output.recomb} -o {output.vcf} {input.vcf_consensus} {input.vcfs_strandseq} > {log} 2>&1'

#rule whatshap_strandseq_10X:
	#input:
		#vcfs_10X=lambda wildcards: ['10X/{sample}/filtered/{chromosome}.vcf.gz'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		#vcfs_strandseq=lambda wildcards: ['strandseq/projected/{sample}/{chromosome}.vcf'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		#vcf_consensus='consensus/freebayes_10X/{family}.{chromosome}.vcf',
		#ped='ped/{family}.ped'
	#output:
		#vcf='whatshap/strandseq-10X/{family}.{chromosome}.vcf',
		#recomb='whatshap/strandseq-10X/{family}.{chromosome}.recomb',
	#log: 'whatshap/strandseq-10X/{family}.{chromosome}.log'
	#shell:
		#'~/scm/whatshap.to-run/bin/whatshap phase --ped {input.ped} --tag PS --recombination-list {output.recomb} -o {output.vcf} {input.vcf_consensus} {input.vcfs_strandseq} {input.vcfs_10X} > {log} 2>&1'


#rule whatshap_pacbio:
	#input:
		#vcf_consensus='consensus/freebayes_10X/{family}.{chromosome}.vcf',
		#ref=ref,
		#bams=lambda wildcards: ['pacbio/{sample}.{chromosome}.bam'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		#bais=lambda wildcards: ['pacbio/{sample}.{chromosome}.bai'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		#ped='ped/{family}.ped'
	#output:
		#vcf='whatshap/pacbio/{family}.{chromosome}.vcf',
		#recomb='whatshap/pacbio/{family}.{chromosome}.recomb',
	#log: 'whatshap/pacbio/{family}.{chromosome}.log'
	#shell:
		#'~/scm/whatshap.to-run/bin/whatshap phase --ped {input.ped} --tag PS --reference {input.ref} --recombination-list {output.recomb} -o {output.vcf} {input.vcf_consensus} {input.bams} > {log} 2>&1'


#rule whatshap_pacbio_noped:
	#input:
		#vcf_consensus='consensus/freebayes_10X/{family}.{chromosome}.vcf',
		#ref=ref,
		#bams=lambda wildcards: ['pacbio/{sample}.{chromosome}.bam'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		#bais=lambda wildcards: ['pacbio/{sample}.{chromosome}.bai'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
	#output:
		#vcf='whatshap/pacbio-noped/{family}.{chromosome}.vcf'
	#log: 'whatshap/pacbio-noped/{family}.{chromosome}.log'
	#shell:
		#'~/scm/whatshap.to-run/bin/whatshap phase --tag PS --reference {input.ref} -o {output.vcf} {input.vcf_consensus} {input.bams} > {log} 2>&1'

# =========== COMPARISON ===========
rule compare_phasing:
	input: 
		vcfs= lambda wildcards: ['whatshap/{}/{}/{}/{}.{}.vcf'.format(wildcards.genotypes,p,wildcards.mode,sample2fam[wildcards.sample],wildcards.chromosome) for p in phaseinputs]
	output: 
		tsv='comparison/all-inputs-{genotypes}-{mode}/{sample}.{chromosome}.tsv'
	log: 'comparison/all-inputs-{genotypes}-{mode}/{sample}.{chromosome}.log'
	run:
		names = ','.join(phaseinputs)
		shell('~/scm/whatshap.to-run/bin/whatshap compare --sample {wildcards.sample} --tsv-pairwise {output.tsv} --names {names} {input.vcfs} > {log} 2>&1')


rule merge_tsvs:
	input: lambda wildcards: ['comparison/all-inputs-{}-{}/{}.{}.tsv'.format(wildcards.genotypes,wildcards.mode,sample,chromosome) for sample in samples[wildcards.family] for chromosome in chromosomes]
	output: 'comparison/all-inputs-{genotypes}-{mode}.{family}.tsv'
	shell:
		'(head -n1 {input[0]} && grep -hv "^#" {input}) > {output}'


#rule compare_phasing:
	#input:
		#pedmec10X=lambda wildcards: 'whatshap/10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		#strandseq='strandseq/raw/{sample}/{chromosome}.vcf'
	#output: 'compare/pedmec10X-strandseq/{sample}/{chromosome}.tsv'
	#log: 'compare/pedmec10X-strandseq/{sample}/{chromosome}.log'
	#shell:
		#'~/scm/whatshap.to-run/bin/whatshap compare --tsv-pairwise {output} --names pedmec10X,strandseq {input.pedmec10X} {input.strandseq} > {log} 2>&1'

#rule compare_phasing_multi5:
	#input:
		#pedmec_strandseq_10X=lambda wildcards: 'whatshap/strandseq-10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		#pedmec_10X=lambda wildcards: 'whatshap/10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		#pedmec_strandseq=lambda wildcards: 'whatshap/strandseq/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		#only_strandseq='strandseq/projected/{sample}/{chromosome}.vcf',
		#only_10X='10X/{sample}/filtered/{chromosome}.vcf.gz'
	#output: 'compare/multi5/{sample}/{chromosome}.tsv'
	#log: 'compare/multi5/{sample}/{chromosome}.log'
	#shell:
		#'~/scm/whatshap.to-run/bin/whatshap compare --tsv-pairwise {output} --names only_10X,only_ss,pedmec_10X,pedmec_ss,pedmec_ss_10X {input.only_10X} {input.only_strandseq} {input.pedmec_10X} {input.pedmec_strandseq} {input.pedmec_strandseq_10X} > {log} 2>&1'

#rule compare_phasing_multi3:
	#input:
		#pedmec_strandseq_10X=lambda wildcards: 'whatshap/strandseq-10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		#pedmec_10X=lambda wildcards: 'whatshap/10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		#pedmec_strandseq=lambda wildcards: 'whatshap/strandseq/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
	#output: 'compare/multi3/{sample}/{chromosome}.tsv'
	#log: 'compare/multi3/{sample}/{chromosome}.log'
	#shell:
		#'~/scm/whatshap.to-run/bin/whatshap compare --tsv-pairwise {output} --sample {wildcards.sample} --names pedmec_10X,pedmec_ss,pedmec_ss_10X {input.pedmec_10X} {input.pedmec_strandseq} {input.pedmec_strandseq_10X} > {log} 2>&1'

#rule merge_tsvs:
	#input: lambda wildcards: ['compare/{}/{}/{}.tsv'.format(wildcards.what, sample,chromosome) for sample in samples[wildcards.family] for chromosome in chromosomes]
	#output: 'compare/{what}/{family}.tsv'
	#shell:
		#'(head -n1 {input[0]} && grep -hv "^#" {input}) > {output}'

# =========== STATISTICS ===========
rule stats_raw_strandseq:
	input: 'strandseq/projected/{sample}/{chromosome}.vcf'
	output: 'stats/raw/strandseq/{sample}/{chromosome}.tsv'
	log: 'stats/raw/strandseq/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap stats --tsv {output} {input} > {log} 2>&1'


rule stats_10X:
	input: '10X/{sample}/filtered/{chromosome}.vcf.gz'
	output: 'stats/raw/10X/{sample}/{chromosome}.tsv'
	log: 'stats/raw/10X/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap stats --tsv {output} {input} > {log} 2>&1'


rule stats_whatshap:
	input: lambda wildcards: 'whatshap/{}/{}/{}/{}.{}.vcf'.format(wildcards.genotypes, wildcards.phaseinput, wildcards.what, sample2fam[wildcards.sample], wildcards.chromosome)
	output: 'stats/whatshap/{genotypes,[^/]+}-{phaseinput,[^/]+}-{what,[^/]+}/{sample,[^/]+}/{chromosome,[^/]+}.tsv'
	log: 'stats/whatshap/{genotypes,[^/]+}-{phaseinput,[^/]+}-{what,[^/]+}/{sample,[^/]+}/{chromosome,[^/]+}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap stats --tsv {output} --sample {wildcards.sample} {input} > {log} 2>&1'


rule merge_stats_tsvs:
	input: lambda wildcards: ['stats/{}/{}/{}/{}.tsv'.format(wildcards.method, wildcards.source,sample,chromosome) for sample in samples[wildcards.family] for chromosome in chromosomes]
	output: 'stats/{method,[^/]+}/{source,[^/]+}/{family,[^/]+}.tsv'
	shell:
		'(head -n1 {input[0]} && tail -q -n1 {input}) > {output}'


rule sequence_dict:
	input: 'ref/{file}.fa'
	output: 'ref/{file}.fa.dict'
	log: 'ref/{file}.fa.dict.log'
	shell:
		'picard CreateSequenceDictionary R={input} O={output} > {log} 2>&1'


rule filter_whatshap:
	input: 'whatshap/{what}/{family}.{chromosome}.vcf'
	output: 'whatshap/{what}/only-phased-sites/{family}.{chromosome}.vcf'
	shell:
		'awk \'($0 ~ /^#/) || ($9 ~ /PS/)\' {input} > {output}'


rule merge_vcfs_for_release:
	input:
		vcf= lambda wildcards: ['whatshap/{}/only-phased-sites/{}.{}.vcf'.format(wildcards.what,wildcards.family, chromosome) for chromosome in chromosomes],
		seqdict='{}.dict'.format(ref)
	output:
		vcf='release/{family}.wgs.whatshap.{what}.{date}.phased-genotypes.vcf.gz',
		tbi='release/{family}.wgs.whatshap.{what}.{date}.phased-genotypes.vcf.gz.tbi'
	log: 'release/{family}.wgs.whatshap.{what}.{date}.phased-genotypes.vcf.log'
	run:
		allinputs = ' '.join('I={}'.format(f) for f in input.vcf)
		shell('picard SortVcf {allinputs} SD={input.seqdict} O={output.vcf} > {log} 2>&1')


rule illumina_cram_to_bam:
	input:
		cram= lambda wildcards: 'ftp/data/{pop}/{sample}/high_cov_alignment/{sample}.alt_bwamem_GRCh38DH.20150715.{pop}.high_coverage.cram'.format(sample=wildcards.sample, pop=fam2pop[sample2fam[wildcards.sample]])
	output:
		bam='illumina/{sample}.bam'
	log: 'illumina/{sample}.bam.log'
	shell:
		'samtools view -h {input.cram} | samtools view -Sb - > {output.bam} 2> {log}'


rule index_bam:
	input: 
		bam='{file}.bam'
	output:
		bai='{file}.bam.bai'
	shell:
		'samtools index {input.bam}'


rule haplotag_illumina:
	input:
		bam='illumina/{sample}.bam',
		bai='illumina/{sample}.bam.bai',
		vcf=lambda wildcards:'release/{family}.wgs.whatshap.strandseq-10X.{date}.phased-genotypes.vcf.gz'.format(family=sample2fam[wildcards.sample], date=wildcards.date)
	output:
		bam='tagged-bam/illumina/{sample}.tagged-{date}.bam'
	log: 'tagged-bam/illumina/{sample}.tagged-{date}.bam.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap haplotag {input.vcf} {input.bam} > {output.bam} 2> {log}'


ruleorder: pacbio_reheader_sourcemanual > pacbio_reheader
ruleorder: pacbio_correct_header_sourcemanual > pacbio_correct_header


rule pacbio_correct_header_sourcemanual:
	input: 'pacbio-blasr-missing/links/{sample}.{chromosome}.bam'
	output: 'pacbio/{sample}.{chromosome}.header'
	shell:
		'samtools view -H {input} | sed -r \'/^@RG/ s|SM:[a-z0-9]+|SM:{wildcards.sample}|g\' > {output}'


rule pacbio_reheader_sourcemanual:
	input:
		header='pacbio/{sample}.{chromosome}.header',
		bam='pacbio-blasr-missing/links/{sample}.{chromosome}.bam'
	output:
		bam='pacbio/{sample}.{chromosome}.bam',
		bai='pacbio/{sample}.{chromosome}.bai',
	log: 'pacbio/{sample}.{chromosome}.bam.log'
	shell:
		'picard ReplaceSamHeader CREATE_INDEX=true CREATE_MD5_FILE=true I={input.bam} O={output.bam} HEADER={input.header} > {log} 2>&1'


rule pacbio_correct_header:
	input: lambda wildcards: 'ftp/working/20160623_chaisson_pacbio_aligns/20160525.{sample}.PacBio.BLASR/{sample}.{pop}.{chromosome}.BLASR.20160525.PacBio.10x-phased.bam'.format(sample=wildcards.sample, chromosome=wildcards.chromosome, pop=fam2pop[sample2fam[wildcards.sample]])
	output: 'pacbio/{sample}.{chromosome}.header'
	shell:
		'samtools view -H {input} | sed -r \'/^@RG/ s|SM:[a-z0-9]+|SM:{wildcards.sample}|g\' > {output}'


rule pacbio_reheader:
	input:
		header='pacbio/{sample}.{chromosome}.header',
		bam=lambda wildcards: 'ftp/working/20160623_chaisson_pacbio_aligns/20160525.{sample}.PacBio.BLASR/{sample}.{pop}.{chromosome}.BLASR.20160525.PacBio.10x-phased.bam'.format(sample=wildcards.sample, chromosome=wildcards.chromosome, pop=fam2pop[sample2fam[wildcards.sample]])
	output:
		bam='pacbio/{sample}.{chromosome}.bam',
		bai='pacbio/{sample}.{chromosome}.bai',
	log: 'pacbio/{sample}.{chromosome}.bam.log'
	shell:
		'picard ReplaceSamHeader CREATE_INDEX=true CREATE_MD5_FILE=true I={input.bam} O={output.bam} HEADER={input.header} > {log} 2>&1'
