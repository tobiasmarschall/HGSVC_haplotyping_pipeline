from itertools import chain
families = ['SH032', 'PR05'] # 'Y117'
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
#chromosomes = ['chr1']

rule master:
	input: 
		expand('whatshap/{source}/{family}.{chromosome}.vcf', family=families, chromosome=chromosomes, source=['10X','strandseq']),
		expand('compare/multi3/{family}.tsv', chromosome=chromosomes, family=families),
		expand('compare/multi5/{family}.tsv', chromosome=chromosomes, family=families),
		expand('stats/{source}/{family}.tsv', family=families, source=['10X-raw', '10X-filtered', 'whatshap-10X', 'whatshap-strandseq', 'whatshap-strandseq-10X', 'strandseq']),
		expand('release/{family}.wgs.whatshap.strandseq-10X.20160622.phased-genotypes.vcf.gz', family=families),
		#expand('compare/pedmec10X-strandseq/{sample}/{chromosome}.tsv', chromosome=chromosomes, sample=samples['SH032']),
		expand('consensus/freebayes_10X/{family}.{chromosome}.vcf', chromosome=chromosomes, family=['Y117']),
		expand('tagged-bam/illumina/{sample}.tagged-20160622.bam.bai', sample=chain(*(samples[f] for f in families)))


rule download_ebi_ftp:
	output: 'ftp/{file}'
	log: 'ftp/{file}.wgetlog'
	shell: 'wget -r -x --no-host-directories --cut-dirs=4 --no-clobber --directory-prefix=ftp --output-file={log} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/{wildcards.file}'
	#shell: 'touch {output}'


rule download_ref:
	output: 'ref/GRCh38_full_analysis_set_plus_decoy_hla.fa'
	log: 'ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.wgetlog'
	shell: 'wget --directory-prefix=ref --output-file={log} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'


rule freebayes:
	input:
		crams=lambda wildcards: ['ftp/data/{pop}/{sample}/high_cov_alignment/{sample}.alt_bwamem_GRCh38DH.20150715.{pop}.high_coverage.cram'.format(sample=s, pop=fam2pop[wildcards.family]) for s in samples[wildcards.family]],
		crais=lambda wildcards: ['ftp/data/{pop}/{sample}/high_cov_alignment/{sample}.alt_bwamem_GRCh38DH.20150715.{pop}.high_coverage.cram.crai'.format(sample=s, pop=fam2pop[wildcards.family]) for s in samples[wildcards.family]]
	output:
		vcf='freebayes/raw/{family}.{chromosome}.vcf'
	shell:
		'samtools merge -R {wildcards.chromosome} - {input.crams} | freebayes -f {ref} --stdin > {output.vcf}'


rule filter_freebayes_vcf:
	input:
		vcf='freebayes/raw/{family}.{chromosome}.vcf'
	output:
		vcf='freebayes/filtered/{family}.{chromosome}.vcf'
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
		vcf_freebayes='freebayes/filtered/{family}.{chromosome}.vcf'
	output:
		vcf='consensus/freebayes_10X/{family}.{chromosome}.vcf'
	log: 'consensus/freebayes_10X/{family}.{chromosome}.log'
	shell:
		'PYTHONPATH=~/scm/whatshap ~/scm/hgsvc/consensus-genotypes-10X-freebayes.py {input.vcf_freebayes} {input.vcfs_10X} > {output.vcf} 2> {log}'


rule extract_chromosome_10X:
	input: 
		vcf='ftp/working/20160513_10XGenomics_data_and_calls/{sample}_WGS/{sample}_WGS_phased_variants.vcf.gz'
	output:
		vcf='10X/{sample}/raw/{chromosome}.vcf.gz'
	shell:
		'zcat {input.vcf} | awk \'($0 ~ /^#/) || ($1=="{wildcards.chromosome}")\' | uniq | gzip > {output.vcf}'


rule project_strandseq:
	input:
		vcf_strandseq='input/strandseq/{sample}/{chromosome}.vcf',
		vcf_consensus=lambda wildcards: 'consensus/freebayes_10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample], chromosome=wildcards.chromosome)
	output:
		vcf='strandseq/{sample}/{chromosome}.vcf'
	log: 'strandseq/{sample}/{chromosome}.log'
	shell:
		'~/scm/hgsvc/project-strandseq-to-consensus.py --max-pq 50 {input.vcf_consensus} {input.vcf_strandseq} > {output.vcf} 2> {log}'

rule extract_ped:
	input: 'ftp/sv_trio_families.ped'
	output: 'ped/{family}.ped'
	shell: 'awk \'$1 == "{wildcards.family}"\' {input} > {output}'

rule whatshap_10X:
	input:
		vcfs_10X=lambda wildcards: ['10X/{sample}/filtered/{chromosome}.vcf.gz'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		vcf_consensus='consensus/freebayes_10X/{family}.{chromosome}.vcf',
		ped='ped/{family}.ped'
	output:
		vcf='whatshap/10X/{family}.{chromosome}.vcf',
		recomb='whatshap/10X/{family}.{chromosome}.recomb',
	log: 'whatshap/10X/{family}.{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap phase --ped {input.ped} --tag PS --recombination-list {output.recomb} -o {output.vcf} {input.vcf_consensus} {input.vcfs_10X} > {log} 2>&1'

rule whatshap_strandseq:
	input:
		vcfs_strandseq=lambda wildcards: ['strandseq/{sample}/{chromosome}.vcf'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		vcf_consensus='consensus/freebayes_10X/{family}.{chromosome}.vcf',
		ped='ped/{family}.ped'
	output:
		vcf='whatshap/strandseq/{family}.{chromosome}.vcf',
		recomb='whatshap/strandseq/{family}.{chromosome}.recomb',
	log: 'whatshap/strandseq/{family}.{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap phase --ped {input.ped} --tag PS --recombination-list {output.recomb} -o {output.vcf} {input.vcf_consensus} {input.vcfs_strandseq} > {log} 2>&1'

rule whatshap_strandseq_10X:
	input:
		vcfs_10X=lambda wildcards: ['10X/{sample}/filtered/{chromosome}.vcf.gz'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		vcfs_strandseq=lambda wildcards: ['strandseq/{sample}/{chromosome}.vcf'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		vcf_consensus='consensus/freebayes_10X/{family}.{chromosome}.vcf',
		ped='ped/{family}.ped'
	output:
		vcf='whatshap/strandseq-10X/{family}.{chromosome}.vcf',
		recomb='whatshap/strandseq-10X/{family}.{chromosome}.recomb',
	log: 'whatshap/strandseq-10X/{family}.{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap phase --ped {input.ped} --tag PS --recombination-list {output.recomb} -o {output.vcf} {input.vcf_consensus} {input.vcfs_strandseq} {input.vcfs_10X} > {log} 2>&1'

rule compare_phasing:
	input:
		pedmec10X=lambda wildcards: 'whatshap/10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		strandseq='input/strandseq/{sample}/{chromosome}.vcf'
	output: 'compare/pedmec10X-strandseq/{sample}/{chromosome}.tsv'
	log: 'compare/pedmec10X-strandseq/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap compare --tsv-pairwise {output} --names pedmec10X,strandseq {input.pedmec10X} {input.strandseq} > {log} 2>&1'

rule compare_phasing_multi5:
	input:
		pedmec_strandseq_10X=lambda wildcards: 'whatshap/strandseq-10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		pedmec_10X=lambda wildcards: 'whatshap/10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		pedmec_strandseq=lambda wildcards: 'whatshap/strandseq/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		only_strandseq='strandseq/{sample}/{chromosome}.vcf',
		only_10X='10X/{sample}/filtered/{chromosome}.vcf.gz'
	output: 'compare/multi5/{sample}/{chromosome}.tsv'
	log: 'compare/multi5/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap compare --tsv-pairwise {output} --names only_10X,only_ss,pedmec_10X,pedmec_ss,pedmec_ss_10X {input.only_10X} {input.only_strandseq} {input.pedmec_10X} {input.pedmec_strandseq} {input.pedmec_strandseq_10X} > {log} 2>&1'

rule compare_phasing_multi3:
	input:
		pedmec_strandseq_10X=lambda wildcards: 'whatshap/strandseq-10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		pedmec_10X=lambda wildcards: 'whatshap/10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		pedmec_strandseq=lambda wildcards: 'whatshap/strandseq/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
	output: 'compare/multi3/{sample}/{chromosome}.tsv'
	log: 'compare/multi3/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap compare --tsv-pairwise {output} --sample {wildcards.sample} --names pedmec_10X,pedmec_ss,pedmec_ss_10X {input.pedmec_10X} {input.pedmec_strandseq} {input.pedmec_strandseq_10X} > {log} 2>&1'

rule merge_tsvs:
	input: lambda wildcards: ['compare/{}/{}/{}.tsv'.format(wildcards.what, sample,chromosome) for sample in samples[wildcards.family] for chromosome in chromosomes]
	output: 'compare/{what}/{family}.tsv'
	shell:
		'(head -n1 {input[0]} && grep -hv "^#" {input}) > {output}'


rule stats_strandseq:
	input: 'strandseq/{sample}/{chromosome}.vcf'
	output: 'stats/strandseq/{sample}/{chromosome}.tsv'
	log: 'stats/strandseq/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap stats --tsv {output} {input} > {log} 2>&1'


rule stats_10X:
	input: '10X/{sample}/{what}/{chromosome}.vcf.gz'
	output: 'stats/10X-{what}/{sample}/{chromosome}.tsv'
	log: 'stats/10X-{what}/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap stats --tsv {output} {input} > {log} 2>&1'


rule stats_whatshap:
	input: lambda wildcards: 'whatshap/{}/{}.{}.vcf'.format(wildcards.what, sample2fam[wildcards.sample], wildcards.chromosome)
	output: 'stats/whatshap-{what}/{sample}/{chromosome}.tsv'
	log: 'stats/whatshap-{what}/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap.to-run/bin/whatshap stats --tsv {output} --sample {wildcards.sample} {input} > {log} 2>&1'


rule merge_stats_tsvs:
	input: lambda wildcards: ['stats/{}/{}/{}.tsv'.format(wildcards.source,sample,chromosome) for sample in samples[wildcards.family] for chromosome in chromosomes]
	output: 'stats/{source}/{family}.tsv'
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
