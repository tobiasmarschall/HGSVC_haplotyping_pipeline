families = ['SH032', 'Y117', 'PR05']
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

rule master:
	input: 
		expand('whatshap/{source}/SH032.{chromosome}.vcf', chromosome=chromosomes, source=['10X','strandseq']),
		expand('compare/multi/{sample}/{chromosome}.tsv', chromosome=chromosomes, sample=samples['SH032']),
		'compare/pedmec10X-strandseq/SH032.tsv',
		expand('stats/{source}/SH032.tsv', source=['10X-raw', '10X-filtered', 'whatshap-10X'])
		#expand('compare/pedmec10X-strandseq/{sample}/{chromosome}.tsv', chromosome=chromosomes, sample=samples['SH032']),

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
		'~/scm/whatshap/bin/whatshap phase --ped {input.ped} --recombination-list {output.recomb} -o {output.vcf} {input.vcf_consensus} {input.vcfs_10X} > {log} 2>&1'

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
		'~/scm/whatshap/bin/whatshap phase --ped {input.ped} --recombination-list {output.recomb} -o {output.vcf} {input.vcf_consensus} {input.vcfs_strandseq} > {log} 2>&1'

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
		'~/scm/whatshap/bin/whatshap phase --ped {input.ped} --recombination-list {output.recomb} -o {output.vcf} {input.vcf_consensus} {input.vcfs_strandseq} {input.vcfs_10X} > {log} 2>&1'

rule compare_phasing:
	input:
		pedmec10X=lambda wildcards: 'whatshap/10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		strandseq='input/strandseq/{sample}/{chromosome}.vcf'
	output: 'compare/pedmec10X-strandseq/{sample}/{chromosome}.tsv'
	log: 'compare/pedmec10X-strandseq/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap/bin/whatshap compare --tsv-pairwise {output} --names pedmec10X,strandseq {input.pedmec10X} {input.strandseq} > {log} 2>&1'

rule compare_phasing_multi:
	input:
		pedmec_strandseq_10X=lambda wildcards: 'whatshap/strandseq-10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		pedmec_10X=lambda wildcards: 'whatshap/10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		pedmec_strandseq=lambda wildcards: 'whatshap/strandseq/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		raw_strandseq='strandseq/{sample}/{chromosome}.vcf'
	output: 'compare/multi/{sample}/{chromosome}.tsv'
	log: 'compare/multi/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap/bin/whatshap compare --tsv-pairwise {output} --names raw_ss,pedmec_10X,pedmec_ss,pedmec_ss_10X {input.raw_strandseq} {input.pedmec_10X} {input.pedmec_strandseq} {input.pedmec_strandseq_10X} > {log} 2>&1'

rule merge_tsvs:
	input: lambda wildcards: ['compare/pedmec10X-strandseq/{}/{}.tsv'.format(sample,chromosome) for sample in samples[wildcards.family] for chromosome in chromosomes]
	output: 'compare/pedmec10X-strandseq/{family}.tsv'
	shell:
		'(head -n1 {input[0]} && tail -q -n1 {input}) > {output}'


rule stats_10X:
	input: '10X/{sample}/{what}/{chromosome}.vcf.gz'
	output: 'stats/10X-{what}/{sample}/{chromosome}.tsv'
	log: 'stats/10X-{what}/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap/bin/whatshap stats --tsv {output} {input} > {log} 2>&1'


rule stats_whatshap:
	input: lambda wildcards: 'whatshap/{}/{}.{}.vcf'.format(wildcards.what, sample2fam[wildcards.sample], wildcards.chromosome)
	output: 'stats/whatshap-{what}/{sample}/{chromosome}.tsv'
	log: 'stats/whatshap-{what}/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap/bin/whatshap stats --tsv {output} --sample {wildcards.sample} {input} > {log} 2>&1'


rule merge_stats_tsvs:
	input: lambda wildcards: ['stats/{}/{}/{}.tsv'.format(wildcards.source,sample,chromosome) for sample in samples[wildcards.family] for chromosome in chromosomes]
	output: 'stats/{source}/{family}.tsv'
	shell:
		'(head -n1 {input[0]} && tail -q -n1 {input}) > {output}'
