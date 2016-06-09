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
		expand('whatshap/10X/SH032.{chromosome}.vcf', chromosome=chromosomes),
		expand('compare/pedmec10X-strandseq/{sample}/{chromosome}.txt', chromosome=chromosomes, sample=samples['SH032'])

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


rule extract_ped:
	input: 'ftp/sv_trio_families.ped'
	output: 'ped/{family}.ped'
	shell: 'awk \'$1 == "{wildcards.family}"\' {input} > {output}'

rule whatshap:
	input:
		vcfs_10X=lambda wildcards: ['10X/{sample}/filtered/{chromosome}.vcf.gz'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		vcf_consensus='consensus/freebayes_10X/{family}.{chromosome}.vcf',
		ped='ped/{family}.ped'
	output:
		vcf='whatshap/10X/{family}.{chromosome}.vcf'
	log: 'whatshap/10X/{family}.{chromosome}.log'
	shell:
		'~/scm/whatshap/bin/whatshap phase --ped {input.ped} -o {output.vcf} {input.vcf_consensus} {input.vcfs_10X} > {log} 2>&1'


rule compare_phasing:
	input:
		pedmec10X=lambda wildcards: 'whatshap/10X/{family}.{chromosome}.vcf'.format(family=sample2fam[wildcards.sample],chromosome=wildcards.chromosome),
		strandseq='input/strandseq/{sample}/{chromosome}.vcf'
	output: 'compare/pedmec10X-strandseq/{sample}/{chromosome}.txt'
	log: 'compare/pedmec10X-strandseq/{sample}/{chromosome}.log'
	shell:
		'~/scm/whatshap/bin/whatshap compare --names pedmec10X,strandseq {input.pedmec10X} {input.strandseq} > {output} 2> {log}'
