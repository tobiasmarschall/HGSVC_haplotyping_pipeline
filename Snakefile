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
ref = 'ref/GRCh38_full_analysis_set_plus_decoy_hla.fa'
chromosomes = ['chr'+str(x) for x in range(1,23)] + ['chrX']

rule master:
	input: 
		expand('freebayes/SH032.{chromosome}.vcf', chromosome=chromosomes),
		expand('freebayes/Y117.{chromosome}.vcf', chromosome=chromosomes),
		expand('freebayes/PR05.{chromosome}.vcf', chromosome=chromosomes),
		expand('whatshap/freebayes_10X/{family}.chr22.vcf', family=families)

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
		vcf='freebayes/{family}.{chromosome}.vcf'
	shell:
		'samtools merge -R {wildcards.chromosome} - {input.crams} | freebayes -f {ref} --stdin > {output.vcf}'


rule extract_chromosome_10X:
	input: 
		vcf='ftp/working/20160513_10XGenomics_data_and_calls/{sample}_WGS/{sample}_WGS_phased_variants.vcf.gz'
	output:
		vcf='10X/{sample}/{chromosome}.vcf.gz'
	shell:
		'zcat {input.vcf} | awk \'($0 ~ /^#/) || ($1=="{wildcards.chromosome}")\' | gzip > {output.vcf}'

rule whatshap:
	input:
		vcfs_10X=lambda wildcards: ['10X/{sample}/{chromosome}.vcf.gz'.format(sample=s, chromosome=wildcards.chromosome) for s in samples[wildcards.family]],
		vcf_freebayes='freebayes/{family}.{chromosome}.vcf'
	output:
		vcf='whatshap/freebayes_10X/{family}.{chromosome}.vcf'
	log: 'whatshap/freebayes_10X/{family}.{chromosome}.log'
	shell:
		'~/scm/whatshap/bin/whatshap --ped ../1kg-sv-trios/chs.ped --genmap ~/scm/whatshap/tests/data/trio.map -o {output.vcf} {input.vcf_freebayes} {input.vcfs_10X} > {log} 2>{log}'
