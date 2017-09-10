#!/bin/bash

zcat consensus-merged/freebayes_10X/Y117.ad.vcf.gz | awk '($0 ~/^#/) || ((substr($10,1,3)=="0/0" && substr($11,1,3)=="0/0" && substr($12,1,3)=="0/1") && ($1!="chrX"))' | bgzip > Y117.denovo-snv.vcf.gz
tabix Y117.denovo-snv.vcf.gz
rm -rf tmp-denovo-snv
bcftools isec whatshap-merged-gt/consensus/Y117.pacbioblasr.single.vcf.gz Y117.denovo-snv.vcf.gz -p tmp-denovo-snv -n =2 -w 1
cat tmp-denovo-snv/0000.vcf | awk '($0 ~/^#/) || (substr($10,1,3)=="0/0" && substr($11,1,3)=="0/0" && (substr($12,1,3)=="0/1" || substr($12,1,3)=="0|1" || substr($12,1,3)=="1|0"))' | bgzip > Y117.denovo-snv.ill-10X-PB.vcf.gz 
