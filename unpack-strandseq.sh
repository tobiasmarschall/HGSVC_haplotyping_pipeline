#!/bin/bash
set -euo pipefail

unzip $1
for f in *.zip ; do 
	d=`echo $f | sed -r 's|VCFfiles_(.*)\.zip|\1|g'`
	echo $d
	mkdir $d
	cd $d
	unzip ../${f}
# 	for f2 in chr*_phased_${d}.vcf ; do
# 		sed -i "s|NA12878|${d}|g" ${f2}
# 	done
	mmv "VCFfiles_${d}/chr*_phased_${d}.vcf" 'chr#1.vcf'
	rmdir "VCFfiles_${d}"
	cd ..
	rm ${f}
done
