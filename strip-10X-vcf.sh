#!/bin/bash
bcftools view -h $1
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t.\tGT:PS[\t%GT:%PS]\n' $1
