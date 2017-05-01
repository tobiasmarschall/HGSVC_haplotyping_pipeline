#!/usr/bin/env Rscript

d = read.table('HG00732.het-rates.tsv', col.names=c('chromosome','start','end','hets','homs'))

pdf('HG00732.het-rates.pdf')
hist(d$hets, breaks = 400, xlab='HET SNVs per 100kb', xlim=c(0,400), main='Histogram of HET rates (HG00732)')
abline(v=10.73, col='red', lwd = 2)
dev.off()
