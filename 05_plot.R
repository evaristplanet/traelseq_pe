#docker
#singularity run docker://registry.c4science.ch/lvg/r-4.2.1:v2

#load libs and funs
library(colorout)
source('~/scripts/R/common/commonFuns.R')
source('funs.R')
mydircreate('../results/ReadPolarity')

#plot
chrom <- 'chr1'
pos.st <- 1.95e+08
pos.en <- 2.3e+08
xlab <- paste0('Chromosome ', chrom, ' position (Mb)')
pos <- paste0(chrom, '-', format(pos.st, scientific=F), '-', format(pos.en, scientific=F))

{
    #t2t 10-250
    rdata <- '../data/RData/t2t_10_250_SCR_S1.RData'
    fileout <- paste0('../results/ReadPolarity/t2t_10_250_SCR_S1_', pos, '.pdf')
    main <- 'SCR_S1'
    sexyPlot(chrom, post.st, pos.en, rdata, xlab, fileout, main)
    rdata <- '../data/RData/t2t_10_250_v6_S2.RData'
    fileout <- paste0('../results/ReadPolarity/t2t_10_250_v6_S2_', pos, '.pdf')
    main <- 'v6_S2'
    sexyPlot(chrom, post.st, pos.en, rdata, xlab, fileout, main)
}
