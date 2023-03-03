#docker
#singularity run docker://registry.c4science.ch/lvg/r-4.2.1:v2

#load libs and funs
library(colorout)
source('funs.R')

#for t2t 10-250
file.fwd <- '../data/bedGraph/t2t_10_250_TRAEL_SCR_S1_fwd.bg'
file.rev <- '../data/bedGraph/t2t_10_250_TRAEL_SCR_S1_rev.bg'
cross.out <- '../data/cross/t2t_10_250_SCR_S1.bed'
cross.out.negpos <- '../data/cross_negpos/t2t_10_250_SCR_S1.bed'
rdata.out <- '../data/RData/t2t_10_250_SCR_S1.RData'
preprocess(file.fwd, file.rev, cross.out, cross.out.negpos, rdata.out)
file.fwd <- '../data/bedGraph/t2t_10_250_TRAEL_v6_S2_fwd.bg'
file.rev <- '../data/bedGraph/t2t_10_250_TRAEL_v6_S2_rev.bg'
cross.out <- '../data/cross/t2t_10_250_v6_S2.bed'
cross.out.negpos <- '../data/cross_negpos/t2t_10_250_v6_S2.bed'
rdata.out <- '../data/RData/t2t_10_250_v6_S2.RData'
preprocess(file.fwd, file.rev, cross.out, cross.out.negpos, rdata.out)
