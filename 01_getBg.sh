#!/bin/bash

#use bamCoverage to get counts of reads in sliding window, once with first in strand being fwd and another time being rev
#bams are in s3

docker="singularity run docker://registry.c4science.ch/lvg/py3:v5"

#for t2t 10-250
cd ~/samba/dropfolder/evarist/filipe/2212_traelseq/mapping_dedup_t2t
$docker bamCoverage --bam TRAEL_SCR_S1.bam -o t2t_10_250_TRAEL_SCR_S1_fwd.bg --binSize 10000 --smoothLength 250000 --outFileFormat bedgraph --samFlagInclude 64 --samFlagExclude 16
$docker bamCoverage --bam TRAEL_SCR_S1.bam -o t2t_10_250_TRAEL_SCR_S1_rev.bg --binSize 10000 --smoothLength 250000 --outFileFormat bedgraph --samFlagInclude 80
$docker bamCoverage --bam TRAEL_v6_S2.bam -o t2t_10_250_TRAEL_v6_S2_fwd.bg --binSize 10000 --smoothLength 250000 --outFileFormat bedgraph --samFlagInclude 64 --samFlagExclude 16
$docker bamCoverage --bam TRAEL_v6_S2.bam -o t2t_10_250_TRAEL_v6_S2_rev.bg --binSize 10000 --smoothLength 250000 --outFileFormat bedgraph --samFlagInclude 80
cd -
mv ~/samba/dropfolder/evarist/filipe/2212_traelseq/mapping_dedup_t2t/*.bg ../data/bedGraph/
