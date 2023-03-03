#!/bin/bash

#Description: get the crossings that are specific from s2

docker="singularity run docker://registry.c4science.ch/lvg/py3:v5"
cd ../data/cross_negpos
genome="/data/samba/resources/data/genome/chromsizes/chm13v2.genome"

#for t2t_1_25
$docker bedtools slop -i t2t_1_25_SCR_S1.bed -g $genome -b 125000 > tmp.bed
$docker bedtools intersect -v -a t2t_1_25_v6_S2.bed -b tmp.bed > t2t_1_25_v6_S2_notInS1.bed
rm tmp.bed
