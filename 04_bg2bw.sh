#!/bin/bash

#get bw of all bg
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
#chmod 775 bedGraphToBigWig

#use it do the conversion
mkdir -p ../data/bigWig
for bg in ../data/bedGraph/t2t_p*.bg; do 
    bw="../data/bigWig/$(basename $bg bg)bw"
    cmd1="sort -k1,1 -k2,2n $bg | grep -v NA > /tmp/tmp.bg"
    cmd2="./bedGraphToBigWig /tmp/tmp.bg ~/samba/resources/data/genome/chromsizes/chm13v2.genome $bw"
    eval "$cmd1; $cmd2"
done
rm /tmp/tmp.bg
