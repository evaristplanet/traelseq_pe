#!/bin/bash

#this script was run on hpc alfred (lvgsrv2.epfl.ch)

docker="singularity run docker://registry.c4science.ch/lvg/py3:v4"
paired="y"

########################################
### preprocess
########################################
#unzip
cd inputdata/fastq
gunzip *.gz
cd -

#put umi in read name
mkdir -p inputdata/fastq_flex1
for fq1 in inputdata/fastq/*R1.fq; do
    fq2="${fq1/_R1.fq/_R2.fq}"
    out=$(basename $fq1 _R1.fq)
    cmd="$docker flexbar -r $fq1 -p $fq2 -b barcode.fa --umi-tags -t inputdata/fastq_flex1/$out"
    sem --jobs 20 $cmd
done
sem --wait
rename _barcode_barcode_ _ inputdata/fastq_flex1/*.fastq
rename _1.fastq _R1.fq inputdata/fastq_flex1/*.fastq
rename _2.fastq _R2.fq inputdata/fastq_flex1/*.fastq

#trim 1, 2 or 3 Ts at the beginning
mkdir -p inputdata/fastq_flex2
for fq1 in inputdata/fastq_flex1/*R1.fq; do
    out=$(basename $fq1 _R1.fq)
    cmd="$docker flexbar -r $fq1 --htrim-left T --htrim-min-length 1 --htrim-max-length 3 -t inputdata/fastq_flex2/$out"
    sem --jobs 20 $cmd
done
sem --wait

#put fq1 and fq2 in one folder
rename .fastq _R1.fq inputdata/fastq_flex2/*.fastq
mv inputdata/fastq_flex2/*.fq inputdata/fastq_flex1/
rm -rf inputdata/fastq_flex2
rm -rf inputdata/fastq_flex1/*.log

#remove adapters
mkdir -p inputdata/fastq_flex2
for fq1 in inputdata/fastq_flex1/*R1.fq; do
    fq2="${fq1/_R1.fq/_R2.fq}"
    out=$(basename $fq1 _R1.fq)
    cmd="$docker flexbar -r $fq1 -p $fq2 --adapter-preset TruSeq -ap ON -t inputdata/fastq_flex2/$out"
    sem --jobs 20 $cmd
done
sem --wait

#put stuff in place
rename _1.fastq _R1.fq inputdata/fastq_flex2/*.fastq
rename _2.fastq _R2.fq inputdata/fastq_flex2/*.fastq
mv inputdata/fastq_flex2 inputdata/fastq_flex

for file in inputdata/fastq_flex/*.fq; do awk '{if ($1 ~ /^@/) {split($0,a,"_| "); print $1 "_" a[3], a[2]} else {print $0}}' $file > $file.cl; done
for file in inputdata/fastq_flex/*.fq; do mv $file.cl $file; done

########################################
### fastqc
########################################
mkdir -p fastqc logs/fastqc
ext=".fq"
for fq in inputdata/fastq_flex/*$ext; do
    sample=$(basename $fq); sample="${sample%$ext}"
    cmd="$docker fastqc -q --outdir fastqc $fq"
    cmd="$cmd 1> logs/fastqc/$sample.out 2>logs/fastqc/$sample.err"
    sem --jobs 20 $cmd; #run max 20 jobs. don't wait for them
done 

########################################
### map t2t
########################################
logs="logs/mapping_t2t"
mkdir -p mapping_t2t $logs 
ncpu=10
index="inputdata/t2t/GCA_009914755.4.chrNames"
fqext="_R1.fq"
for fq1 in inputdata/fastq_flex/*$fqext; do
    fq2="${fq1/_R1.fq/_R2.fq}"
    reads="-1 $fq1 -2 $fq2"
    cmd="bowtie2"
    sample=$(basename $fq1); sample="${sample%$fqext}"
    #cmd="$docker $cmd --sensitive-local -p $ncpu -x $index $reads"
    cmd="$docker $cmd --local -D 15 -R 2 -N 1 -L 25 -i S,1,0.75 -p $ncpu -x $index $reads"
    cmd="$cmd | $docker samtools view -bSh - > mapping_t2t/$sample.bam"
    if [ ! -s mapping_t2t/$sample.bam ]; then
        #run max 10 jobs (3 hisat +1 samtools = 4 cpu each)
        sem --jobs 9 $cmd 1>$logs/$sample.out 2>$logs/$sample.err & 
    fi
done
sem --wait

tail -n 15 logs/mapping_t2t/*.err > summary_alignment_t2t.txt

########################################
### sort and index bams
########################################
ncpu="4"
mkdir -p tmp
for bam in mapping_t2t/*.bam; do
    sample=$(basename ${bam%.bam})
    cmd1="$docker samtools sort -T tmp -@ $ncpu -m 1G $bam > tmp/${sample}_s.bam"
    cmd2="mv -v tmp/${sample}_s.bam $bam"
    cmd3="$docker samtools index $bam"
    cmd="$cmd1; $cmd2; $cmd3"
    sem --jobs 5 $cmd & #run max 10 jobs (4 cpu each)
done

########################################
### use umi to remove duplicated reads t2t
########################################
conda activate umi_tools
logs="logs/umi_tools_t2t"
mkdir -p mapping_dedup_t2t $logs 
for bam in mapping_t2t/*.bam; do
    sample=$(basename $bam .bam)
    cmd="umi_tools dedup -I $bam -L $logs/${sample}.log -S mapping_dedup_t2t/${sample}.bam --paired --extract-umi-method read_id --method unique"
    cmd="$cmd; $docker samtools index mapping_dedup_t2t/${sample}.bam"
    if [ ! -s mapping_dedup_t2t/$sample.bam ]; then
        sem --id umi_tools --jobs 9 $cmd 1>$logs/$sample.out 2>$logs/$sample.err & 
    fi
done
conda deactivate

########################################
### Peak calling t2t
########################################
chips=("" "TRAEL_SCR_S1" "TRAEL_v6_S2")

logs="logs/macs2_t2t"
mkdir -p $logs MACS_t2t
ncomp=$((${#chips[@]}-1))
for i in `seq 1 $ncomp`; do
    chipbam=mapping_dedup_t2t/${chips[$i]}.bam
    sample=$(basename ${chipbam%.bam})
    if [ $paired = "n" ]; then
        name=${chips[$i]}
        cmd="$docker macs2 callpeak -t $chipbam -f BAM -g hs -n MACS_t2t/$name -B -q 0.01"
    elif [ $paired = "y" ]; then
        name=${chips[$i]}
        cmd="$docker macs2 callpeak -t $chipbam -f BAM -g hs -n MACS_t2t/$name -B -q 0.01 --format BAMPE"
    fi
    sem --id macs2 --jobs 20 $cmd 1>$logs/$sample.out 2>$logs/$sample.err & #run max 5 jobs (4 cpu each)
done
sem --wait

#make bed files of peaks
mkdir -p peaks_t2t
find MACS_t2t/ -name '*peaks.bed' -print0 | xargs -0 -I{} cp -v {} peaks_t2t
for f in MACS_t2t/*.narrowPeak; do
  cut -f 1-6 $f > "peaks_t2t/$(basename $f .narrowPeak).bed"
  awk -v OFS='\t' '{if ($5 > 50) print $1,$2,$3,$4,$5,$6}' $f > "peaks_t2t/$(basename $f .narrowPeak)_stringent.bed"
done

wc -l peaks_t2t/*.bed > summary_peaks_t2t.txt

########################################
### Peak calling hg19
########################################
chips=("" "TRAEL_SCR_S1" "TRAEL_v6_S2")

logs="logs/macs2_hg19"
mkdir -p $logs MACS_hg19
ncomp=$((${#chips[@]}-1))
for i in `seq 1 $ncomp`; do
    chipbam=mapping_dedup_hg19/${chips[$i]}.bam
    sample=$(basename ${chipbam%.bam})
    if [ $paired = "n" ]; then
        name=${chips[$i]}
        cmd="$docker macs2 callpeak -t $chipbam -f BAM -g hs -n MACS_hg19/$name -B -q 0.01"
    elif [ $paired = "y" ]; then
        name=${chips[$i]}
        cmd="$docker macs2 callpeak -t $chipbam -f BAM -g hs -n MACS_hg19/$name -B -q 0.01 --format BAMPE"
    fi
    sem --id macs2 --jobs 20 $cmd 1>$logs/$sample.out 2>$logs/$sample.err & #run max 5 jobs (4 cpu each)
done
sem --wait

#make bed files of peaks
mkdir -p peaks_hg19
find MACS_hg19/ -name '*peaks.bed' -print0 | xargs -0 -I{} cp -v {} peaks_hg19
for f in MACS_hg19/*.narrowPeak; do
  cut -f 1-6 $f > "peaks_hg19/$(basename $f .narrowPeak).bed"
  awk -v OFS='\t' '{if ($5 > 50) print $1,$2,$3,$4,$5,$6}' $f > "peaks_hg19/$(basename $f .narrowPeak)_stringent.bed"
done

wc -l peaks_hg19/*.bed > summary_peaks_hg19.txt

########################################
### QC report
########################################
$docker multiqc -f .
rm -rf fastqc/*.zip

########################################
### clean up
########################################
rmdir tmp
cd ..
