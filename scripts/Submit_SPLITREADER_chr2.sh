#!/bin/bash

echo -e "\n\n"

workdir="/data/a2e/pbaduel/SPLITREADER"
project_dir='/workspaces/a2e/pbaduel/EPICTEA'

bamext='dupl_fixed_paired'

cd $workdir
if [ ! -d fastq ]; then
    mkdir ./fastq
fi
if [ ! -d BAMs ]; then
    mkdir ./BAMs
fi
if [ ! -d localtmp ]; then
    mkdir ./localtmp
fi

cd $workdir/Reference/
echo -en "reference genome...\t"
if [ ! -e TAIR10.fa.bwt ]; then
    ## build indexes for reference genome
    bowtie2-build TAIR10.fa TAIR10
    samtools faidx TAIR10.fa
    bwa index -a is TAIR10.fa
fi
echo -en "done\n"
cd $workdir/TE_sequence/
echo -en "TE library...\t"
if [ ! -e TE_all_Athaliana.fa.bwt ]; then
    ## build indexes for TE library
    bowtie2-build TE_all_Athaliana.fa TE_all_Athaliana
    bwa index -a is TE_all_Athaliana.fa
fi
echo -en "done\n"

echo -e "\n\n"

cd $workdir

setName='qlife'
if [ ! -d $setName ]; then
    mkdir ./$setName
fi

cd $workdir/fastq

fqNames=($(ls *.1.fastq | sed -e 's/\.Chr2.1.fastq$//'))
fqNB=${#fqNames[@]}
fqNB=1
echo -e "$fqNB samples: ${fqNames[*]}\n\n"

for (( fq=0; fq<$fqNB; fq++ ))  
    do
    popname=${fqNames[$fq]}
    echo -e -n "\t$popname...\t"
    
    echo "processing"
    cd $workdir/BAMs

    # align on reference genome
    bowtie2 --mp 13 --rdg 8,5 --rfg 8,5 -x $workdir/Reference/TAIR10 -1 $workdir/fastq/$popname.Chr2.1.fastq -2 $workdir/fastq/$popname.Chr2.2.fastq -S $workdir/BAMs/$popname.Chr2.sam -p 2

    # convert in bam, sort and index
    samtools view -bS $workdir/BAMs/$popname.Chr2.sam > $workdir/BAMs/$popname.Chr2.bam
    samtools sort -O bam -o $workdir/BAMs/$popname.Chr2.sort.bam $workdir/BAMs/$popname.Chr2.bam
    samtools index $workdir/BAMs/$popname.Chr2.sort.bam

    cd $workdir/$setName
    if [ ! -d $popname ]; then
        mkdir ./$popname
    fi
    cd $workdir/$setName/$popname
    if [ ! -d part1 ]; then
        mkdir ./part1
    fi
    cd $workdir/$setName/$popname/part1
    # run SPLITREADER part 1
    $workdir/scripts/SPLITREADER-beta1.5_part1.sh $popname part1 $workdir/BAMs .Chr2.sort $setName $workdir TE_all_Athaliana #> $popname.SR1.e

    
    cd $workdir/$setName/$popname
    if [ ! -d part2 ]; then
        mkdir ./part2
    fi
    cd $workdir/$setName/$popname/part2
    $workdir/scripts/SPLITREADER-beta1.5_part2.sh $popname $setName $workdir Reference/TAIR10 TAIR10_Quesneville_GFF3_transposons_only

done

echo -e "\n\n\n\n"
