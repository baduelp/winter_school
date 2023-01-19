#!/bin/bash 

#############################################
#                                           #  
#                                           #
#           metaTE beta1.5  part 2          #
#                                           #
#                                           # 
#############################################

##Questions or comments to quadrana(ar)biologie.ens.fr

#$1 -> basename of the bam file (i.e. if the file is "sequencing.bam", you should enter only "sequencing")

in=$1
cohort=$2
workdir=$3
GenomeFile=$4
TEannotfile=$5

#############################################################
#IN THE FOLLOWING VARIABLE YOU CAN PROVIDE THE MINIMUM READ LEGTH. By default this is 100nt
# LENGTH=100
LENGTH=50 # capture
#### If not specified, the program will calculate the longest read


#############################################################
#IN THE FOLLOWING VARIABLE YOU CAN PROVIDE THE MINIMUM NUMBER OF SPLIT-READS IN EACH EXTRIMITY OF THE TE. By default is 5 reads
READS=2
maxcov=100
# maxcov=1000000 # for capture
#LS=$2
LS=400
pe=TRUE
#### If not specified, the program estimate them. To this end, the program calculates the minimum number of reads as 3 standard deviation under the mean whole genome coverage
####This value should be at least 3, if not, it is forced to be 3
 
#############################################################

#############################################################
#IN THE FOLLOWING VARIABLE YOU CAN EXPLICITE THE NUMBER OF THREADS YOU WANT TO USE FOR MAPPING. By default is 2 threads
CORES=2
#############################################################


########################### edit the following paths !!! ####################

InputDir=/$workdir/$cohort/${in}/part1 #TE-sam from part1 

OutputDir=/$workdir/$cohort/${in}/part2
# if [ -e $OutputDir/log.txt ]
# then
#   rm $OutputDir/log.txt
# fi
# if [ -e $OutputDir/$in-insertion-sites.bed ]
# then
#   rm $OutputDir/$in-insertion-sites.bed
# fi
# echo '' > $OutputDir/$in-insertion-sites.bed
#path to list of TEs to analyze (TE-information-all.txt)
#TE-> should be indicated in the first column of the TE-information.txt file located in listDir
#TSD -> should be indicated in the second column of the TE-information.txt file located in listDir
listDir=/$workdir/TE_sequence
TE_info='TE-information-all'
# TE_info='TE-information-TECAP' #to match ASfilt2 TE capture list

TEannot=/$workdir/TE_sequence/$TEannotfile.gff
# TEannot=/groups/a2e/TAIR9/TAIR10_Quesneville_GFF3_transposons_only.gff

#Bowtie2 index for reference genome
GenomeIndexFile=/$workdir/$GenomeFile
# GenomeIndexFile=/workspaces/a2e/pbaduel/EPICTEA/Ath_reference/TAIR10

#Bowtie2 executable path
Bowtie2Dir=/usr/local/bin
# #samtools executable path
samtoolsDir=/usr/local/bin
#bedtools executable path
bedtoolsdir=/usr/local/bin
#picard tools executable path
# picardDir=/import/bc_users/bioinfo/vincens/pkg/te-tracker/picard-tools-1.119
# picardDir=/usr/share/java/picard.jar
picardDir=/import/localsrc/picard-2.26.2/bin/
# PID of this batch
IDPID=$$

# Temporary directory
TmpDir=/$workdir/localtmp/PB-$in
# TmpDir=/localtmp/PB-$IDPID
mkdir -p $TmpDir
##
if [ ! -d "$OutputDir" ]; then
  mkdir -p $OutputDir
fi


############################################################################



function cleanup {
  rm -f $TmpDir/*
  rm -r -f $TmpDir
  }
 trap cleanup EXIT
  
 
echo "["$(date +"%y-%m-%d %T")"] Running SPLITREADER beta1.5 part 2" | tee -a $OutputDir/log.txt
echo '' > $OutputDir/$in-insertion-sites.bed
if [ ! -e $InputDir/$in-TE.bam ]
then  
  if [ -e $InputDir/$in-TE.sam ]
  then
    echo "convert sam to bam" | tee -a $OutputDir/log.txt
     $samtoolsDir/samtools view -Sb $InputDir/$in-TE.sam > $InputDir/$in-TE.bam
     rm $InputDir/$in-TE.sam
  fi
fi
  
#Extracting superfamily reads from remapped TE sam
echo "["$(date +"%y-%m-%d %T")"] Extracting superfamily reads from ${in}" | tee -a $OutputDir/log.txt

end=`wc -l $listDir/$TE_info.txt | awk '{print $1}'`
# end=1
###Starting the SPLITREADER pipeline for each TE in the TE-information.txt file
for ((l=1; $l<=$end; l=$l+1)); do

  STARTTIME=$(date +%s)


  TE=`sed -n "${l}p" $listDir/$TE_info.txt | awk '{print $1}'`
  TSD=`sed -n "${l}p" $listDir/$TE_info.txt | awk '{print $2}'`
  #  TE=ATCOPIA78
  #  TSD=5
  echo -e "\n"
  if [ ! -s $OutputDir/$in-$TE-split.bam ]; then
  
    TmpResultsDir=$OutputDir/$TE
    mkdir -p $TmpResultsDir
    echo "["$(date +"%y-%m-%d %T")"] ##### RUNNING SPLIT-READ ANALYSIS ON $TE (TSD size = $TSD bp)######"  | tee -a $OutputDir/log.txt  
    echo ""
    ############# 

    # Selecting split-reads by mapping the unmapped reads over TE extremities
      
    #############################################################
      ###filter soft-clipped reads with at least 20nt softclipped at 5' or 3' read's end 
    echo "["$(date +"%y-%m-%d %T")"] Selecting split-reads"  | tee -a $OutputDir/log.txt    
    $samtoolsDir/samtools view -F 4 $InputDir/$in-TE.bam | grep -w $TE | awk '$6~/^[2-8][0-9]S/ || $6~/[2-8][0-9]S$/ || $6~/^1[0-9][0-9]S/ || $6~/1[0-9][0-9]S$/ {print $1"\t"$10"\t"$11}' | sort -k1,1 -k2,2 -u | awk '{print "@"$1"\n"$2"\n+\n"$3}' > $TmpResultsDir/$in-$TE-split.fastq 2>> $OutputDir/log.txt

    #############################################################
      ###filtering pair ends of unmapped reads  
    echo "["$(date +"%y-%m-%d %T")"] Selecting discordant-reads"  | tee -a $OutputDir/log.txt

    $samtoolsDir/samtools view -F 8 $InputDir/$in-TE.bam | grep -w $TE | awk '{print $1}' >> $TmpResultsDir/reads.name.disc  2>> $OutputDir/log.txt
    
    Nsplitreads=`wc -l $TmpResultsDir/$in-$TE-split.fastq | awk '{print $1/4}'` 
    Ndiscreads=`sort -k1,1 $TmpResultsDir/reads.name.disc | uniq | wc -l | awk '{print $1}'`  
    
    echo "$Nsplitreads splitreads and $Ndiscreads discordant reads identified on $TE" | tee -a $OutputDir/log.txt  

    if [ $(($Nsplitreads+$Ndiscreads)) -eq 0 ]
    then
      echo "No reads identified: end analysis on $TE" | tee -a $OutputDir/log.txt 

    else
      
      if [ $Ndiscreads -gt 0 ]
      then

        # grab first in pair
        $samtoolsDir/samtools view -f 64 -u $InputDir/$in-TE.bam > $TmpResultsDir/$in-TE-disc-first.bam 2>> $OutputDir/log.txt
        # grab second in pair
        $samtoolsDir/samtools view -f 128 -u $InputDir/$in-TE.bam > $TmpResultsDir/$in-TE-disc-second.bam 2>> $OutputDir/log.txt
        
        # java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp -cp $picardDir net.sf.picard.sam.FilterSamReads INPUT=$TmpResultsDir/$in-TE-disc-first.bam FILTER=includeReadList READ_LIST_FILE=$TmpResultsDir/reads.name.disc OUTPUT=$TmpResultsDir/$in-$TE-selected-disc-first.sam TMP_DIR=$TmpDir/javatemp  2>> $OutputDir/log.txt

        # # ## cannot specify memory requirement > use .jar instead
        picard-tools FilterSamReads INPUT=$TmpResultsDir/$in-TE-disc-first.bam FILTER=includeReadList READ_LIST_FILE=$TmpResultsDir/reads.name.disc OUTPUT=$TmpResultsDir/$in-$TE-selected-disc-first.sam TMP_DIR=$TmpDir/javatemp  2>> $OutputDir/log.txt
        
        # java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp -jar $picardDir/picard.jar FilterSamReads -INPUT $TmpResultsDir/$in-TE-disc-first.bam -FILTER includeReadList -READ_LIST_FILE $TmpResultsDir/reads.name.disc -OUTPUT $TmpResultsDir/$in-$TE-selected-disc-first.sam -TMP_DIR $TmpDir/javatemp  2>> $OutputDir/log.txt
    
        # java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp -cp $picardDir net.sf.picard.sam.FilterSamReads INPUT=$TmpResultsDir/$in-TE-disc-second.bam FILTER=includeReadList READ_LIST_FILE=$TmpResultsDir/reads.name.disc OUTPUT=$TmpResultsDir/$in-$TE-selected-disc-second.sam TMP_DIR=$TmpDir/javatemp  2>> $OutputDir/log.txt
              
        # # # ## cannot specify memory requirement > use .jar instead
        picard-tools FilterSamReads INPUT=$TmpResultsDir/$in-TE-disc-second.bam FILTER=includeReadList READ_LIST_FILE=$TmpResultsDir/reads.name.disc OUTPUT=$TmpResultsDir/$in-$TE-selected-disc-second.sam TMP_DIR=$TmpDir/javatemp  2>> $OutputDir/log.txt

        # java -Xmx10g -XX:+UseSerialGC -Djava.io.tmpdir=$TmpDir/javatemp -jar $picardDir/picard.jar FilterSamReads -INPUT $TmpResultsDir/$in-TE-disc-second.bam -FILTER includeReadList -READ_LIST_FILE $TmpResultsDir/reads.name.disc -OUTPUT $TmpResultsDir/$in-$TE-selected-disc-second.sam -TMP_DIR $TmpDir/javatemp  2>> $OutputDir/log.txt
        
        
        cat $TmpResultsDir/$in-$TE-selected-disc-first.sam | awk '$1!~/^@/ {print $1"\t"$10"\t"$11}' | sort -u -k1,1 -k2,2 |  awk '{print "@"$1"|1\n"$2"\n+\n"$3}' > $TmpResultsDir/$in-$TE-disc.fastq 2>> $OutputDir/log.txt
        
        cat $TmpResultsDir/$in-$TE-selected-disc-second.sam | awk '$1!~/^@/ {print $1"\t"$10"\t"$11}' | sort -u -k1,1 -k2,2 | awk '{print "@"$1"|2\n"$2"\n+\n"$3}' >> $TmpResultsDir/$in-$TE-disc.fastq 2>> $OutputDir/log.txt
      fi
      

      
      rm -f $TmpResultsDir/reads.name
        
      ################################
    
    
      ###Estimating max read size (If necessary)
      
      if [ -z "$LENGTH" ]
      then
        LENGTH=`head -5000 $TmpResultsDir/$in-$TE-split.fastq | awk 'NR%4 == 2 {print length($0)}' | sort | tail -1 `  
        length=$((LENGTH-20))
        echo "["$(date +"%y-%m-%d %T")"] Maximum Read length: $LENGTH [Estimated] " | tee -a $OutputDir/log.txt
        else
        length=$((LENGTH-20))
        echo "["$(date +"%y-%m-%d %T")"] Maximum Read length: $LENGTH [User defined] " | tee -a $OutputDir/log.txt
      fi
      
          
      ###Recursive split-reads mapping
      # step 1 for 3' read extremity: begining the loop.
      
        
      echo "["$(date +"%y-%m-%d %T")"] Analyzing split-reads" | tee -a $OutputDir/log.txt
      
      
      if [ $Nsplitreads -gt 0 ]
      then
        $Bowtie2Dir/bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$in-$TE-split.fastq -S $TmpResultsDir/$in-$TE-local.sam --local --very-sensitive --threads $CORES --quiet 

        $samtoolsDir/samtools view -H -S $TmpResultsDir/$in-$TE-local.sam > $TmpResultsDir/$in-$TE-split-local-up.sam 
        cat $TmpResultsDir/$in-$TE-split-local-up.sam > $TmpResultsDir/$in-$TE-split-local-down.sam 
        
        ##extracting split reads upstream insertion
        $samtoolsDir/samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-local.sam | awk '$6~/^[0-9][0-9]S/ || $6~/^1[0-9][0-9]S/ {print $0}' >> $TmpResultsDir/$in-$TE-split-local-down.sam 
        $samtoolsDir/samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-local.sam | awk '$6~/[0-9][0-9]S$/ || $6~/1[0-9][0-9]S$/ {print $0}' >> $TmpResultsDir/$in-$TE-split-local-up.sam 

        $samtoolsDir/samtools view -Sbu $TmpResultsDir/$in-$TE-split-local-down.sam | $samtoolsDir/samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-down.bam 
        $samtoolsDir/samtools view -Sbu $TmpResultsDir/$in-$TE-split-local-up.sam | $samtoolsDir/samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-up.bam 
        
        rm -f $TmpResultsDir/$in-$TE-split-local-up.sam
        rm -f $TmpResultsDir/$in-$TE-split-local-down.sam
          
        ############################################
        
        ##Refining insertion sites

        length=$(($((LENGTH/2))-1))
        echo "["$(date +"%y-%m-%d %T")"] Refining insertion sites" | tee -a $OutputDir/log.txt
        echo -e "\n"
        echo -n "Progression: ["
        
        $Bowtie2Dir/bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$in-$TE-split.fastq -S $TmpResultsDir/$in-$TE-splitjunction-5-$length.sam --un $TmpResultsDir/$in-$TE-split-5-$length -5 $length --very-sensitive --threads $CORES --local --quiet 2>> $OutputDir/log.txt
      
        $samtoolsDir/samtools view -H $TmpResultsDir/$in-$TE-splitjunction-5-$length.sam > $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam 2>> $OutputDir/log.txt
        cat $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam > $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.sam 2>> $OutputDir/log.txt
        cat $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam > $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.sam 2>> $OutputDir/log.txt
        cat $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam > $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.sam 2>> $OutputDir/log.txt
        
        $samtoolsDir/samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-splitjunction-5-$length.sam | awk '$6~/^[0-9][0-9]S/ || $6~/^1[0-9][0-9]S/ {print $0}' >> $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam 2>> $OutputDir/log.txt

        $samtoolsDir/samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-splitjunction-5-$length.sam | awk '$6~/[0-9][0-9]S$/ || $6~/^1[0-9][0-9]S/ {print $0}' >> $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.sam 2>> $OutputDir/log.txt

        $samtoolsDir/samtools view -Sbu $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam  | $samtoolsDir/samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.bam 2>> $OutputDir/log.txt
        $samtoolsDir/samtools view -Sbu $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.sam  | $samtoolsDir/samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.bam 2>> $OutputDir/log.txt

        rm -f $TmpResultsDir/$in-$TE-splitjunction-5-$length.sam
        rm -f $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.sam 
        rm -f $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.sam 
      
            
        $Bowtie2Dir/bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$in-$TE-split.fastq -S $TmpResultsDir/$in-$TE-splitjunction-3-$length.sam --un $TmpResultsDir/$in-$TE-split-3-$length -3 $length --very-sensitive --threads $CORES --local --quiet 2>> $OutputDir/log.txt

        $samtoolsDir/samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-splitjunction-3-$length.sam | awk '$6~/^[0-9][0-9]S/ || $6~/^1[0-9][0-9]S/ {print $0}' >> $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.sam 2>> $OutputDir/log.txt
        $samtoolsDir/samtools view -F 4 -q 5 -S $TmpResultsDir/$in-$TE-splitjunction-3-$length.sam | awk '$6~/[0-9][0-9]S$/ || $6~/^1[0-9][0-9]S/ {print $0}' >> $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.sam 2>> $OutputDir/log.txt

        $samtoolsDir/samtools view -Sbu $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.sam  | $samtoolsDir/samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.bam 2>> $OutputDir/log.txt
        $samtoolsDir/samtools view -Sbu $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.sam  | $samtoolsDir/samtools sort - -o $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.bam 2>> $OutputDir/log.txt
        
        rm -f $TmpResultsDir/$in-$TE-splitjunction-3-$length.sam
        rm -f $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.sam 
        rm -f $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.sam 
        
        echo -n "]"
        echo -e "\n"
      fi

      ###########################
      if [ -n "$pe" ]
      then
        if [ $Ndiscreads -gt 0 ]
        then
          echo "["$(date +"%y-%m-%d %T")"] Analyzing discordant reads" | tee -a $OutputDir/log.txt
          
          # increase mapping penalities for mismatches and indels to prevent mismapped discordant reads to break-down the clusters
          $Bowtie2Dir/bowtie2 -x $GenomeIndexFile -U $TmpResultsDir/$in-$TE-disc.fastq -S $TmpResultsDir/$in-$TE-endtoend.sam --very-sensitive --threads $CORES --quiet  --mp 13 --rdg 8,5 --rfg 8,5
        
          $samtoolsDir/samtools view -Sbu -f 16 -q 5 $TmpResultsDir/$in-$TE-endtoend.sam | $samtoolsDir/samtools sort - -o $TmpResultsDir/$in-$TE-disc-down.bam 
          $samtoolsDir/samtools view -Sbu -F 16 -q 5 $TmpResultsDir/$in-$TE-endtoend.sam | $samtoolsDir/samtools sort - -o $TmpResultsDir/$in-$TE-disc-up.bam 
        
          rm -f $TmpResultsDir/$in-$TE-local.sam
          rm -f $TmpResultsDir/$in-$TE-endtoend.sam
        fi
      fi


      ############################################
      # Post-treatment:
      
      echo "Merge the 5' and 3' clusters to create the downstream and upstream cluster" | tee -a $OutputDir/log.txt

      $samtoolsDir/samtools merge -f -u $TmpResultsDir/$in-$TE-up.bam $TmpResultsDir/$in-$TE-splitjunction-5-$length-up.bam $TmpResultsDir/$in-$TE-splitjunction-3-$length-up.bam $TmpResultsDir/$in-$TE-splitjunction-up.bam  2>> $OutputDir/log.txt

      $samtoolsDir/samtools merge -f -u $TmpResultsDir/$in-$TE-down.bam $TmpResultsDir/$in-$TE-splitjunction-5-$length-down.bam $TmpResultsDir/$in-$TE-splitjunction-3-$length-down.bam $TmpResultsDir/$in-$TE-splitjunction-down.bam  2>> $OutputDir/log.txt

      $samtoolsDir/samtools sort $TmpResultsDir/$in-$TE-down.bam -o $TmpResultsDir/$in-$TE-down 2>> $OutputDir/log.txt
      $samtoolsDir/samtools sort $TmpResultsDir/$in-$TE-up.bam -o $TmpResultsDir/$in-$TE-up 2>> $OutputDir/log.txt
    
    

      echo "Calculate the coverage over mapped regions - filter regions according to minimum and maximum read-depth" | tee -a $OutputDir/log.txt
      if [ -s $TmpResultsDir/$in-$TE-up.bam ]; then
        $samtoolsDir/samtools depth $TmpResultsDir/$in-$TE-up.bam | awk -v M=$maxcov '$3<(M) {print $1 "\t" $2 "\t"$2"\t"$3}' | sort -k 1,1 -k2,2n > $TmpResultsDir/$in-$TE-up.bed 
        if [ -s $TmpResultsDir/$in-$TE-up.bed ]; then
          $bedtoolsdir/mergeBed -i $TmpResultsDir/$in-$TE-up.bed  -c 4 -o max | awk -v l=$LENGTH '($3-$2)<=l {print $0}' > $TmpResultsDir/$in-$TE-up-merge.bed 2>> $OutputDir/log.txt
        else
          echo '' > $TmpResultsDir/$in-$TE-up-merge.bed
        fi
      else
        echo '' > $TmpResultsDir/$in-$TE-up-merge.bed
        echo '' > $TmpResultsDir/$in-$TE-up.bam
      fi
      
      if [ -s $TmpResultsDir/$in-$TE-down.bam ]; then
        $samtoolsDir/samtools depth $TmpResultsDir/$in-$TE-down.bam | awk -v M=$maxcov ' $3<(M) {print $1 "\t" $2 "\t"$2"\t"$3}' | sort -k 1,1 -k2,2n > $TmpResultsDir/$in-$TE-down.bed 
        if [ -s $TmpResultsDir/$in-$TE-down.bed ]; then
        $bedtoolsdir/mergeBed -i $TmpResultsDir/$in-$TE-down.bed -c 4 -o max | awk -v l=$LENGTH '($3-$2)<=l {print $0}' > $TmpResultsDir/$in-$TE-down-merge.bed 2>> $OutputDir/log.txt
        else
          echo '' > $TmpResultsDir/$in-$TE-down-merge.bed
        fi
      else
        echo '' > $TmpResultsDir/$in-$TE-down-merge.bed
        echo '' > $TmpResultsDir/$in-$TE-down.bam
      fi
      #grep -w $TE $SequencesDir/TAIR10_TE.bed | awk '{print $1"\t"$2-5000"\t"$3+5000}' > $TmpResultsDir/disc-excluding.tmp
    
      teid=`echo $TE | sed 's/\@/\t/' | awk '{print $1}'`
      grep -w $teid $TEannot | awk '{print $1"\t"$4"\t"$5}' > $TmpResultsDir/disc-excluding.tmp

      echo "merge cluster of discordant-reads" | tee -a $OutputDir/log.txt
      if [ -n "$pe" ]
      then
        echo "["$(date +"%y-%m-%d %T")"] Searching for discordant-reads clusters..." | tee -a $OutputDir/log.txt
        
        $samtoolsDir/samtools depth $TmpResultsDir/$in-$TE-disc-down.bam | awk -v M=$maxcov '$3<(M) {print $1 "\t" $2 "\t"$2 "\t"$3}' | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-down.bed 
        if [ -s $TmpResultsDir/$in-$TE-disc-down.bed ]; then
          $bedtoolsdir/mergeBed -i $TmpResultsDir/$in-$TE-disc-down.bed  -d $LS -c 4 -o max | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-down.tmp  2>> $OutputDir/log.txt
        else
          echo '' > $TmpResultsDir/$in-$TE-disc-down.tmp
        fi

        $samtoolsDir/samtools depth $TmpResultsDir/$in-$TE-disc-up.bam | awk -v M=$maxcov '$3<(M) {print $1 "\t" $2 "\t"$2 "\t"$3}' | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-up.bed 
        if [ -s $TmpResultsDir/$in-$TE-disc-up.bed ]; then
          $bedtoolsdir/mergeBed -i $TmpResultsDir/$in-$TE-disc-up.bed  -d $LS -c 4 -o max | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-up.tmp 2>> $OutputDir/log.txt
        else
          echo '' > $TmpResultsDir/$in-$TE-disc-up.tmp
        fi
      
        #filtering discordat-read clusters that overlap by a legth >than TSD+15bp
        $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-disc-up.tmp -b $TmpResultsDir/$in-$TE-disc-down.tmp -wa -wb | awk -v l=$LENGTH '($3-$6)>l {print $1"\t"$2"\t"$7}' >> $TmpResultsDir/disc-excluding.tmp  2>> $OutputDir/log.txt

        sort -o $TmpResultsDir/disc-excluding.tmp -k1,1 -k2,2n $TmpResultsDir/disc-excluding.tmp  2>> $OutputDir/log.txt

        # $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-disc-up.tmp -b $TmpResultsDir/$in-$TE-disc-down.tmp -wa -wb | awk -v l=$LENGTH '($3-$6)>l {print $1"\t"$2"\t"$7}' >> $TmpResultsDir/disc-excluding.tmp 2>> $OutputDir/log.txt
        
        $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-disc-up.tmp -b $TmpResultsDir/disc-excluding.tmp	 -sorted -v | awk '{print $1 "\t" $2 "\t"$3}' > $TmpResultsDir/$in-$TE-disc-up.tmp2 2>> $OutputDir/log.txt
        
        $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-disc-down.tmp -b $TmpResultsDir/disc-excluding.tmp -sorted -v | awk '{print $1 "\t" $2 "\t"$3}' > $TmpResultsDir/$in-$TE-disc-down.tmp2 2>> $OutputDir/log.txt
        
        $bedtoolsdir/intersectBed -b $TmpResultsDir/$in-$TE-disc-up.bam -a $TmpResultsDir/$in-$TE-disc-up.tmp2 -bed -c | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-up.bed 2>> $OutputDir/log.txt

        $bedtoolsdir/intersectBed -b $TmpResultsDir/$in-$TE-disc-down.bam -a $TmpResultsDir/$in-$TE-disc-down.tmp2 -bed -c  | sort -k1,1 -k2,2n > $TmpResultsDir/$in-$TE-disc-down.bed 2>> $OutputDir/log.txt
        
        rm -f $TmpResultsDir/$in-$TE-disc-down.tmp*
        rm -f $TmpResultsDir/$in-$TE-disc-up.tmp*
      fi 

      # echo "Defining TE insertions based on splitreads only: searching for overlapping clusters meeting the expected TSD size and number of supporting reads
      echo "["$(date +"%y-%m-%d %T")"] Defining insertions..." | tee -a $OutputDir/log.txt
      
      # # Taking into account Helitrons with TSD=0
      # if [ "$TSD" -eq "0" ]
      # then
      #   awk '{print $1"\t"$2"\t"$3+1"\t"$4}' $TmpResultsDir/$in-$TE-up-merge.bed > $TmpResultsDir/$in-$TE-up-merge1.bed  2>> $OutputDir/log.txt
      #   awk '{print $1"\t"$2-1"\t"$3"\t"$4}' $TmpResultsDir/$in-$TE-down-merge.bed > $TmpResultsDir/$in-$TE-down-merge1.bed 2>> $OutputDir/log.txt
      #         $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-up-merge1.bed -b $TmpResultsDir/$in-$TE-down-merge1.bed -wo | awk -v tsd=$TSD -v te=$TE '($6-$2)>10 && ($7-$3)>10 && $9>=tsd && $9<(tsd+15) {print $1 "\t" $6 "\t" $3 "\t"te"\t" ($9-1)"\t"$4"\t"$8}' | $bedtoolsdir/intersectBed -a stdin -b $TmpResultsDir/disc-excluding.tmp -v >> $TmpResultsDir/$in-insertion-sites.1 2>> $OutputDir/log.txt
    
      # else
        $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-up-merge.bed -b $TmpResultsDir/$in-$TE-down-merge.bed -wo | awk -v tsd=$TSD -v te=$TE '($6-$2)>10 && ($7-$3)>10 && $9>=tsd && $9<(tsd+15) {print $1 "\t" $6 "\t" $3 "\t"te"\t" ($9-1)"\t"$4"\t"$8}' | $bedtoolsdir/intersectBed -a stdin -b $TmpResultsDir/disc-excluding.tmp -v >> $TmpResultsDir/$in-insertion-sites.1 2>> $OutputDir/log.txt
      # fi
      
      
      if [ -s $TmpResultsDir/$in-insertion-sites.1 ]
      then
      awk -v l=$LS '{print $1"\t"$2"\t"$3+l"\t"$4}' $TmpResultsDir/$in-$TE-disc-up.bed | $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-insertion-sites.1 -b stdin -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11}' > $TmpResultsDir/$in-insertion-sites.2  2>> $OutputDir/log.txt
      awk -v l=$LS '{print $1"\t"$2"\t"$3+l"\t"$4}' $TmpResultsDir/$in-$TE-disc-up.bed | $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-insertion-sites.1 -b stdin -v | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t""0"}' >> $TmpResultsDir/$in-insertion-sites.2  2>> $OutputDir/log.txt
      
      awk -v l=$LS ' $2-l>0 {print $1"\t"$2-l"\t"$3"\t"$4}' $TmpResultsDir/$in-$TE-disc-down.bed | $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-insertion-sites.2 -b stdin -wa -wb | awk '($6+$7)>=r {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$12}' >> $OutputDir/$in-insertion-sites.bed  2>> $OutputDir/log.txt
      
      # -v r=$READS '($6+$7)>=r (non informative filter)
      awk -v l=$LS '$2-l>0 {print $1"\t"$2-l"\t"$3"\t"$4}' $TmpResultsDir/$in-$TE-disc-down.bed | $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-insertion-sites.2 -b stdin -v | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t""0"}' >> $OutputDir/$in-insertion-sites.bed  2>> $OutputDir/log.txt
      
      else
      echo "No split-reads clusters identified... skipping this step" | tee -a $OutputDir/log.txt
      fi
      
      ################################
        
      #Defining TE insertions based in discordant reads and supported for at least one split-read:
      if [ -n "$pe" ]
        then
        echo "selecting cluster of discordant reads that were not called by split-reads"  | tee -a $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        awk -v l=$LS '{print $1"\t"$2"\t"$3+l"\t"$4}' $TmpResultsDir/$in-$TE-disc-up.bed | $bedtoolsdir/intersectBed -b $TmpResultsDir/$in-insertion-sites.1 -a stdin -v > $TmpResultsDir/$in-$TE-disc-up-outer.bed 2>> $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        awk -v l=$LS '$2>l {print $1"\t"$2-l"\t"$3"\t"$4}' $TmpResultsDir/$in-$TE-disc-down.bed | $bedtoolsdir/intersectBed -b $TmpResultsDir/$in-insertion-sites.1 -a stdin -v > $TmpResultsDir/$in-$TE-disc-down-outer.bed 2>> $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        #selecting overlapping upstream and downstream clusters of discordant reads   
        $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-disc-up-outer.bed -b $TmpResultsDir/$in-$TE-disc-down-outer.bed -wa -wb | awk -v l=$LS '$3>l {print $1"\t"$2"\t"$7"\t"$4"\t"$8"\t"$3-l"\t"$6+l}' > $TmpResultsDir/$in-$TE-disc-cluster.0 2>> $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        awk '$6<=($7+8) {print $1"\t"$6-4"\t"$7+4"\t"$4"\t"$5"\t"$2"\t"$3}' $TmpResultsDir/$in-$TE-disc-cluster.0 > $TmpResultsDir/$in-$TE-disc-cluster.1 2>> $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        awk '$6>($7+8) {print $1"\t"$7-4"\t"$6+4"\t"$4"\t"$5"\t"$2"\t"$3}' $TmpResultsDir/$in-$TE-disc-cluster.0 >> $TmpResultsDir/$in-$TE-disc-cluster.1 2>> $OutputDir/log.txt
        echo -n "."  >> $OutputDir/log.txt
        if [[ -s $TmpResultsDir/$in-$TE-disc-cluster.1 ]] ; then
          echo -n "."  >> $OutputDir/log.txt
          $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-disc-cluster.1 -b $TmpResultsDir/$in-$TE-up-merge.bed -wa -wb | awk '{print $1"\t"$2"\t"$10"\t"$4"\t"$5"\t"$11"\t"$6"\t"$7}' | sort -k1,1 -k2,2n | $bedtoolsdir/mergeBed -i stdin -c 4,5,6,7,8 -o max,max,max,max,max > $TmpResultsDir/$in-$TE-disc-cluster.3 2>> $OutputDir/log.txt
          echo -n "."  >> $OutputDir/log.txt
          $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-disc-cluster.1 -b $TmpResultsDir/$in-$TE-up-merge.bed -v | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""0""\t"$6"\t"$7}'  >> $TmpResultsDir/$in-$TE-disc-cluster.3 2>> $OutputDir/log.txt
          echo -n "."  >> $OutputDir/log.txt
          $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-disc-cluster.3 -b $TmpResultsDir/$in-$TE-down-merge.bed -wa -wb | awk -v te=$TE '{print $1"\t"$10"\t"$3"\t"te"\t""NA""\t"$6"\t"$12"\t"$4"\t"$5"\t"$7"\t"$8}' | sort -k1,1 -k2,2n | $bedtoolsdir/mergeBed -i stdin -c 4,5,6,7,8,9,10,11 -o distinct,distinct,distinct,max,max,max,min,max >> $OutputDir/$in-insertion-sites.bed 2>> $OutputDir/log.txt
          echo -n "."  >> $OutputDir/log.txt
          $bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-disc-cluster.3 -b $TmpResultsDir/$in-$TE-down-merge.bed -v | awk -v te=$TE '{print $1"\t"$2"\t"$3"\t"te"\t""NA""\t"$6"\t""0""\t"$4"\t"$5"\t"$7"\t"$8}' | sort -k1,1 -k2,2n | $bedtoolsdir/mergeBed -i stdin -c 4,5,6,7,8,9,10,11 -o distinct,distinct,distinct,max,max,max,min,max >> $OutputDir/$in-insertion-sites.bed 2>> $OutputDir/log.txt
          echo -n "."  >> $OutputDir/log.txt
        else

          echo "No discording-reads clusters identified... skiping this step" | tee -a $OutputDir/log.txt
        fi
      fi
        
      ###only if excluding inner pericentromeres and donnors TEs:
      #awk '$1!="" {print $0}' $SequencesDir/$TE.txt > $TmpResultsDir/$TE-donnor.txt
      #$bedtoolsdir/intersectBed -a $TmpResultsDir/$in-$TE-up-merge.bed -b $TmpResultsDir/$in-$TE-down-merge.bed -wo | awk -v tsd=$TSD -v te=$TE '($6-$2)>10 && ($7-$3)>10 && $9>=tsd {print $1 "\t" $6 "\t" $3 "\t"te"\t" ($9-1)"\t"$4"\t"$8}' | intersectBed -a stdin -b $TmpResultsDir/$TE-donnor.txt -v | intersectBed -a stdin -b $Cent -v >> $OutputDir/$in-insertion-sites.bed

      INSERTIONS=`grep -w $TE $OutputDir/$in-insertion-sites.bed | wc -l | awk '{print $1}'` 
      
      echo "["$(date +"%y-%m-%d %T")"] Split-read analyis done: $INSERTIONS putative insertions identified..." | tee -a $OutputDir/log.txt

      
      # if [ $INSERTIONS -gt 0 ]
      # then
        ###merging bam files and moving them to the output folder 
        $samtoolsDir/samtools merge -f $TmpResultsDir/$in-$TE-split.bam $TmpResultsDir/$in-$TE-up.bam $TmpResultsDir/$in-$TE-down.bam $TmpResultsDir/$in-$TE-disc-down.bam $TmpResultsDir/$in-$TE-disc-up.bam 
        $samtoolsDir/samtools sort $TmpResultsDir/$in-$TE-split.bam -o $TmpResultsDir/$in-$TE-split-sorted.bam 
        $samtoolsDir/samtools view -H $TmpResultsDir/$in-$TE-split-sorted.bam | grep ^@SQ > $TmpResultsDir/$in-$TE-split-sorted.sam 
        $samtoolsDir/samtools view -H $TmpResultsDir/$in-$TE-split-sorted.bam | grep ^@PG | head -1 >> $TmpResultsDir/$in-$TE-split-sorted.sam 
        $samtoolsDir/samtools view $TmpResultsDir/$in-$TE-split-sorted.bam >> $TmpResultsDir/$in-$TE-split-sorted.sam 
        $samtoolsDir/samtools view -hb $TmpResultsDir/$in-$TE-split-sorted.sam > $OutputDir/$in-$TE-split.bam 
        $samtoolsDir/samtools index $OutputDir/$in-$TE-split.bam 
      # fi
      ENDTIME=$(date +%s)
      echo "["$(date +"%y-%m-%d %T")"] It takes $((ENDTIME-STARTTIME)) seconds to analyse $TE" | tee -a $OutputDir/log.txt
    fi
    rm -rf $TmpResultsDir/*
    rmdir $TmpResultsDir
  else
  
    echo -e "["$(date +"%y-%m-%d %T")"] ##### FOUND EXISTING SPLIT-BAM ON $TE (TSD size = $TSD bp)######"  | tee -a $OutputDir/log.txt  
    echo ""
  fi
done
echo "["$(date +"%y-%m-%d %T")"] SPLITREADER PART2 COMPLETED" | tee -a $OutputDir/log.txt
# cp $OutputDir/log.txt $OutputDir/
