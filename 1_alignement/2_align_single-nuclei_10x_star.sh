#!/bin/bash
# SNIC 2022/5-263
#SBATCH -A snic2022-5-263
#SBATCH -p node
#SBATCH -t 1-12:00:00
#SBATCH -n 1
#SBATCH -J align_single-nuclei_10x_star.sh
#SBATCH -e align_single-nuclei_10x_star.err
#SBATCH -o align_single-nuclei_10x_star.o
#SBATCH --mail-user emilio4ag@gmail.com
#SBATCH --mail-type=ALL
 
module load bioinfo-tools
module load star/2.7.9a


declare -A acc_link
acc_link[frontal]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0000686/CNS0127223/CNX0111986/CNR0134259/*
acc_link[hypothalamus1]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0000686/CNS0127224/CNX0111987/CNR0134260/*
acc_link[hypothalamus2]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0000686/CNS0127225/CNX0111988/CNR0134261/*
acc_link[occipital]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0000686/CNS0127226/CNX0111989/CNR0134262/*
acc_link[parietal]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0000686/CNS0127227/CNX0111990/CNR0134263/*
acc_link[temporal]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0000686/CNS0127228/CNX0111991/CNR0134264/*

declare -A acc_sample
acc_sample[frontal]=V300015611B
acc_sample[hypothalamus1]=CL100128991
acc_sample[hypothalamus2]=CL100128991
acc_sample[occipital]=CL100132063
acc_sample[parietal]=CL100132063
acc_sample[temporal]=V300015611B

VERSION=109
reference_directory=/proj/rnaatlas/nobackup/private/EmilioTemp/ref/Sus_scrofa_genome_109_STAR
homedir=/proj/rnaatlas/nobackup/private/EmilioTemp/pig_sc_109


cd $homedir

if [[ ! -d single_nuclei_10x ]]
then
        mkdir single_nuclei_10x
fi

for acc in "${!acc_link[@]}"
do
        cd $homedir"/single_nuclei_10x"
        test=2
        echo "*******************************************"
        echo "**************"$acc"**********************"
        echo "*******************************************"
        if [[ ! -d $acc"_fastq" ]]
        then
                mkdir $acc"_fastq"
                cd $acc"_fastq"
                wget -q ${acc_link[$acc]} && echo "downloaded" $acc

##Check .md5 fiels for integrity
                echo "md5 check"
                test=0
                add_by=1
                for md5 in *.md5
                do
                  echo $md5
                  temp_var=`cat $md5`
                  set -- $temp_var
                  og_hash=$2
                  file_name=$1
                  md5sum $file_name > hash.md5
                  temp_var=`cat hash.md5`
                  set -- $temp_var
                  new_hash=$1
                  if [[ $og_hash == $new_hash ]]
                  then
                        echo $md5 "fastaq correctly downloaded"
                        test=$(( $test + $add_by ))
                  fi
                done
                rm hash.md5
                for i in *
                do
                        if [[ ! "$i" == *S1_L00* ]]
                        then
                                temp_name="${i/L0/S1_L00}"
                                temp_name="${temp_name/read_/R}"
                                temp_name="${temp_name/.fq.gz/_001.fastq.gz}"
                                mv $i $temp_name
                        fi
                done
        fi

  cd $homedir"/single_nuclei_10x_star"
  echo "*******************************************"
  echo "**************"$acc"**********************"
  echo "*******************************************"

  input_fastqs="/proj/rnaatlas/private/pig/pig_sc/single_nuclei_10X_v2/"${acc}"_fastq/"*.fastq.gz
  set $input_fastqs


  if [[ ! -d $acc && $test -eq 2 ]]
  then
    STAR --genomeDir $reference_directory \
               --soloType CB_UMI_Simple \
               --readFilesCommand zcat \
               --readFilesIn "$2" "$1" \
               --soloCBwhitelist /proj/rnaatlas/nobackup/private/EmilioTemp/whitelist_10x.txt \
               --soloCBstart 1 \
               --soloCBlen 16 \
               --soloUMIstart 17 \
               --soloUMIlen 10 \
               --outFileNamePrefix $homedir"/single_nuclei_10x_star/"$acc"/" \
         --limitOutSJcollapsed 2000000 \
         --runThreadN 20 \
         --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
         --soloUMIfiltering MultiGeneUMI_CR \
         --soloUMIdedup 1MM_CR
         rm $homedir"/single_nuclei_10x_star/"$acc"/Aligned.out.sam"
         gzip $homedir"/single_nuclei_10x_star/"$acc"/Solo.out/Gene/filtered/"*
         gzip $homedir"/single_nuclei_10x_star/"$acc"/Solo.out/Gene/raw/"*
       else
     echo did not run STAR solo align for $acc
   fi
done
