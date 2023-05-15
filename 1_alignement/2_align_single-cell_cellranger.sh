#!/bin/bash
# SNIC 2022/5-263
#SBATCH -A snic2022-5-263
#SBATCH -p node
#SBATCH -t 7-0:00:00
#SBATCH -n 1
#SBATCH -J align_single-cell.sh
#SBATCH -e align_single-cell.err
#SBATCH -o align_single-cell.o
#SBATCH --mail-user emilio4ag@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load cellranger/6.1.2


VERSION=109
reference_directory=/proj/rnaatlas/nobackup/private/EmilioTemp/ref/Sus_scrofa_genome_"$VERSION"_cellranger
homedir=/proj/rnaatlas/nobackup/private/EmilioTemp/pig_sc_109

declare -A acc_link
acc_link[Adipose-S]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417941/CNX0347167/CNR0427020/*
acc_link[Adipose-V]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417942/CNX0347168/CNR0427021/*
acc_link[brain]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417943/CNX0347169/CNR0427022/*
acc_link[Intestine]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417944/CNX0347170/CNR0427023/*
acc_link[liver]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417945/CNX0347171/CNR0427024/*
acc_link[Lung]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417946/CNX0347172/CNR0427025/*
acc_link[PBMC]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417947/CNX0347173/CNR0427026/*
acc_link[retina]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417948/CNX0347174/CNR0427027/*
acc_link[spleen]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417949/CNX0347175/CNR0427028/*


cd $homedir

if [[ ! -d single_cell ]]
then
        mkdir single_cell
fi
cd single_cell

for acc in "${!acc_link[@]}"
do
        cd $homedir"/single_cell"
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
                  add_by=1
                  if [[ $og_hash == $new_hash ]]
                  then
                    echo $md5 "fastq correctly downloaded"
                    test=$(( $test + $add_by ))
                  fi
                done
                rm hash.md5
        fi
  cd $homedir"/single_cell"
  #Allignment 
  if [[ ! -d $acc && $test -eq 2 ]]
  then
    cellranger count \
      --id=$acc \
      --transcriptome=$reference_directory \
      --fastqs="$homedir"/single_cell/"$acc"_fastq \
      --sample=$acc \
      --no-bam
  fi
done
