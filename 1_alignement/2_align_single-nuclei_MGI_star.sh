#!/bin/bash
# SNIC 2022/5-263
#SBATCH -A snic2022-5-263
#SBATCH -p node
#SBATCH -t 1-12:00:00
#SBATCH -n 1
#SBATCH -J align_single-nuclei.sh
#SBATCH -e align_single-nuclei.err
#SBATCH -o align_single-nuclei.o
#SBATCH --mail-user emilio4ag@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load star/2.7.9a



##Links to access the libraries for each tissues
##NA, means there was no data openly accessile online, data shared privately

declare -A acc_link
acc_link[Area_Postrema_1]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417903/CNX0347134/CNR0426987/*
acc_link[Area_Postrema_2]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417904/CNX0347135/CNR0426988/*
acc_link[Area_Postrema_3]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417905/CNX0347136/CNR0426989/*
acc_link[Cerebellum_1]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417906/CNX0347137/CNR0426990/*
acc_link[Cerebellum_2]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417907/CNX0347138/CNR0426991/*
acc_link[Cerebellum_3]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417908/CNX0347139/CNR0426992/*
acc_link[Cerebellum_4]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417909/CNX0347140/CNR0426993/*
acc_link[Cerebellum_5]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417910/CNX0347141/CNR0426994/*
acc_link[Cerebellum_6]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417911/CNX0347142/CNR0426995/*
acc_link[Cerebellum_7]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417912/CNX0347143/CNR0426996/*
acc_link[Cerebellum_8]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417913/CNX0347144/CNR0426997/*
acc_link[Heart_1]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417914/CNX0347145/CNR0426998/*
acc_link[Heart_2]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417915/CNX0347146/CNR0426999/*
acc_link[Heart_3]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417916/CNX0347147/CNR0427000/*
acc_link[Heart_4]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417917/CNX0347148/CNR0427001/*
acc_link[Heart_5]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417918/CNX0347149/CNR0427002/*
acc_link[Kidney_1]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417919/CNX0347150/CNR0427003/*
acc_link[Kidney_2]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417920/CNX0347151/CNR0427004/*
acc_link[Kidney_3]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417921/CNX0347152/CNR0427005/*
acc_link[Kidney_4]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417922/CNX0347153/CNR0427006/*
acc_link[Liver_1]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417923/CNX0347154/CNR0427007/*
acc_link[Liver_2]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417924/CNX0347155/CNR0427008/*
acc_link[Liver_3]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417925/CNX0347156/CNR0427009/*
acc_link[Liver_4]=NA
acc_link[OVoLT_1]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417929/CNX0347157/CNR0427010/*
acc_link[OVoLT_2]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417930/CNX0347158/CNR0427011/*
acc_link[OVoLT_3]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417931/CNX0347159/CNR0427012/*
acc_link[OVoLT_4]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417932/CNX0347160/CNR0427013/*
acc_link[Retina_1]=NA
acc_link[Retina_2]=NA
acc_link[Spleen_1]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417933/CNX0347161/CNR0427014/*
acc_link[Spleen_2]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417934/CNX0347162/CNR0427015/*
acc_link[Spleen_3]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417935/CNX0347163/CNR0427016/*
acc_link[SubfomicalOrgan_1]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417936/CNX0347164/CNR0427017/*
acc_link[SubfomicalOrgan_2]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417937/CNX0347165/CNR0427018/*
acc_link[SubfomicalOrgan_3]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0002165/CNS0417938/CNX0347166/CNR0427019/*

acc_link[Lung_1]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0001486/CNS0305922/CNX0280404/CNR0347245/*
acc_link[Lung_2]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0001486/CNS0305923/CNX0280405/CNR0347246/*
acc_link[Lung_3]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0001486/CNS0305924/CNX0280406/CNR0347247/*
acc_link[Lung_4]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0001486/CNS0305925/CNX0280407/CNR0347248/*
acc_link[Lung_5]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0001486/CNS0305926/CNX0280408/CNR0347249/*
acc_link[Lung_6]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0001486/CNS0305929/CNX0280409/CNR0347250/*
acc_link[Lung_7]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0001486/CNS0305927/CNX0280410/CNR0347251/*
acc_link[Lung_8]=ftp://ftp.cngb.org/pub/CNSA/data3/CNP0001486/CNS0305928/CNX0280411/CNR0347252/*

reference_directory=/proj/rnaatlas/nobackup/private/EmilioTemp/ref/Sus_scrofa_genome_109_STAR
homedir=/proj/rnaatlas/nobackup/private/EmilioTemp/pig_sc_109


cd $homedir

if [[ ! -d single_nuclei ]]
then
        mkdir single_nuclei
fi
cd single_nuclei

for acc in "${!acc_link[@]}"
do
        cd $homedir"/single_nuclei"
  test=2
        echo *******************************************
        echo "**************"$acc"**********************"
        echo *******************************************
        if [[ ! -d $acc"_fastq" ]]
        then
                mkdir $acc"_fastq"
                cd $acc"_fastq"
                wget ${acc_link[$acc]} && echo "downloaded" $acc

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
        test=$(( $test + $add_by))
                  fi
                done
                rm hash.md5
  fi


        cd $homedir"/single_nuclei"

        #Allignment
        input_fastqs=$homedir"/single_nuclei/"${acc}"_fastq/"*.fq.gz
        set $input_fastqs

  if [[ ! -d $acc && $test -eq 2 ]]
        then
    STAR --genomeDir $reference_directory \
           --soloType CB_UMI_Simple \
           --readFilesCommand zcat \
           --readFilesIn "$2" "$1" \
           --soloCBwhitelist /proj/rnaatlas/nobackup/private/EmilioTemp/BC_whitelist.txt \
           --soloCBstart 1 \
     			 --soloCBlen 20 \
           --soloUMIstart 21 \
           --soloUMIlen 10 \
           --outFileNamePrefix $homedir"/single_nuclei/"$acc"/" \
     --limitOutSJcollapsed 2000000 \
     --runThreadN 20 \
     --soloCBmatchWLtype 1MM multi Nbase pseudocounts \
     --soloUMIfiltering MultiGeneUMI_CR \
     --soloUMIdedup 1MM_CR \
		  rm $homedir"/single_nuclei/"$acc"/Aligned.out.sam"
      gzip $homedir"/single_nuclei/"$acc"/Solo.out/Gene/filtered/"*
			gzip $homedir"/single_nuclei/"$acc"/Solo.out/Gene/raw/"*

   else
     echo did not run STAR solo align for $acc
   fi


done
