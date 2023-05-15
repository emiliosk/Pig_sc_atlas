#!/bin/bash
# SNIC 2022/5-263
#SBATCH -A snic2022-5-263
#SBATCH -p node
#SBATCH -t 0-4:00:00
#SBATCH -n 1
#SBATCH -J mk_ref_cellranger.sh
#SBATCH -e mk_ref_cellranger.err
#SBATCH -o mk_ref_cellranger.o
#SBATCH --mail-user emilio4ag@gmail.com
#SBATCH --mail-type=ALL

#### Script for setting up the pig reference to use with cellranger

module load bioinfo-tools
module load cellranger/6.1.2

homedir=/proj/rnaatlas/nobackup/private/EmilioTemp
VERSION=109

ref_name=Sus_scrofa_genome_"$VERSION"_cellranger
gtf_file=Sus_scrofa.Sscrofa11.1."$VERSION".gtf
fasta_file=Sus_scrofa.Sscrofa11.1.dna.toplevel.fa

#Download and make cellranger reference

cd $homedir

if [[ ! -d ref ]]
then
  mkdir ref
fi

cd ref

#download gtf
if [[ ! -f Sus_scrofa.Sscrofa11.1."$VERSION".gtf ]]
then
  wget https://ftp.ensembl.org/pub/release-"$VERSION"/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1."$VERSION".gtf.gz
  gunzip Sus_scrofa.Sscrofa11.1."$VERSION".gtf.gz
fi

#download fasta
if [[ ! -f Sus_scrofa.Sscrofa11.1.dna.toplevel.fa ]]
then
  wget https://ftp.ensembl.org/pub/release-"$VERSION"/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
  gunzip Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
fi

#No filtering done at all on the gtfs
#make reference
if [[ ! -d $ref_name ]]
then
  cellranger mkref \
    --genome=$ref_name \
    --fasta=$fasta_file \
    --genes=$gtf_file\
    --ref-version=$VERSION
fi
