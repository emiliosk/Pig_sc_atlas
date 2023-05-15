#!/bin/bash
# SNIC 2022/5-263

#SBATCH -A snic2022-5-263
#SBATCH -p node
#SBATCH -t 0-5:00:00
#SBATCH -n 1
#SBATCH -J mk_ref_STAR.sh
#SBATCH -e mk_ref_STAR.err
#SBATCH -o mk_ref_STAR.o
#SBATCH --mail-user emilio4ag@gmail.com
#SBATCH --mail-type=ALL

##Script used to setup the reference to be used with STARsolo

module load bioinfo-tools
module load star/2.7.9a

VERSION=109
ref_name=Sus_scrofa_genome_"$VERSION"_STAR
gtf_file=Sus_scrofa.Sscrofa11.1."$VERSION".gtf
fasta_file=Sus_scrofa.Sscrofa11.1.dna.toplevel.fa

cd /proj/rnaatlas/nobackup/private/EmilioTemp/ref/

STAR --runMode genomeGenerate \
 --genomeDir $ref_name \
 --genomeFastaFiles $fasta_file \
 --sjdbGTFfile $gtf_file \
 --sjdbOverhang 99
