#!/bin/bash
#SBATCH -p medium
#SBATCH -a 1-3
#SBATCH -C scratch
#SBATCH -n 20
#SBATCH -t 12:00:00
module load gatk

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
gatk --java-options "-Xmx6G" SelectVariants -R /scratch/users/lhe/complex/genome/brachista.fasta -V nCDS.nMis.passGTDP.DPINFO.BIA.passVSVchr01_19ex15.vcf.gz -sn new.CladeE.args -O /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/00.raw_data/CladeE.NMS.vcf.gz
fi

if [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then
gatk --java-options "-Xmx6G" SelectVariants -R /scratch/users/lhe/complex/genome/brachista.fasta -V nCDS.nMis.passGTDP.DPINFO.BIA.passVSVchr01_19ex15.vcf.gz -sn new.CladeW1.args -O /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/00.raw_data/CladeW1.NMS.vcf.gz
fi

if [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then
gatk --java-options "-Xmx6G" SelectVariants -R /scratch/users/lhe/complex/genome/brachista.fasta -V nCDS.nMis.passGTDP.DPINFO.BIA.passVSVchr01_19ex15.vcf.gz -sn CladeW2.args -O /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/00.raw_data/CladeW2.NMS.vcf.gz
fi