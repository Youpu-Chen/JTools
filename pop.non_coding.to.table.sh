#!/bin/bash
#SBATCH -p medium
#SBATCH -C scratch
#SBATCH -a 1-3
#SBATCH -n 20
#SBATCH -t 24:00:00

module load gatk

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
gatk VariantsToTable -V  CladeE.NMS.vcf.gz -F CHROM -F POS -F REF -F AN -F DP -GF GT -O CladeE.NMS.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then
gatk VariantsToTable -V  CladeW1.NMS.vcf.gz -F CHROM -F POS -F REF -F AN -F DP -GF GT -O CladeW1.NMS.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then
gatk VariantsToTable -V  CladeW2.NMS.vcf.gz -F CHROM -F POS -F REF -F AN -F DP -GF GT -O CladeW2.NMS.table
fi