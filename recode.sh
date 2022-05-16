#!/bin/bash
#SBATCH -p medium
#SBATCH -C scratch
#SBATCH -a 1-3
#SBATCH -n 1
#SBATCH -t 24:00:00
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
perl -e 'while (<>) { s/chr/scaffold\_/; print; }' /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/00.raw_data/CladeE.NMS.table > CladeE_raw.table

cat CladeE_raw.table | tail -n+2 > CladeE.table

python3 /scratch1/users/lhe/complex/complex/23.Sweep/software/ScanTools_ProtEvol/recode012.py -i /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/01.ScanTools_ProtEvol/CladeE.table -pop CladeE -o /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/01.ScanTools_ProtEvol/
fi


if [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then
perl -e 'while (<>) { s/chr/scaffold\_/; print; }' /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/00.raw_data/CladeW1.NMS.table > CladeW1_raw.table

cat CladeW1_raw.table | tail -n+2 > CladeW1.table

python3 /scratch1/users/lhe/complex/complex/23.Sweep/software/ScanTools_ProtEvol/recode012.py -i /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/01.ScanTools_ProtEvol/CladeW1.table -pop CladeW1 -o /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/01.ScanTools_ProtEvol/
fi


if [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then
perl -e  'while (<>) { s/chr/scaffold\_/; print; }' /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/00.raw_data/CladeW2.NMS.table > CladeW2_raw.table

cat CladeW2_raw.table | tail -n+2 > CladeW2.table

python3 /scratch1/users/lhe/complex/complex/23.Sweep/software/ScanTools_ProtEvol/recode012.py -i /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/01.ScanTools_ProtEvol/CladeW2.table -pop CladeW2 -o /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/01.ScanTools_ProtEvol/
fi

rm -rf *table