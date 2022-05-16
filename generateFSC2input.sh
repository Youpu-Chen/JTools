#!/bin/bash
#SBATCH -p medium
#SBATCH -C scratch
#SBATCH -n 1
#SBATCH -t 24:00:00
python3 /scratch/users/lhe/complex/complex/23.Sweep/software/ScanTools_ProtEvol/FSC2input.py -i /scratch/users/lhe/complex/complex/23.Sweep/03.fastsimcoal2/pure_individual_non_coding/01.ScanTools_ProtEvol/E.W1.W2.recode.concat.txt -o ./ -prefix E.W1.W2 -ws 10000 -bs 0 -np 3 -alpha true