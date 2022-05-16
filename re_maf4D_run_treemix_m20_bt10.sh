#!/bin/bash
#SBATCH -p medium
#SBATCH -a 1-10
#SBATCH -C scratch
#SBATCH -n 4
#SBATCH --qos long
#SBATCH -t 120:00:00

module load anaconda3
source activate
conda activate test


if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_1

fi


if [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_2

fi


if [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_3

fi


if [ $SLURM_ARRAY_TASK_ID -eq 4 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_4

fi


if [ $SLURM_ARRAY_TASK_ID -eq 5 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_5

fi


if [ $SLURM_ARRAY_TASK_ID -eq 6 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_6

fi


if [ $SLURM_ARRAY_TASK_ID -eq 7 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_7

fi


if [ $SLURM_ARRAY_TASK_ID -eq 8 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_8

fi


if [ $SLURM_ARRAY_TASK_ID -eq 9 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_9

fi


if [ $SLURM_ARRAY_TASK_ID -eq 10 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_10

fi