#!/bin/bash
#SBATCH -p medium
#SBATCH -a 1-30
#SBATCH -C scratch
#SBATCH -n 4
#SBATCH --qos long
#SBATCH -t 120:00:00

module load anaconda3
source activate
conda activate test


if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run1

fi


if [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run2

fi


if [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_3

fi


if [ $SLURM_ARRAY_TASK_ID -eq 4 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_4

fi


if [ $SLURM_ARRAY_TASK_ID -eq 5 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_5

fi


if [ $SLURM_ARRAY_TASK_ID -eq 6 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_6

fi


if [ $SLURM_ARRAY_TASK_ID -eq 7 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_7

fi


if [ $SLURM_ARRAY_TASK_ID -eq 8 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_8

fi


if [ $SLURM_ARRAY_TASK_ID -eq 9 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_9

fi


if [ $SLURM_ARRAY_TASK_ID -eq 10 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_10

fi


if [ $SLURM_ARRAY_TASK_ID -eq 11 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_11

fi


if [ $SLURM_ARRAY_TASK_ID -eq 12 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_12

fi


if [ $SLURM_ARRAY_TASK_ID -eq 13 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_13

fi


if [ $SLURM_ARRAY_TASK_ID -eq 14 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_14

fi


if [ $SLURM_ARRAY_TASK_ID -eq 15 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_15

fi


if [ $SLURM_ARRAY_TASK_ID -eq 16 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_16

fi


if [ $SLURM_ARRAY_TASK_ID -eq 17 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_17

fi


if [ $SLURM_ARRAY_TASK_ID -eq 18 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_18

fi


if [ $SLURM_ARRAY_TASK_ID -eq 19 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_19

fi


if [ $SLURM_ARRAY_TASK_ID -eq 20 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_20

fi


if [ $SLURM_ARRAY_TASK_ID -eq 21 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_21

fi


if [ $SLURM_ARRAY_TASK_ID -eq 22 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_22

fi


if [ $SLURM_ARRAY_TASK_ID -eq 23 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_23

fi


if [ $SLURM_ARRAY_TASK_ID -eq 24 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_24

fi


if [ $SLURM_ARRAY_TASK_ID -eq 25 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_25

fi


if [ $SLURM_ARRAY_TASK_ID -eq 26 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_26

fi


if [ $SLURM_ARRAY_TASK_ID -eq 27 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_27

fi


if [ $SLURM_ARRAY_TASK_ID -eq 28 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_28

fi


if [ $SLURM_ARRAY_TASK_ID -eq 29 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_29

fi


if [ $SLURM_ARRAY_TASK_ID -eq 30 ]
then

treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -o treemix_results_m17_run30/migration_17_run_30

fi