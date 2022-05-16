#!/bin/bash
#SBATCH -p medium
#SBATCH -C scratch
#SBATCH -a 1-210
#SBATCH -n 1
#SBATCH --qos long
#SBATCH -t 120:00:00
module load anaconda3
source activate
conda activate test
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 0 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_0_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 0 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_0_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 0 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_0_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 4 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 0 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_0_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 5 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 0 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_0_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 6 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 0 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_0_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 7 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 0 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_0_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 8 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 0 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_0_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 9 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 0 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_0_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 10 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 0 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_0_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 11 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 1 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_1_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 12 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 1 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_1_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 13 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 1 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_1_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 14 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 1 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_1_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 15 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 1 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_1_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 16 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 1 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_1_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 17 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 1 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_1_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 18 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 1 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_1_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 19 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 1 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_1_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 20 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 1 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_1_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 21 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 2 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_2_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 22 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 2 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_2_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 23 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 2 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_2_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 24 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 2 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_2_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 25 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 2 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_2_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 26 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 2 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_2_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 27 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 2 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_2_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 28 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 2 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_2_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 29 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 2 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_2_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 30 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 2 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_2_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 31 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 3 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_3_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 32 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 3 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_3_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 33 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 3 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_3_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 34 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 3 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_3_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 35 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 3 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_3_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 36 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 3 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_3_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 37 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 3 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_3_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 38 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 3 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_3_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 39 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 3 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_3_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 40 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 3 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_3_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 41 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 4 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_4_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 42 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 4 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_4_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 43 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 4 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_4_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 44 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 4 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_4_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 45 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 4 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_4_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 46 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 4 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_4_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 47 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 4 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_4_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 48 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 4 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_4_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 49 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 4 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_4_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 50 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 4 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_4_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 51 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 5 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_5_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 52 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 5 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_5_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 53 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 5 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_5_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 54 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 5 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_5_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 55 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 5 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_5_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 56 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 5 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_5_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 57 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 5 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_5_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 58 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 5 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_5_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 59 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 5 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_5_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 60 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 5 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_5_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 61 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 6 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_6_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 62 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 6 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_6_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 63 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 6 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_6_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 64 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 6 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_6_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 65 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 6 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_6_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 66 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 6 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_6_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 67 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 6 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_6_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 68 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 6 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_6_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 69 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 6 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_6_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 70 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 6 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_6_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 71 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 7 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_7_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 72 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 7 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_7_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 73 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 7 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_7_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 74 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 7 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_7_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 75 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 7 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_7_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 76 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 7 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_7_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 77 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 7 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_7_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 78 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 7 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_7_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 79 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 7 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_7_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 80 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 7 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_7_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 81 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 8 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_8_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 82 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 8 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_8_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 83 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 8 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_8_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 84 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 8 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_8_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 85 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 8 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_8_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 86 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 8 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_8_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 87 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 8 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_8_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 88 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 8 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_8_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 89 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 8 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_8_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 90 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 8 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_8_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 91 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 9 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_9_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 92 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 9 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_9_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 93 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 9 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_9_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 94 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 9 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_9_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 95 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 9 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_9_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 96 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 9 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_9_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 97 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 9 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_9_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 98 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 9 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_9_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 99 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 9 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_9_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 100 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 9 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_9_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 101 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 10 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_10_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 102 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 10 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_10_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 103 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 10 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_10_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 104 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 10 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_10_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 105 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 10 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_10_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 106 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 10 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_10_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 107 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 10 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_10_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 108 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 10 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_10_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 109 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 10 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_10_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 110 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 10 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_10_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 111 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 11 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_11_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 112 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 11 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_11_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 113 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 11 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_11_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 114 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 11 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_11_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 115 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 11 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_11_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 116 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 11 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_11_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 117 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 11 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_11_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 118 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 11 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_11_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 119 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 11 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_11_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 120 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 11 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_11_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 121 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 12 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_12_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 122 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 12 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_12_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 123 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 12 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_12_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 124 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 12 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_12_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 125 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 12 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_12_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 126 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 12 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_12_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 127 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 12 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_12_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 128 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 12 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_12_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 129 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 12 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_12_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 130 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 12 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_12_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 131 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 13 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_13_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 132 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 13 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_13_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 133 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 13 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_13_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 134 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 13 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_13_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 135 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 13 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_13_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 136 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 13 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_13_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 137 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 13 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_13_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 138 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 13 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_13_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 139 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 13 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_13_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 140 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 13 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_13_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 141 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 14 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_14_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 142 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 14 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_14_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 143 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 14 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_14_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 144 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 14 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_14_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 145 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 14 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_14_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 146 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 14 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_14_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 147 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 14 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_14_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 148 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 14 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_14_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 149 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 14 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_14_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 150 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 14 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_14_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 151 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 15 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_15_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 152 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 15 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_15_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 153 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 15 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_15_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 154 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 15 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_15_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 155 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 15 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_15_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 156 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 15 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_15_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 157 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 15 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_15_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 158 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 15 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_15_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 159 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 15 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_15_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 160 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 15 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_15_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 161 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 16 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_16_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 162 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 16 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_16_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 163 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 16 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_16_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 164 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 16 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_16_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 165 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 16 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_16_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 166 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 16 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_16_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 167 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 16 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_16_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 168 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 16 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_16_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 169 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 16 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_16_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 170 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 16 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_16_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 171 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 172 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 173 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 174 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 175 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 176 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 177 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 178 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 179 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 180 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 17 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_17_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 181 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 18 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_18_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 182 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 18 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_18_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 183 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 18 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_18_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 184 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 18 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_18_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 185 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 18 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_18_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 186 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 18 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_18_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 187 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 18 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_18_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 188 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 18 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_18_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 189 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 18 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_18_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 190 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 18 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_18_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 191 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 19 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_19_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 192 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 19 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_19_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 193 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 19 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_19_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 194 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 19 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_19_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 195 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 19 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_19_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 196 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 19 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_19_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 197 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 19 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_19_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 198 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 19 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_19_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 199 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 19 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_19_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 200 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 19 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_19_bt_10
fi

if [ $SLURM_ARRAY_TASK_ID -eq 201 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 20 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_20_bt_1
fi

if [ $SLURM_ARRAY_TASK_ID -eq 202 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 20 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_20_bt_2
fi

if [ $SLURM_ARRAY_TASK_ID -eq 203 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 20 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_20_bt_3
fi

if [ $SLURM_ARRAY_TASK_ID -eq 204 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 20 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_20_bt_4
fi

if [ $SLURM_ARRAY_TASK_ID -eq 205 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 20 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_20_bt_5
fi

if [ $SLURM_ARRAY_TASK_ID -eq 206 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 20 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_20_bt_6
fi

if [ $SLURM_ARRAY_TASK_ID -eq 207 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 20 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_20_bt_7
fi

if [ $SLURM_ARRAY_TASK_ID -eq 208 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 20 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_20_bt_8
fi

if [ $SLURM_ARRAY_TASK_ID -eq 209 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 20 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_20_bt_9
fi

if [ $SLURM_ARRAY_TASK_ID -eq 210 ]
then
treemix -i allele_table/treemix_input.table.gz -root OG -k 500 -se -seed $RANDOM -m 20 -bootstrap -o rerun_randomseed_treemix_results_m20_bt10/migration_20_bt_10
fi

