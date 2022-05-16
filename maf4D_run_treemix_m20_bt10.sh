#!/bin/bash
#SBATCH -p medium
#SBATCH -a 1-21
#SBATCH -C scratch
#SBATCH -n 4
#SBATCH --qos long
#SBATCH -t 120:00:00

module load anaconda3
source activate
conda activate test


if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 0 -bootstrap -o treemix_results_m20_bt10/migration_0_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 1 -bootstrap -o treemix_results_m20_bt10/migration_1_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 2 -bootstrap -o treemix_results_m20_bt10/migration_2_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 4 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 3 -bootstrap -o treemix_results_m20_bt10/migration_3_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 5 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 4 -bootstrap -o treemix_results_m20_bt10/migration_4_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 6 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 5 -bootstrap -o treemix_results_m20_bt10/migration_5_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 7 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 6 -bootstrap -o treemix_results_m20_bt10/migration_6_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 8 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 7 -bootstrap -o treemix_results_m20_bt10/migration_7_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 9 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 8 -bootstrap -o treemix_results_m20_bt10/migration_8_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 10 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 9 -bootstrap -o treemix_results_m20_bt10/migration_9_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 11 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 10 -bootstrap -o treemix_results_m20_bt10/migration_10_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 12 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 11 -bootstrap -o treemix_results_m20_bt10/migration_11_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 13 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 12 -bootstrap -o treemix_results_m20_bt10/migration_12_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 14 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 13 -bootstrap -o treemix_results_m20_bt10/migration_13_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 15 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 14 -bootstrap -o treemix_results_m20_bt10/migration_14_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 16 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 15 -bootstrap -o treemix_results_m20_bt10/migration_15_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 17 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 16 -bootstrap -o treemix_results_m20_bt10/migration_16_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 18 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 17 -bootstrap -o treemix_results_m20_bt10/migration_17_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 19 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 18 -bootstrap -o treemix_results_m20_bt10/migration_18_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 20 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 19 -bootstrap -o treemix_results_m20_bt10/migration_19_bt_${i}
done

fi


if [ $SLURM_ARRAY_TASK_ID -eq 21 ]
then

for i in {1..10}
do
	treemix -i allele_table/treemix_input.table.gz -root OG -k 500   -se -m 20 -bootstrap -o treemix_results_m20_bt10/migration_20_bt_${i}
done

fi