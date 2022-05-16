#!/bin/bash
#SBATCH -p medium
#SBATCH -C scratch
#SBATCH -a 1-30
#SBATCH -n 4
#SBATCH -t 24:00:00

module load gatk

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
gatk VariantsToTable -V ./di_BS.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_BS_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 2 ]
then
gatk VariantsToTable -V ./di_BSJ.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_BSJ_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 3 ]
then
gatk VariantsToTable -V ./di_DS.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_DS_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 4 ]
then
gatk VariantsToTable -V ./di_GQW.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_GQW_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 5 ]
then
gatk VariantsToTable -V ./di_HYXG.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_HYXG_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 6 ]
then
gatk VariantsToTable -V ./di_JZG.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_JZG_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 7 ]
then
gatk VariantsToTable -V ./di_LJSa.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_LJSa_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 8 ]
then
gatk VariantsToTable -V ./di_PQG.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_PQG_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 9 ]
then
gatk VariantsToTable -V ./di_QA.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_QA_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 10 ]
then
gatk VariantsToTable -V ./di_SWP.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_SWP_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 11 ]
then
gatk VariantsToTable -V ./di_TBS.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_TBS_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 12 ]
then
gatk VariantsToTable -V ./di_TS.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_TS_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 13 ]
then
gatk VariantsToTable -V ./di_WL.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_WL_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 14 ]
then
gatk VariantsToTable -V ./di_WTS.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_WTS_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 15 ]
then
gatk VariantsToTable -V ./di_XLS.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_XLS_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 16 ]
then
gatk VariantsToTable -V ./di_XM.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_XM_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 17 ]
then
gatk VariantsToTable -V ./di_XWT.args.vcf.gz -F CHROM -F POS -F AC -F AN -O di_XWT_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 18 ]
then
gatk VariantsToTable -V ./tetra_CSBT.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_CSBT_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 19 ]
then
gatk VariantsToTable -V ./tetra_CZPT.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_CZPT_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 20 ]
then
gatk VariantsToTable -V ./tetra_GST.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_GST_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 21 ]
then
gatk VariantsToTable -V ./tetra_GTST.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_GTST_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 22 ]
then
gatk VariantsToTable -V ./tetra_LJSb.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_LJSb_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 23 ]
then
gatk VariantsToTable -V ./tetra_QA.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_QA_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 24 ]
then
gatk VariantsToTable -V ./tetra_SNJ.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_SNJ_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 25 ]
then
gatk VariantsToTable -V ./tetra_TB.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_TB_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 26 ]
then
gatk VariantsToTable -V ./tetra_TBS.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_TBS_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 27 ]
then
gatk VariantsToTable -V ./tetra_WL.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_WL_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 28 ]
then
gatk VariantsToTable -V ./tetra_XJS.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_XJS_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 29 ]
then
gatk VariantsToTable -V ./tetra_XLS.args.vcf.gz -F CHROM -F POS -F AC -F AN -O tetra_XLS_tm.table
fi

if [ $SLURM_ARRAY_TASK_ID -eq 30 ]
then
gatk VariantsToTable -V ./OG.args.vcf.gz -F CHROM -F POS -F AC -F AN -O OG_tm.table
fi