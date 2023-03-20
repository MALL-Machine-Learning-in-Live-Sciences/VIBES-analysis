#!/bin/bash
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH --mem=64G
#SBATCH --error="/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/BV_Microbiome/out_data/16S_error.txt"
#SBATCH --output="/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/BV_Microbiome/out_data/16S_salida.txt"

TMP=/mnt/netapp2/Store_uni/home/ulc/co/dfe/TEMPR Rscript /mnt/netapp2/Store_uni/home/ulc/co/dfe/git/BV_Microbiome/01_pipeline_16S.r