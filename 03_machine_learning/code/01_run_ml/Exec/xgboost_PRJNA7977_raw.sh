#!/bin/bash 
#SBATCH -p short 
#SBATCH -t 06:00:00 
#SBATCH --mem=120GB 
#SBATCH -N 1 
#SBATCH -n 24 
name=PRJNA7977_raw
data=03_machine_learning/data/PRJNA7977_raw.rds
outdir=03_machine_learning/res/aprox-2/
Rscript 03_machine_learning/code/ml//modelsxgboost.r $data $name $outdir