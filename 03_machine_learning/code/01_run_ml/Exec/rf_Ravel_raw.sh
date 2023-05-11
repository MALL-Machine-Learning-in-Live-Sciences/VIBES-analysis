#!/bin/bash 
#SBATCH -p short 
#SBATCH -t 06:00:00 
#SBATCH --mem=120GB 
#SBATCH -N 1 
#SBATCH -n 24 
name=Ravel_raw
data=03_machine_learning/data/Ravel_raw.rds
outdir=03_machine_learning/res/aprox-2/
Rscript 03_machine_learning/code/ml//modelsrf.r $data $name $outdir