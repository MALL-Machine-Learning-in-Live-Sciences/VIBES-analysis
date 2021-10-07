#!/bin/bash
#SBATCH -N 1
#SBATCH -p thinnodes
#SBATCH -t 03:00:00
#SBATCH --mem=32G
#SBATCH --error="/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/ML3C_error.txt"
#SBATCH --output="/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/ML3C_salida.txt"

module load gcc/6.4.0 R/4.0.2
TMP=/mnt/netapp2/Store_uni/home/ulc/co/dfe/TEMPR Rscript /mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy2/ML3C.R