#!/bin/bash
#SBATCH -N 1
#SBATCH -p thinnodes
#SBATCH -t 00:30:00
#SBATCH --mem=32G
#SBATCH --error="/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/NC_error.txt"
#SBATCH --output="/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/NC_salida_R.txt"

module load gcc/6.4.0 R/4.0.2
TMP=/mnt/netapp2/Store_uni/home/ulc/co/dfe/TEMPR Rscript /home/ulc/co/dfe/git/Entropy/BatchEffect_Cluster/NO_CLUST_MMPUHI.R