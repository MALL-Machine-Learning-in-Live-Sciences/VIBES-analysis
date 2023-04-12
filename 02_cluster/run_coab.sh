#SBATCH -N 1
#SBATCH -c 48
#SBATCH -t 02:00:00
#SBATCH --mem=128G
#SBATCH --error="/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/BV_Microbiome/out_data/coab_error.txt"
#SBATCH --output="/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/BV_Microbiome/out_data/coab_salida.txt"

TMP=/mnt/netapp2/Store_uni/home/ulc/co/dfe/TEMPR Rscript /mnt/netapp2/Store_uni/home/ulc/co/dfe/git/BV_Microbiome/02_coab_spieceasi.r