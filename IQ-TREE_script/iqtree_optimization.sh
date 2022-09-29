#!/bin/bash -e
#SBATCH -A uoo03666
#SBATCH -J iqTreeModSelec
#SBATCH --time 5-00:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --mem=64G
#SBATCH --mail-user=c.adams@otago.ac.nz
#SBATCH --output iqTreeModSelec.%j.out # CHANGE each run
#SBATCH --error iqTreeModSelec.%j.err # CHANGE each run


module load IQ-TREE/2.2.0.5-gimpi-2022a
cd /nesi/nobackup/uoo03666/Adenocystis/trimmed/trimmed/70bp/LargeFiles/lil_m2/SNPfiltering/otherformats/iqtree

iqtree2 -s final.recode.min4.phy -m MFP -AIC -mtree -cmax 8

#iqtree calls iqtree software
#-s is your phylip file
#-m is the model; MFP stands for ModelFinderPlus
#AIC is AIC rather than BIC for model choosing
#cmax is 16 cores
