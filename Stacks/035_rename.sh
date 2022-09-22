#!/bin/bash -e
#SBATCH -A uoo03666
#SBATCH -J ParamTest
#SBATCH --time 08:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --mem=64G
#SBATCH --partition=bigmem,infill
#SBATCH --mail-user=c.adams@otago.ac.nz
#SBATCH --output optomization.%j.out # CHANGE each run
#SBATCH --error optomization.%j.err # CHANGE each run


#module load Stacks/2.61-gimkl-2022a


cd /nesi/nobackup/uoo03666/Adenocystis/trimmed/trimmed/70bp/LargeFiles/

for sample in *.70bp.1.fq.gz
do
base=$(basename ${sample} .70bp.1.fq.gz)
echo "${base}"
cp "${sample}" "${base}".1.fq.gz
done

for sample in *.70bp.2.fq.gz
do
base=$(basename ${sample} .70bp.2.fq.gz)
echo "${base}"
cp "${sample}" "${base}".2.fq.gz
done
