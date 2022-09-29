#!/bin/bash -e
#SBATCH -A uoo03666
#SBATCH -J CutAdapt
#SBATCH --time 2:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH --partition=large
#SBATCH --mail-user=c.adams@otago.ac.nz
#SBATCH --output cutAdapt.%j.out # CHANGE each run
#SBATCH --error cutAdapt.%j.err # CHANGE each run

module load cutadapt/3.5-gimkl-2020a-Python-3.8.2

cd /nesi/nobackup/uoo03666/Adenocystis/trimmed/

for sample in *.1.fq.gz
do
base=$(basename ${sample} .1.fq.gz)
echo "${base}"
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 8 -o /nesi/nobackup/uoo03666/Adenocystis/trimmed/"${base}".1.trim.fq.gz -p //nesi/nobackup/uoo03666/Adenocystis/trimmed/"${base}".2.trim.fq.gz "${base}".1.fq.gz "${base}".2.fq.gz -m 25:20 -q 10
done


#mkdir trimmed
#mv *trim*gz trimmed/
