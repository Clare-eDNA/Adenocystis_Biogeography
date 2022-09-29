#!/bin/bash -e
#SBATCH -A uoo03666
#SBATCH -J CutAdapt
#SBATCH --time 1:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH --partition=large
#SBATCH --mail-user=c.adams@otago.ac.nz
#SBATCH --output cutAdapt.%j.out # CHANGE each run
#SBATCH --error cutAdapt.%j.err # CHANGE each run

module load cutadapt/3.5-gimkl-2020a-Python-3.8.2

cd /nesi/nobackup/uoo03666/Adenocystis/trimmed/trimmed/

echo "loaded CutAdapt, now starting trim of Fs to 70bp"

for sample in *.1.trim.fq.gz
do
base=$(basename ${sample} .1.trim.fq.gz)
echo "${base}"
cutadapt -l 70 -o "${base}".70bp.1.fq.gz "${base}".1.trim.fq.gz -m 70 -j 8
done

echo "done with forwards trimming, moving onto reverse"
# l trims down the length to 70 bp
# m discards any reads below 70 bp
# j is the number of threads/cores
# all reads should now be 70 bp

for sample in *.2.trim.fq.gz
do
base=$(basename ${sample} .2.trim.fq.gz)
echo "${base}"
cutadapt -l 70 -o "${base}".70bp.2.fq.gz "${base}".2.trim.fq.gz -m 70 -j 8
done

echo "done with trimming to 70bp for all"

mkdir 70bp
mv *70bp*gz 70bp/

