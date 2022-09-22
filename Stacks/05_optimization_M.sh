#!/bin/bash -e
#SBATCH -A uoo03666
#SBATCH -J ParamTest
#SBATCH --time 10:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --mem=64G
#SBATCH --partition=bigmem,infill
#SBATCH --mail-user=c.adams@otago.ac.nz
#SBATCH --output optomization.%j.out # CHANGE each run
#SBATCH --error optomization.%j.err # CHANGE each run


module load Stacks/2.61-gimkl-2022a

cd /nesi/nobackup/uoo03666/Adenocystis/trimmed/trimmed/70bp/LargeFiles/

#from source folder
for i in 2 3 4 5 6 7 8
do
mkdir -p M$i
echo '#!/bin/sh' > runM$i.sh
echo "module load Stacks/2.61-gimkl-2022a" >> runM$i.sh
echo "cd /nesi/nobackup/uoo03666/Adenocystis/trimmed/trimmed/70bp/LargeFiles/"
echo "denovo_map.pl --samples /nesi/nobackup/uoo03666/Adenocystis/trimmed/trimmed/70bp/LargeFiles/ --popmap [pop map]  -o M$i -p 0.8 -M $i -n $i -m 3 -T 8 -X \"populations: --write-random-snp\"" >> runM$i.sh
sbatch -A uoo03666 -t 10:00:00 -J M$i -c 8 --mem=64G runM$i.sh # specific to mahuika and user
done
