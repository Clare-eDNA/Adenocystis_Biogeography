#!/bin/bash -e
#SBATCH -A uoo03666
#SBATCH -J lil-m-test
#SBATCH --time 02:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --mem=64G
#SBATCH --partition=bigmem,infill
#SBATCH --mail-user=c.adams@otago.ac.nz
#SBATCH --output lil-m-test.%j.out # CHANGE each run
#SBATCH --error lil-m-test.%j.err # CHANGE each run
module load Stacks/2.61-gimkl-2022a
cd /nesi/nobackup/uoo03666/Adenocystis/trimmed/trimmed/70bp/LargeFiles/
#from source folder
for i in 2 3 4 5
do
mkdir -p lil_m$i
echo '#!/bin/sh' > run_m$i.sh
echo "module load Stacks/2.61-gimkl-2022a" >> run_m$i.sh
echo "denovo_map.pl --samples /nesi/nobackup/uoo03666/Adenocystis/trimmed/trimmed/70bp/LargeFiles --popmap samples.txt  -o lil_m$i -p 0.8 -M 3 -n 3 -m $i -T 8 -X \"populations:--vcf --write-random-snp\"" >> run_m$i.sh
sbatch -A uoo03666 -t 08:00:00 -J m$i -c 8 --mem=64G run_m$i.sh # specific to mahuika and ludovic.dutoit
done
