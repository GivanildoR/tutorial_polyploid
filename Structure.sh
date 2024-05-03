#!/bin/bash
#SBATCH --job-name=Structure/array.%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
#SBATCH --mem=150GB
#SBATCH --qos=mresende-b
#SBATCH --account=mresende
#SBATCH -t 96:00:00
#SBATCH --output=1.tmp/imp.%a.array.%A.out
#SBATCH --error=2.error/imp.%a.array.%A.err
#SBATCH --mail-user=givanido.rodrigu@ufl.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-10

module load ufrc
module load structure


# input
K=$SLURM_ARRAY_TASK_ID


# Define input file containing allele count data
input_file="/blue/mresende/share/Givanildo/Structure_files/RUN_STRUCTURE/structure_output.txt"

# Define main parameter file
main_params="/blue/mresende/share/Givanildo/Structure_files/RUN_STRUCTURE/mainparams.txt"

# Define extra parameter file
extra_params="/blue/mresende/share/Givanildo/Structure_files/RUN_STRUCTURE/extraparams.txt"


# Define a function to run STRUCTURE

for rep in 1 2 3
do

  structure -i "$input_file" -o "output_K${K}_replicate${rep}" -m "$main_params" -K "$K" -e "$extra_params"

done
