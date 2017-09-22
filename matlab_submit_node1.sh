#!/bin/bash -x
#SBATCH --partition=ivy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40 # 48 (no HT)
#SBATCH --ntasks-per-core=2  # 1 (no HT)
#SBATCH --job-name=Copula_downscaling
#SBATCH --mail-user=christof.lorenz@kit.edu
#SBATCH --mail-type=ALL

echo "I ran on:"
cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST

ulimit -s unlimited

module load app/matlab/2015b

matlab -nodisplay -nodesktop -r "cd('/home/lorenz-c/SMOS_Downscaling'); copula_downscaling('/home/lorenz-c/SMOS_Downscaling/SMOS_val_downscaling_settings.mat', 1, 1, false)"
