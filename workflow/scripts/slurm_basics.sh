#!/bin/bash

#SBATCH --job-name=hg5SV       # Job name
#SBATCH --output=/path/to/output/sv_hg5_t24def.%J.out         # Output file name
#SBATCH --error=/path/to/error/sv_hg5_t24def.%J.err          # Error file name
#SBATCH --partition=all,highmem               # Partition short max 1h, use all for more
#SBATCH --time=24:00:00                 # Time limit
# SBATCH --partition=gpu
# SBATCH --gres=gpu:1
# SBATCH --ntasks=1
# SBATCH --cpus-per-task=16

#SBATCH --mem=190GB						# RAM memory
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks-per-node=24             # MPI processes per node, only 1 node is reserved !
# SBATCH -n 72 				# reserves n tasks/cores on multiple shared nodes (24 x 4 nodes)


echo "Running on hosts: $SLURM_JOB_NODELIST"
echo "Running on $SLURM_JOB_NUM_NODES nodes."
echo "Running $SLURM_NTASKS tasks."
echo "Account: $SLURM_JOB_ACCOUNT"
echo "Job ID: $SLURM_JOB_ID"
echo "Job name: $SLURM_JOB_NAME"
echo "Node running script: $SLURMD_NODENAME"

hostname
echo 'workdir:'
echo 'optional Parameters:' #-s=5 -maxsv=1   maxsv=2 results in empty outputs file error...

