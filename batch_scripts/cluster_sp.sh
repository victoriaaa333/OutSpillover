#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=baobao_job
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01-00:00:00
#SBATCH --mail-type=ALL

# do something

ml load "R/4.1.0-foss-2020b"

#cd /home/yb254/project/OutSpillover

cd /gpfs/ysm/project/forastiere/yb254/OutSpillover

Rscript simulation_scripts/cluster_sp_2.R
