#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=baobao_job
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01-00:00:00
#SBATCH --mail-type=ALL

# do something

ml load "R/4.1.0-foss-2020b"

#cd /home/yb254/project/OutSpillover

cd /gpfs/ysm/project/forastiere/yb254/OutSpillover

Rscript batch_scripts/no_parellel/n1.R