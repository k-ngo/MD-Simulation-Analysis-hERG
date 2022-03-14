#!/bin/bash

#SBATCH --job-name=SimAnalysis
#SBATCH -p bigmem
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem=140G
#SBATCH -o slurm-%A_%a.out
#SBATCH -e slurm-%A_%a.err
#SBATCH --mail-user=khoango@ucdavis.edu
#SBATCH --mail-type=ALL

module load vmd
python simulation_analysis.py -d *.trr -t 8090 -s 3 --drug HETB