#!/bin/bash
#SBATCH -J scRNA-seq-analysis
#SBATCH -p nonfs
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=2
#SBATCH --mem=2000
#SBATCH -t 0-5:00
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH --qos=low
