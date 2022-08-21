#!/bin/bash
#SBATCH --job-name=test1              # create a short name for your job
#SBATCH --nodes=1                           # node count
#SBATCH --ntasks=1                          # total number of tasks across all nodes
#SBATCH --cpus-per-task=4                   # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=5G                    # memory per cpu-core
#SBATCH --time=23:59:99                     # total run time limit (HH:MM:SS)
#SBATCH --output="test.out"
#SBATCH --error="test.err"
#SBATCH --mail-type=fail                   # notifications for job done & fail
#SBATCH --mail-user=fangweic@princeton.edu  # send-to address


module add julia/1.5.0
module add gurobi
julia Run.jl

date
