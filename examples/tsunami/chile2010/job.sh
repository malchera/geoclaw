#!/bin/bash
#SBATCH -J CHILE2010        # job name
#SBATCH -o CHILE2010_OUT_%j # output file name (%j expands to jobID)
#SBATCH -e CHILE2010_ERR_%j # error file name (%j expands to jobID)
#SBATCH -N 1                # total number of nodes requested (16 cores/node)
#SBATCH -n 1                # 1 task
#SBATCH -p serial           # queue (partition) -- normal, development, etc.
#SBATCH -t 00:10:00         # run time (hh:mm:ss)
#SBATCH --mail-user=am4491@columbia.edu
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes

make topo
make data
make output
