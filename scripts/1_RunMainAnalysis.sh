#!/bin/bash
#SBATCH --qos castles
#SBATCH --partition=icelake-castles
#SBATCH --mail-type ALL
#SBATCH --nodes 1
#SBATCH --ntasks 70
#SBATCH --time 3:0:0
#SBATCH --mem 400G
#SBATCH --job-name Blackquick



module purge;
module load bluebear
module load bear-apps/2021b
module load R/4.2.0-foss-2021b

Rscript /castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/1_MainAnalysis.R
#before July2
