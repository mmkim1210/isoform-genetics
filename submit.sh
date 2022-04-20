#!/bin/bash
#$ -cwd
#$ -e ./log/joberr.$JOB_ID.$TASK_ID
#$ -o ./log/joblog.$JOB_ID.$TASK_ID 
#$ -j y 
#$ -l h_rt=24:00:00,h_data=8G,arch=intel*,highp
#$ -pe shared 4
#$ -M $USER@mail 
#$ -m a
#$ -t 1-25000:100
#$ -tc 20

. /u/local/Modules/default/init/modules.sh
module load julia/1.7

echo ${SGE_TASK_ID}
julia --project=. ./src/vc.jl ${SGE_TASK_ID} > ./log/output.$JOB_ID.${SGE_TASK_ID} 2>&1
