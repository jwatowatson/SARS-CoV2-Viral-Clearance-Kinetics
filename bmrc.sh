#!/bin/bash
#$ -P moru-batty.prjc
#$ -wd /well/moru-batty/users/gka292/SARS-CoV2-Viral-Clearance-Kinetics
#$ -N res
#$ -pe shmem 20
#$ -o /well/moru-batty/users/gka292/SARS-CoV2-Viral-Clearance-Kinetics/o_and_e_files
#$ -e /well/moru-batty/users/gka292/SARS-CoV2-Viral-Clearance-Kinetics/o_and_e_files
#$ -q short.qf
#$ -t 1-1000
#$ -tc 32

echo started=`date`
module purge
module load R/4.1.2-foss-2021b

echo "job=$JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`
Rscript /well/moru-batty/users/gka292/SARS-CoV2-Viral-Clearance-Kinetics/Sample_sizes_sims_stan.R ${SGE_TASK_ID} --no-save --no-restore
echo "finished="`date`
exit 0
