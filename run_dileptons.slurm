#!/bin/bash
#SBATCH --error=ok_%A_%a.log
#SBATCH --output=qxpt_out_2_%A_%a.log
#SBATCH --job-name=Skyrme-soft
#SBATCH --partition=general1
#SBATCH --nodes=1
#SBATCH --ntasks=20	
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1024
#SBATCH --time=40:00:00
#SBATCH --mail-type=FAIL



#####################
### Configuration ###
#####################
# Blocks must be 32 CPU
# Here : 10 blocks of 32 processors = 320 jobs
# Blocks are allocated one by one (groups of 32 CPU)
# Single CPU jobs
# 700 Mb of memory/CPU
#########################
### Global parameters ###
#########################

MAINDIR=$PWD
SCRATCHDIR='/scratch/deeplearn/savchuk/potential/v2-dilepton/123b66'
PROJECTNAME='Dileptons'



INPUTFILE='inputfile'
EVENTNUM=1250
CORES=40



    FLDR=${SCRATCHDIR}/$PROJECTNAME
    echo $FLDR
    ./dilepton_v2 $FLDR ${2} ${3} 20


wait
