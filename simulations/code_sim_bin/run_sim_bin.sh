#!/bin/bash

#PBS -N sim_stab
#PBS -l walltime=500:00:00
#PBS -l nodes=1:ppn=10
#PBS -l mem=20gb
#PBS -V
#PBS -j oe
#PBS -d .
#PBS -t 1-4%4
#PBS -o messages_outputs/
#PBS -e messages_errors/

set -e
cpus=$PBS_NUM_PPN

export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER/Stability_2020
[ ! -d $TMPDIR ] && mkdir $TMPDIR
export TMPDIR=$TMPDIR/sim_data
[ ! -d $TMPDIR ] && mkdir $TMPDIR
#tmp=$(mktemp -d --tmpdir)
#export TMPDIR=$tmp
#trap "rm -r $tmp; unset TMPDIR" EXIT

# do something
source activate r-env
Rscript ind_results_binary_update.R $TMPDIR # t = 1-16%10
Rscript toe_results_binary_update.R $TMPDIR # t = 1-80%10
Rscript block_results_binary_update.R $TMPDIR # t = 1-80%10
Rscript boot_sim_binary.R $TMPDIR # t = 1-4%4
source deactivate r-env

#mv $tmp/outdir ./outdir
