#!/bin/bash

#PBS -N ind_cv
#PBS -l walltime=500:00:00
#PBS -l nodes=1:ppn=10
#PBS -l mem=50gb
#PBS -V
#PBS -j oe
#PBS -d .

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
Rscript ind_results.R $TMPDIR
Rscript block_results.R $TMPDIR
Rscript toe_results.R $TMPDIR
Rscript boot_RF_testing.R $TMPDIR
source deactivate r-env

#mv $tmp/outdir ./outdir
