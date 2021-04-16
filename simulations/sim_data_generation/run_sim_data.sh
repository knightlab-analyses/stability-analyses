#!/bin/bash

#PBS -N sim_ind_toe
#PBS -l walltime=500:00:00
#PBS -l nodes=1:ppn=4
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
Rscript sim_dat_ind_toeplitz.R $TMPDIR
Rscript sim_dat_block.R $TMPDIR
source deactivate r-env

#mv $tmp/outdir ./outdir
