#!/bin/bash

#PBS -N 88soils_ph
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=10gb
#PBS -V
#PBS -j oe
#PBS -d .

set -e
cpus=$PBS_NUM_PPN

export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER/Stability_2020
[ ! -d $TMPDIR ] && mkdir $TMPDIR
export TMPDIR=$TMPDIR/data_applications
[ ! -d $TMPDIR ] && mkdir $TMPDIR
#tmp=$(mktemp -d --tmpdir)
#export TMPDIR=$tmp
#trap "rm -r $tmp; unset TMPDIR" EXIT

# do something
source activate r-c-env
Rscript 88soils_stab_application.R $TMPDIR
source deactivate r-c-env

#mv $tmp/outdir ./outdir
