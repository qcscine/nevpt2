#!/bin/bash

export LSB_DJOB_NUMPROC=1
nproc=$LSB_DJOB_NUMPROC
mycurrdir=$PWD
templatedir=$PWD/../../template
chkp=*.checkpoint_state.*.*.h5
chkptemplate=$templatedir/$chkp

resfile=*.result*.*.h5

# Load necessary modules here:
#module load gcc/4.9.2 szip hdf5 boost python/2.7.6 intel

cp -rv $chkptemplate $mycurrdir
OMP_NUM_THREADS=$nproc $MOLCAS/qcmaquis/bin/dmrg_meas $mycurrdir/meas-4rdm.in > $mycurrdir/meas-4rdm.log
rm -r $chkp
