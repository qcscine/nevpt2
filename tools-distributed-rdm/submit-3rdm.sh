#!/bin/bash

export LSB_DJOB_NUMPROC=1
nproc=$LSB_DJOB_NUMPROC
mycurrdir=$PWD
templatedir=$PWD/../../template
chkp=*.checkpoint_state.*.*.h5
chkptemplate=$templatedir/$chkp

resfile=*.result*.*.h5

rsync -aL $templatedir/$fcid $chkptemplate $mycurrdir
OMP_NUM_THREADS=$nproc $MOLCAS/qcmaquis/bin/dmrg_meas $mycurrdir/meas-3rdm.in > $mycurrdir/meas-3rdm.log

rm -r $chkp
