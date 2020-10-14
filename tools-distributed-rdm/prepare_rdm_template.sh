#!/bin/bash
#set -x
t3rdmmode=0
spin=0
states=()
submitfile=$MOLCAS/Tools/distributed-4rdm/submit-4rdm.sh
jobmanager=$MOLCAS/Tools/distributed-4rdm/jobmanager.py

command_usage()
{
  echo "Usage: "
	echo " $0 [ -3 ] [ -s <SPIN> ] <state range> "
	echo " $0 prepares the scratch directory for the distributed 4-RDM or t-3-RDM evaluation."
	echo " -3 -- prepare t-3RDM evaluation instead of 4-RDM"
	echo " -s <SPIN> -- use the checkpoint of the wavefunction with the total spin S=SPIN, if multiple are available (0 -- singlet, 1 -- doublet, 2 -- triplet etc.) Default -- singlet."
	echo " <state range> (optional) list of states for which the directories should be prepared. Default: all states"
}

if [ $# -ne 0 ] && [ ${1} == "-h" ] && [ ${1} == "--help" ] ; then
  command_usage
  exit 1
fi

if [ $# -gt 0 ]; then
  while [ x${1} != "x" ]; do
    if [ ${1} == "-3" ]; then
      t3rdmmode=1
    elif [ ${1} == "-s" ]; then
      shift
      spin=${1}
    else
      states+=(${1})
    fi
    shift
  done
fi

# if no states are specified, populate the states array with all states. the states present are extracted from the meas-4rdm.X.in files present in the directory
if [ ${#states[@]} -eq 0 ]; then
  states=($(eval 'for i in meas-4rdm.*; do b=${i%%.in}; echo ${b##meas-4rdm.}; done | sort -n'))
fi

# Extract the MOLCAS project name from the name of OneInt file. we will need it later
oneint=$(ls *.OneInt)
if [ x${oneint} != "x" ]; then
  molcasproject=${oneint%%.OneInt}
else
  echo "Cannot extract MOLCAS project. Are you sure you are in the MOLCAS scratch directory?"
  exit 1
fi

## for each state, populate a scratch directory
if [ $t3rdmmode -eq 0 ]; then
  echo -n "Creating scratch files for distributed 4-RDM evaluation for state "
  for state in ${states[@]}; do
    echo -n "${state} "
    dir="4rdm-scratch."${state}

    # obtain checkpoint name from the resultfile parameter in the QCMaquis input file template
    checkpointstate=$(grep 'chkpfile' meas-4rdm.${state}.in | awk '{print $3}')
    if [ ! -e $checkpointstate ]; then
	echo "Checkpoint $checkpointstate not found in scratch directory."
	exit 1
    fi

    # if it exists, copy the necessary files into the scratch directory
    mkdir -p $dir/template
    mkdir -p $dir/parts
    cp -r $checkpointstate $dir/template
    cp meas-4rdm.${state}.in $dir/template/meas-4rdm-template.in
    cp $submitfile $dir/template/
    cd $dir/template/ && sed -i "s/state=0/state=${state}/g" submit-4rdm.sh && cd ../..
    cp $jobmanager $dir/
  done
  echo "Done."
else
  echo -n "Creating scratch files for distributed t-3-RDM evaluation for states "
  submitfile=$MOLCAS/Tools/distributed-4rdm/submit-3rdm.sh
  for statei in ${states[@]}; do
    for statej in ${states[@]}; do
      if [ ${statei} -lt ${statej} ]; then
        echo -n "${statei}-${statej} "
        dir="3rdm-scratch."${statei}"."${statej}
        chkpi=$(grep 'chkpfile' meas-3tdm.${statei}.${statej}.in | awk '{print $3}')
        chkpj=$(grep 'MEASURE\[trans3rdm' meas-3tdm.${statei}.${statej}.in | awk '{print $3}' | awk -F\; '{gsub("\"",""); print $1}')
	if [ ! -e $chkpi ]; then
	   echo "Checkpoint $chkpi not found in scratch directory."
	   exit 1
	fi
	if [ ! -e $chkpj ]; then
	   echo "Checkpoint $chkpj not found in scratch directory."
	   exit 1
	fi
        mkdir -p $dir/template
        mkdir -p $dir/parts
        cp -r $chkpi $chkpj $dir/template
        cp meas-3tdm.${statei}.${statej}.in $dir/template/meas-3rdm-template.in
        cp $submitfile $dir/template/
        cd $dir/template/ && sed -i "s/statei=0/statei=${statei}/g" submit-3rdm.sh && sed -i "s/statej=0/statej=${statej}/g" submit-3rdm.sh && cd ../..
        cp $jobmanager $dir/
      fi
    done
  done
  echo "Done."
fi
