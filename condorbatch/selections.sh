#!/bin/bash

#This is the first LCG environment with uproot3 and boost_histogram together
LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt 

#untar your crap
echo "Untaring  directory with selections code"
tar -xf selectionsForCondor.tar.gz
cd selection_jobs

#source the environment
echo "Sourcing the environment"
source $LCG/setup.sh

#debug
echo $PYTHONPATH

#Arguments taken
#1 - the eos directory to drop the output
#2 - the topiaries
#3 - the channel 
#4 - the special scenario
#5 - the event weight systematic (command line options to selector)

#Running the selection maker
echo "Beginning analysis"
python runSelections.py -s $2 -c $3 -conf $4 -evntwsyst $5

for FILE in analysis_output_ZpAnomalon/*/*
do
    echo ${FILE}
    echo "copying ${FILE} to eos $1"
    xrdcp ${FILE} $1/${FILE}
    XRDEXIT=$?
    if [[ $XRDEXIT -ne 0 ]]; then
	rm ${FILE}
	echo "failure in xrdcp, exit code $XRDEXIT"
	exit $XRDEXIT
    fi
    rm ${FILE}
done


#cleanup
cd ${_CONDOR_SCRATCH_DIR_}
rm -rf selection_jobs
rm selectionsForCondor.tar.gz

