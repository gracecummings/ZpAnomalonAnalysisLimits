#!/bin/bash

#untar your crap
echo "Untaring  directory with topiary code"
tar -xf topiaryForCondor.tar.gz
cd topiary_jobs


#source the environment
echo "Sourcing the environment"
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
source setup_RestFrames.sh
echo ${RESTFRAMESSYS}

#Arguments taken
#1 - the eos directory to drop the topiary
#2 - the sample name
#3 - the files to put into the topiary
#4 - the channel 
#5 - the systematic selection

#Running topiary 
echo "Making topiary"
python runTopiary.py -s $2 -c $4 -l $3 -syst $5

for FILE in ./*.root
do
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
rm -rf topiary_jobs
rm topiaryForCondor.tar.gz
