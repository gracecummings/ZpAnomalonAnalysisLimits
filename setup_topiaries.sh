#!/usr/bin/env bash

NAME=physpy
LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt

source $LCG/setup.sh

if [[ -f $NAME/bin/activate ]]; then
  echo "$NAME already installed. Run \`source $NAME/bin/activate\` to activate"
else
    python -m venv --copies $NAME
    source $NAME/bin/activate
    python -m pip install boost-histogram
fi


