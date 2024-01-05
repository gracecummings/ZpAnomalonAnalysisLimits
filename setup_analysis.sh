#!/usr/bin/env bash

NAME=physanal
LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc11-opt

if [[ -f $NAME/bin/activate ]]; then
  echo "$NAME already installed. Run \`source $NAME/bin/activate\` to activate"
else
    source $LCG/setup.sh
    python -m venv --copies $NAME
    source $NAME/bin/activate
    LOCALPATH=$(python -c 'import sys; print(f"{sys.prefix}/lib/python{sys.version_info.major}.{sys.version_info.minor}/site-packages")')
    export PYTHONPATH=${LOCALPATH}:$PYTHONPATH
    python -m pip install fsspec-xrootd
    sed -i "2a source ${LCG}/setup.sh" $NAME/bin/activate
    sed -i "3a export PYTHONPATH=${LOCALPATH}:\$PYTHONPATH" $NAME/bin/activate
    ipython kernel install --user --name=$NAME
fi
