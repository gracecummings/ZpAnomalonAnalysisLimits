#!/usr/bin/env bash

for f in $1/fitDiagnosticsZp4000ND800NS200_expectedsignal*;
do
    echo $f
    python makeSignalInjectionTestPlots.py -f $f
done
   
