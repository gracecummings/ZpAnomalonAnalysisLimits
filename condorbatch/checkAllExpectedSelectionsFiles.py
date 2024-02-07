import glob
import json
import sys
import os
import argparse
from datetime import date
import subprocess

parser = argparse.ArgumentParser()

if __name__=='__main__':
    parser.add_argument("-j","--sampleJson",help=".json file with the skims to be made into topiary")
    parser.add_argument("-d","--date",type=str,help="when you submitted to see if all worked")
    parser.add_argument("-c","--channel",type=str,help="channel: mumu,ee,emu")
    parser.add_argument("-r","--region_config",type=str,help="special region to be considerd, ie 'unblind', 'alpham-basic','alpham-validation'")
    #parser.add_argument("-syst","--syststr",type=str,help="the systematic flag you used, again, to check")
    args = parser.parse_args()

    #Check json
    expsamps = []
    if (args.sampleJson is None):
        print("No list of skims provided, please provide appropriate json")
        fsjson = {}
    else:
        if (not os.path.exists(args.sampleJson)):
            print("Invalid json")
        fstotop = open(args.sampleJson,'r')
        fsjson = json.load(fstotop)

    if len(fsjson.keys()) > 0:
        for key in fsjson.keys():
            sampleName = key.split("_topiary")[0]
            expsamps.append(sampleName)
    
   #Get the files that were made
    eospath = "/store/user/lpcboostres/leptophobic_selections-_"+args.date+"/analysis_output_ZpAnomalon/"+args.date+"/"
    fs = subprocess.check_output("eos root://cmseos.fnal.gov ls "+eospath,shell=True).decode(sys.stdout.encoding).split()


    #What were you expecting to make?
    #right now for data only
    ftypes = [["upout",".root"],["totalevents",".npy"],["selected_errors",".pkl"]]
    regs   = ["signalr","sideband","sideband_alphat"]
    cuts   = "systnominal_btagnom_muidnom_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8"
    expregfs = []
    for ft in ftypes:
        for reg in regs:
            expregfs.append(ft[0]+"_"+reg+"_"+cuts+ft[1])

    #do the checks
    for samp in expsamps:
        expfs = [samp+"_"+freg for freg in expregfs]
        obsfs = [f for f in fs if samp in f]
        missing = list(set(obsfs) - set(expfs))
        print("Checking sample: ",samp)
        print("    number missing: ",len(missing))
