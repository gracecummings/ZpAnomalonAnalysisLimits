import argparse
import glob
import gecorg_test as go
import numpy as np
import configparser
import subprocess
import sys

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    args = parser.parse_args()

    #Get systematic info                                                                                         
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()

    #naming conventions
    regions = ['totalr','sideband','signalr','totalr_alphat','signalr_alphat','sideband_alphat']
    ftype   = 'upout'
    cutstring = 'DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(args.zptcut)+'_Hptcut'+str(args.hptcut)+'_metcut'+str(args.metcut)+'_btagwp'+str(args.btagwp)

    #open list with samples
    f = open('list_all_background_mc_samples.txt','r')
    samples = f.readlines()
    samples = [x[:-1] for x in samples]

    #are they all there?
    for region in regions:
        print("Checking the {} samples".format(region))
        for syst in systs[1:]:
            print("    Checking the {} samples".format(syst))
            realpathup = config.get(syst,'pathup')
            realpathdn = config.get(syst,'pathdwn')
            logpathup = realpathup.split("root://cmseos.fnal.gov/")[-1]
            logpathdn = realpathdn.split("root://cmseos.fnal.gov/")[-1]
            syststrup = config.get(syst,'strup')
            syststrdn = config.get(syst,'strdwn')
            upstr = ftype+"_"+region+"_"+syststrup+"_"+cutstring
            dnstr = ftype+"_"+region+"_"+syststrdn+"_"+cutstring
            upquery = "eos root://cmseos.fnal.gov/ ls "+logpathup+" | grep "+upstr
            dnquery = "eos root://cmseos.fnal.gov/ ls "+logpathdn+" | grep "+dnstr
            #print(upquery)
            filesup  = subprocess.check_output(upquery,shell=True).decode(sys.stdout.encoding).split()
            filesdn  = subprocess.check_output(dnquery,shell=True).decode(sys.stdout.encoding).split()
            
            expup = [x+"_"+upstr+".root" for x in samples]
            expdn = [x+"_"+dnstr+".root" for x in samples]
            
            missup = sorted(list(set(expup) - set(filesup)))
            missdn = sorted(list(set(expdn) - set(filesdn)))
            
            if len(missup) > 0:
                print("      Missing the following:" )
                for m in missup:
                    print("        ",m)
            if len(missdn) > 0:
                print("      Missing the following:" )
                for m in missdn:
                    print("        ",m)
            else:
                print("Not Missing anything!")
        for syst in systs[:1]:
            print("    Checking the {} samples".format(syst))
            realpath = config.get(syst,'pathnom')
            logpath = realpath.split("root://cmseos.fnal.gov/")[-1]
            syststr = config.get(syst,'strnom')
            fullstr = ftype+"_"+region+"_"+syststr+"_"+cutstring
            query = "eos root://cmseos.fnal.gov/ ls "+logpath+" | grep "+fullstr
            files = subprocess.check_output(query,shell=True).decode(sys.stdout.encoding).split()
            
            exp = [x+"_"+fullstr+".root" for x in samples]
            miss = sorted(list(set(exp) - set(files)))
            
            if len(miss) > 0:
                print("      Missing the following:" )
                for m in miss:
                    print("        ",m)
            else:
                print("Not Missing anything!")

