import subprocess
import sys
import json
import argparse

parser = argparse.ArgumentParser()

if __name__=="__main__":
    parser.add_argument("-f","--txtfile",help="text file with list of directories in lpcboostres")
    parser.add_argument("-n","--namein",type=str,help="name to append to file")
    args = parser.parse_args()

    dirs = open(args.txtfile)

    dirnames = dirs.readlines()
    dirnames = list(map(lambda x : x.split("\n")[0],dirnames))
    paths = ['/store/group/lpcboostres/'+x for x in dirnames]
    print(paths)

    samples = {}
    for path in paths:
        print("Checking ............................. ",path)
        fs = subprocess.check_output("eos root://cmseos.fnal.gov ls "+path,shell=True).decode(sys.stdout.encoding).split()
        if "Run20" not in path:
            fnames = [f.split('pythia')[0] for f in fs]
        else:
            fnames = [f.split('_')[0] for f in fs]
        fnamesclean = list(set(fnames))
        samps = {}
        for name in fnamesclean:
            print("   Adding list of files for: ",name)
            samps[name] = [f for f in fs if name in f]
        samples[path] = samps

    jsonname = "skim_locations_"+args.namein+".json"
    with open(jsonname,'w') as js:
        json.dump(samples,js)

