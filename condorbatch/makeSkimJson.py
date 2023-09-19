import subprocess
import sys
import json


if __name__=="__main__":
    dirs = open('skim_directories.txt')

    dirnames = dirs.readlines()
    dirnames = list(map(lambda x : x.split("\n")[0],dirnames))
    paths = ['/store/group/lpcboostres/'+x for x in dirnames]
    print(paths)

    samples = {}
    for path in paths:
        print("Checking ............................. ",path)
        fs = subprocess.check_output("eos root://cmseos.fnal.gov ls "+path,shell=True).decode(sys.stdout.encoding).split()
        fnames = [f.split('pythia')[0] for f in fs]
        fnamesclean = list(set(fnames))
        samps = {}
        for name in fnamesclean:
            samps[name] = [f for f in fs if name in f]
            print(samps[name])
        samples[path] = samps

    with open('skim_locations.json','w') as js:
        json.dump(samples,js)

