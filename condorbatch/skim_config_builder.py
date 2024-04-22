import argparse
import os

if __name__=='__main__':
    infile = 'contents_ZpAnomalon_OffSig_2018_2021-08-07.txt'
    f = open(infile)
    lines = f.readlines()
    samples = [x.split('-madgraph-')[0] for x in lines]
    names = list(set(samples))

    for name in names:
        skimconf = "samplesToSkim_"+name+"_ZpAnomalon_OffSig_2018_2021-08-07.txt"
        cmd = "eos root://cmseos.fnal.gov ls /store/group/lpcboostres/ZpAnomalon_OffSig_2018_2021-08-07 | grep "+name+" > "+skimconf
        os.system(cmd)
    
    
