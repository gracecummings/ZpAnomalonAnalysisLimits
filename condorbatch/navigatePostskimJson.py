import json
import sys
import os
import argparse
from datetime import date

#def getSampleName(sampstr):
#    name = "badsamp"
#    if "WZ" in sampstr:
#        name = sampstr.split(")#

parser = argparse.ArgumentParser()

if __name__=='__main__':
    parser.add_argument("-j","--sampleJson",help=".json file with the skims to be made into topiary")
    args = parser.parse_args()

    #Check json
    if (args.sampleJson is None):
        print("No list of skims provided, please provide appropriate json")
        fsjson = {}
    else:
        if (not os.path.exists(args.sampleJson)):
            print("Invalid json")
        fstotop = open(args.sampleJson,'r')
        fsjson = json.load(fstotop)

    #Print info
    print(fsjson.keys())
    gensets = list(set([f.split("_topiary")[0] for f in fsjson.keys()]))
    print(gensets)
    
    #Get full list grouped by sample
    sampgrouped = {}
    for samp in gensets:
        sampgrouped[samp] = [fsjson[f]["path"]+"/"+fsjson[f]["topiaries"][0] for f in fsjson.keys() if samp in f]
        #sampgrouped[samp] = [f for f in fsjson.keys() if samp in f]
        print(samp)
        print(sampgrouped[samp])
    
    #if len(fsjson.keys()) > 0:
    #    for key in fsjson.keys():
            #print(key)
            #print("   ",fsjson[key]['path'])
            #for f in fsjson[key]['topiaries']:
            #    print("          ",f)
