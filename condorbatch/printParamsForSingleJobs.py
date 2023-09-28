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
    if len(fsjson.keys()) > 0:
        for key in fsjson.keys():
            for samp in fsjson[key]:
                sampleName = samp.split("_13TeV")[0]
                print("    ",sampleName)
                fullpaths = ["root://cmseos.fnal.gov/"+key+"/"+x for x in fsjson[key][samp]]
                samplistasstr = str(fullpaths).replace(' ','')
                print(samplistasstr)
