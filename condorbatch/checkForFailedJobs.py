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
    parser.add_argument("-r","--resub",type=bool,help="resubmit the jobs?")
    parser.add_argument("-c","--channel",type=str,help="channel: mumu,ee,emu")
    parser.add_argument("-syst","--syststr",type=str,help="the systematic flag you used, again, to check")
    args = parser.parse_args()

    #Check json
    if (args.sampleJson is None):
        print("No list of samples provided, please provide appropriate json")
        fsjson = {}
    else:
        if (not os.path.exists(args.sampleJson)):
            print("Invalid json")
        fstotop = open(args.sampleJson,'r')
        fsjson = json.load(fstotop)

    totalsamps = []
    if len(fsjson.keys()) > 0:
        for key in fsjson.keys():
            for samp in fsjson[key]:
                sampleName = samp.split("_13TeV")[0]
                #print("    ",sampleName)
                totalsamps.append(sampleName)

    #Get the files that were made
    systname = ""
    if args.syststr:
        systname = args.syststr
    eospath = "/store/user/lpcboostres/topiaries_systematics-"+systname+"_"+args.date
    fs = subprocess.check_output("eos root://cmseos.fnal.gov ls "+eospath,shell=True).decode(sys.stdout.encoding).split()
    fnames = [f.split('_topiary')[0] for f in fs]
    fnamesclean = list(set(fnames))


    #Find the missed ones 
    Ls = list(set(totalsamps) - set(fnamesclean))

    if len(Ls) > 0:
        print("You are missing the following topiaries")
        print(Ls)
        if args.resub:
            print("Resubmitting jobs")
            for samp in Ls:
                subprocess.run(["python","submitTopiaryJobs.py","-j",args.sampleJson,"-c",args.channel,"-syst",args.syststr,"-s",samp+"_13TeV"])
    else:
        print("No topiaries missing, ready for next process")

        
        


