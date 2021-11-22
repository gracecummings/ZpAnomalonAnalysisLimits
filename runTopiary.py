import ROOT
import glob
import argparse
import os
import sys
import numpy as np
import gecorg_test as go
from datetime import date

def channelEncoding(string):
    #channel code comes from a 3 bit binary
    #there are 3 possible canidates in an event
    #z->mumu, z->ee, z->emu
    #the number in decimal is this flag
    #ex. Z->mumu, no others: 100, or 4
    if "mumu" == string:
        channel = 4
        message = "Using the Z->mumu selections"
    elif "ee" == string:
        channel = 2
        message = "Using the Z->ee selections"
    elif "emu" == string:
        channel = 1
        message = "Usiing the ttbar background emu selections"
    else:
        channel = -1
    return channel,message

parser = argparse.ArgumentParser()

if __name__=="__main__":
    parser.add_argument("-s","--sample",help="sample name")
    parser.add_argument("-c","--channel",help="string for channel: mumu,ee, or emu")
    parser.add_argument("-upjec","--upsystematics",type=bool,help="if you want jec up systematics output")
    parser.add_argument("-dwnjec","--downsystematics",type=bool,help="if you want jec down systematics output")
    args = parser.parse_args() 
    samp = args.sample
    samptype = -1
    extrajeccrap = True
    #Check what you are working with
    samptype,checkedyear = go.sampleType(samp,extrajeccrap)
    channel,channelprint = channelEncoding(args.channel)

    topyear = checkedyear #has era encoding

    #This is account for era encoding
    #Eras add a digit to the year, i.e, 18A -> 180
    if checkedyear > 100:
        tmp = checkedyear/10
        checkedyear = round(tmp)
    year = "20"+str(checkedyear)
    
    if samptype < 0:
        print("You have a problem, we do not undertand the sample coding")
        sys.exit()
    if channel <= 0:
        print("You have a problem, no channel given")
        sys.exit()
    origevnts = 0

    #Prepare your TChain
    if samptype != 1 and ".root" not in samp:
    #if ".root" not in samp:#should also take multi file signal
        #for non-signal samples
        inChain = ROOT.TChain("PreSelection")
        inputs  = glob.glob("../dataHandling/"+year+"/"+samp+"*.root")
        for f in inputs:
            inChain.Add(f)
            tf = ROOT.TFile.Open(f)
            origevnts += tf.Get("hnevents").GetBinContent(1)
    elif ".root" in samp:
        #for debug, and ntuple in working directory
        inChain = ROOT.TChain("TreeMaker2/PreSelection")
        inChain.Add(samp)
    elif samptype == 1 and ".root" not in samp:
        inChain = ROOT.TChain("TreeMaker2/PreSelection")
        inputs  = glob.glob("../dataHandling/"+year+"/"+samp+"*.root")
        print(inputs)
        for f in inputs:
            inChain.Add(f)
            tf = ROOT.TFile.Open(f)
            #origevnts += tf.Get("hnevents").GetBinContent(1)
    else:
        #for signal
        inChain = ROOT.TChain("TreeMaker2/PreSelection")
        inputs = glob.glob("../dataHandling/"+year+"/"+samp+"*.root")
        inChain.Add("../dataHandling/"+year+"/"+samp+"*.root")
        origevnts = inChain.GetEntries()


    #Systematics?
    syststring = "systnominal"
    systind    = 0
    if args.upsystematics:
        syststring = "systjecup"
        systind = 1
    if args.downsystematics:
        syststring = "systjecdwn"
        systind = -1

    sysvec = ROOT.TVector(6)
    sysvec[0] = 0
    sysvec[1] = systind
    #sysvecfill = np.array([0,systind])
    #sysvec.Use(sysvecfill)
        
    outFile = go.makeOutFile(samp,'topiary_'+args.channel+'_'+syststring,'.root','0.0','250.0','0.0','0.0')#Needs to become dynamic with cuts
    print( "Making topiary of ",samp)
    print("     Sample type ",samptype)
    print("     Sample Year ",year)
    print("     Topiary Year ",topyear)
    print("    ",channelprint)
    print("     Systematics ",syststring)
    print("     Events in TChain: ",inChain.GetEntries())
    print(("     Original data set had {0} events in type.").format(origevnts))
    print("    Saving topiary in ",outFile)

    ROOT.gSystem.Load("../UHH2/JetMETObjects/obj/libSUHH2JetMETObjects.so")
    ROOT.gROOT.ProcessLine(".include ../UHH2/JetMETObjects/interface")
    ROOT.gSystem.Load("TreeMakerTopiary.so")
    ROOT.gInterpreter.Declare('#include "TreeMakerTopiary.h"')


    #topiary = ROOT.TreeMakerTopiary(inChain,samptype,topyear,channel,systind)
    #topiary.Loop(outFile,origevnts,samptype,topyear,channel,systind)

    topiary = ROOT.TreeMakerTopiary(inChain,samptype,topyear,channel,sysvec)
    topiary.Loop(outFile,origevnts,samptype,topyear,channel,sysvec)



    


