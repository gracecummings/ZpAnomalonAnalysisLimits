import ROOT
import glob
import argparse
import os
import sys
import numpy as np
#import gecorg_test as go
from datetime import date

def sampleType(sampstring,givejecs=False):
    #Make numerical code for type of sample
    if "Run" in sampstring:
        if "SingleMuon" in sampstring:
            samptype = -1
        elif ("SingleElectron" in sampstring) or ("EGamma" in sampstring):
            samptype = -2
        else:
            samptype = 0
    elif "ZpAnomalon" in sampstring:
        samptype = 1
    elif "DYJetsToLL" in sampstring:
         samptype = 2
    elif "TTTo" in sampstring:
        samptype = 3
    elif "WZTo" in sampstring:
        samptype = 4
    elif "ZZTo" in sampstring:
        samptype = 5
    elif "WWTo" in sampstring:
        samptype = 6

    else:
        samptype = -1000

    if "2018" in sampstring:
        year = 18
        if givejecs:
            if "Run2018A" in sampstring:
                year = 180
            if "Run2018B" in sampstring:
                year = 181
            if "Run2018C" in sampstring:
                year = 182
            if "Run2018D" in sampstring:
                year = 183
    if "Autumn18" in sampstring:
        year = 18
    if "2017" in sampstring:
        year = 17
        if givejecs:
            if "Run2017B" in sampstring:
                year = 170
            if "Run2017C" in sampstring:
                year = 171
            if "Run2017D" in sampstring:
                year = 172
            if "Run2017E" in sampstring:
                year = 172
            if "Run2017F" in sampstring:
                year = 173
    if "Fall17" in sampstring:
        year = 17
    if "2016" in sampstring:
        year = 16
        if givejecs:
            if "Run2016B" in sampstring:
                year = 160
            if "Run2016C" in sampstring:
                year = 160
            if "Run2016D" in sampstring:
                year = 160
            if "Run2016E" in sampstring: 
                year = 161
            if "Run2016F" in sampstring:
                year = 161
            if "Run2016G" in sampstring:
                year = 162
            if "Run2016H" in sampstring:
                year = 162
    if "Summer16" in sampstring:
        year = 16
    return samptype,year

def makeOutFile(sampstring,descrip,ftype,zptcut,hptcut,metcut,btagwp):
    if not os.path.exists("analysis_output_ZpAnomalon/"+str(date.today())+"/"):
        os.makedirs("analysis_output_ZpAnomalon/"+str(date.today())+"/")
    outFile = "analysis_output_ZpAnomalon/"+str(date.today())+"/"+sampstring+"_"+descrip+"_Zptcut"+zptcut+"_Hptcut"+hptcut+"_metcut"+metcut+"_btagwp"+btagwp+ftype
    return outFile

def channelEncoding(string):
    #channel code comes from a 3 bit binary
    #there are 3 possible canidates in an event
    #z->mumu, z->ee, z->emu
    #the number in decimal is this flag
    #ex. Z->mumu, no others: 100, or 4
    message = "Who the hell sent this command? Add a channel."
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
    parser.add_argument("-upjec","--upjecsystematics",type=bool,help="if you want jec up systematics output")
    parser.add_argument("-dwnjec","--downjecsystematics",type=bool,help="if you want jec down systematics output")
    parser.add_argument("-upjer","--upjersystematics",type=bool,help="if you want jer up systematics output")
    parser.add_argument("-dwnjer","--downjersystematics",type=bool,help="if you want jer down systematics output")
    parser.add_argument("-upuncl","--upunclsystematics",type=bool,help="if you want unclustered met  up systematics output")
    parser.add_argument("-dwnuncl","--downunclsystematics",type=bool,help="if you want unclustered met down systematics output")
    args = parser.parse_args() 
    samp = args.sample
    samptype = -1000
    extrajeccrap = True
    #Check what you are working with
    samptype,checkedyear = sampleType(samp,extrajeccrap)
    channel,channelprint = channelEncoding(args.channel)

    topyear = checkedyear #has era encoding

    #This is account for era encoding
    #Eras add a digit to the year, i.e, 18A -> 180
    if checkedyear > 100:
        tmp = checkedyear/10
        checkedyear = round(tmp)
    year = "20"+str(checkedyear)
    
    if samptype < -10:
        print "You have a problem, we do not undertand the sample coding"
        sys.exit()
    if channel <= 0:
        print "You have a problem, no channel given"
        sys.exit()
    origevnts = 0

    #First set of print statements
    print "Making topiary of ",samp
    print "     Sample type ",samptype
    print "     Sample Year ",year
    print "     Topiary Year ",topyear #This is a debug
    print "    ",channelprint

    #Prepare your TChain
    if samptype != 1 and ".root" not in samp:
    #if ".root" not in samp:#should also take multi file signal
        #for non-signal samples
        inChain = ROOT.TChain("PreSelection")
        #inputs  = glob.glob("../dataHandling/"+year+"/"+samp+"*.root")
        inputs  = glob.glob("../dataHandling/"+year+"_new/"+samp+"*.root")
        #inputs  = glob.glob("../dataHandling/"+year+"_old/"+samp+"*.root")
        for f in inputs:
            #print(f)
            inChain.Add(f)
            tf = ROOT.TFile.Open(f)
            origevnts += tf.Get("hnevents").GetBinContent(1)
    elif ".root" in samp:
        #for debug, and ntuple in working directory
        #inChain = ROOT.TChain("TreeMaker2/PreSelection")#ntuple
        inChain = ROOT.TChain("PreSelection")#skim
        inChain.Add(samp)
    elif samptype == 1 and ".root" not in samp and ".txt" not in samp:
        inChain = ROOT.TChain("TreeMaker2/PreSelection")
        #inputs  = glob.glob("../dataHandling/"+year+"/"+samp+"*.root")
        inputs  = glob.glob("../dataHandling/"+year+"_old/"+samp+"*.root")
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

    if ".txt" in samp:
        #for EOS signal sample
        #print("Ami eikhane")
        inChain = ROOT.TChain("PreSelection")
        with open(samp) as ftxt:
            lines = ftxt.readlines()
            inputs = [i.strip() for i in lines if i.startswith("root")]
        print(inputs)
        for f in inputs:
            inChain.Add(f)
            tf = ROOT.TFile.Open(f)
            origevnts = inChain.GetEntries()


    #Systematics?
    syststring  = "systnominal"
    systjerind  = 0
    systjecind  = 0
    systunclind = 0
    if args.upjecsystematics:
        syststring = "systjecup"
        systjecind = 1
    if args.downjecsystematics:
        syststring = "systjecdwn"
        systjecind = -1
    if args.upjersystematics:
        syststring = "systjerup"
        systjerind = 1
    if args.downjersystematics:
        syststring = "systjerdwn"
        systjerind = -1
    if args.upunclsystematics:
        syststring = "systunclup"
        systunclind = 1
    if args.downunclsystematics:
        syststring = "systuncldwn"
        systunclind = -1

    sysvec = ROOT.TVector(6)
    sysvec[0] = systjerind
    sysvec[1] = systjecind
    sysvec[5] = systunclind

    if sysvec.Norm1() > 1.0:
        print "Too many systematic flags at one time, stopping this maddness"
        exit()

        
    #outFile = go.makeOutFile(samp,'topiary_'+args.channel+'_'+syststring,'.root','No','Req','On','Reco')#Needs to become dynamic with cuts
    #outFile = go.makeOutFile(samp,'topiary_'+args.channel+'_'+syststring,'.root','0.0','250.0','0.0','0.0')#normal name
    outFile = makeOutFile(samp,'topiary_'+args.channel+'_'+syststring,'.root','0.0','0.0','0.0','0.0')#nemu

    print "     Systematics ",syststring
    print "     Events in TChain: ",inChain.GetEntries()
    print "     "+("Original data set had {0} events in type.").format(origevnts)
    print "    Saving topiary in ",outFile

    ROOT.gSystem.Load("RestFrames/lib/libRestFrames.so")
    ROOT.gSystem.Load("../UHH2/JetMETObjects/obj/libSUHH2JetMETObjects.so")
    ROOT.gROOT.ProcessLine(".include ../UHH2/JetMETObjects/interface")
    ROOT.gSystem.Load("TreeMakerTopiary.so")#having problems grabbing restframes
    ROOT.gInterpreter.Declare('#include "TreeMakerTopiary.h"')

    topiary = ROOT.TreeMakerTopiary(inChain,samptype,topyear,channel,sysvec)#Trying to access something that is not there
    topiary.Loop(outFile,origevnts,samptype,topyear,channel,sysvec)



    


