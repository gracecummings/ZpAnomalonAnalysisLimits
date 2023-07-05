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
    samptype,checkedyear = go.sampleType(samp,extrajeccrap)
    channel,channelprint = channelEncoding(args.channel)

    topyear = checkedyear #has era encoding

    #This is account for era encoding
    #Eras add a digit to the year, i.e, 18A -> 180
    if checkedyear > 100:
        tmp = checkedyear/10
        checkedyear = round(tmp)
    year = "20"+str(checkedyear)
    
    if samptype < -10:
        print("You have a problem, we do not undertand the sample coding")
        sys.exit()
    if channel <= 0:
        print("You have a problem, no channel given")
        sys.exit()
    origevnts = 0
    trimmedevents = 0
    passedfilter = 0
    passedzcandnum = 0
    passedfat = 0
    passedfilandfat = 0
    passedfilandz = 0

    #First set of print statements
    print( "Making topiary of ",samp)
    print("     Sample type ",samptype)
    print("     Sample Year ",year)
    print("     Topiary Year ",topyear)#This is a debug
    print("    ",channelprint)

    #Prepare your TChain
    if samptype != 1 and ".root" not in samp:
    #if ".root" not in samp:#should also take multi file signal
        #for non-signal samples
        inChain = ROOT.TChain("PreSelection")
        #inputs  = glob.glob("../dataHandling/"+year+"/"+samp+"*.root")
        #inputs  = glob.glob("../dataHandling/"+year+"_new/"+samp+"*.root")
        inputs  = glob.glob("../dataHandling/"+year+"_newee/"+samp+"*.root")
        #inputs  = glob.glob("../dataHandling/"+year+"_newee/Autumn18_DYJetsToLL_2022-06-17/"+samp+"*.root")o
        #inputs  = glob.glob("../dataHandling/"+year+"_old/"+samp+"*.root")
        #inputs  = glob.glob("../dataHandling/"+year+"_newpu/"+samp+"*.root")
        for f in inputs:
            #print(f)
            inChain.Add(f)
            tf = ROOT.TFile.Open(f)
            origevnts += tf.Get("hnevents").GetBinContent(1)
            #trimmedevents += tf.Get("hpassall").GetBinContent(1)
            #passedfilter += tf.Get("hpassfilter").GetBinContent(1)
            #passedzcandnum += tf.Get("hpassZ").GetBinContent(1)
            #passedfat += tf.Get("hpassfat").GetBinContent(1)
            #passedfilandfat += tf.Get("hpassfilandfat").GetBinContent(1)
            #passedfilandz += tf.Get("hpassfilandz").GetBinContent(1)
    elif ".root" in samp:
        #for debug, and ntuple in working directory
        inChain = ROOT.TChain("TreeMaker2/PreSelection")
        inChain.Add(samp)
    elif samptype == 1 and ".root" not in samp:
        if "2018" not in year:
            inChain = ROOT.TChain("PreSelection")
            inputs  = glob.glob("../dataHandling/"+year+"_newee/"+samp+"*.root")
        else:
            inChain = ROOT.TChain("TreeMaker2/PreSelection")
            inputs  = glob.glob("../dataHandling/"+year+"_old/"+samp+"*.root")
        #inputs   = glob.glob("../dataHandling/2018_data_muonunc/Autumn18_ZpAnomalonHZ_2023-01-26/"+samp+"*.root")
        #PU debug
        #inChain = ROOT.TChain("PreSelection")
        #inputs  = glob.glob("../dataHandling/"+year+"_newpu/"+samp+"*.root")
        #
        #print(inputs)
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
        print("Too many systematic flags at one time, stopping this maddness")
        exit()
        

    if (origevnts < 1):
        print('      ***Could not get events from the trims')
        checkevnts = 0        
        if "MiniAOD" in inputs[0]:
            print('      ***  getting number from file name, borked skim')
            for f in inputs:
                nproc = int((f.split("MiniAOD")[0]).split('skim_')[-1])
                checkevnts +=nproc
            origevnts = checkevnts
        elif "RA" in inputs[0]:
            print('      ***  getting number from tree, it is the ntuple')
            origevnts = inChain.GetEntries()
        else:
            "cannot figure out what went wrong with the counting, so something is wrong"
            
        
    #outFile = go.makeOutFile(samp,'topiary_'+args.channel+'_'+syststring,'.root','No','Req','On','Reco')#Needs to become dynamic with cuts
    outFile = go.makeOutFile(samp,'topiary_'+args.channel+'_'+syststring,'.root','0.0','250.0','0.0','0.0')#normal name
    #outFile = go.makeOutFile(samp,'topiary_'+args.channel+'_'+syststring,'.root','0.0','0.0','0.0','0.0')#nemu
    print("     Systematics ",syststring)
    print("     Events in TChain: ",inChain.GetEntries())
    print("     Number of files in the TChain: ",len(inputs))
    print(("     Original data set had {0} events in type.").format(origevnts))
    print("    Saving topiary in ",outFile)
    #print("    trimmed events ",trimmedevents)
    #print("    passed the filters ",passedfilter)
    #print("    passed the Z number ",passedzcandnum)
    #print("    passed the fat number ",passedfat)
    #print("    passed filters and fat ",passedfilandfat)
    #print("    passed filters and Z ",passedfilandz)


    ROOT.gSystem.Load("../UHH2/JetMETObjects/obj/libSUHH2JetMETObjects.so")
    ROOT.gROOT.ProcessLine(".include ../UHH2/JetMETObjects/interface")
    ROOT.gSystem.Load("TreeMakerTopiary.so")
    ROOT.gInterpreter.Declare('#include "TreeMakerTopiary.h"')

    topiary = ROOT.TreeMakerTopiary(inChain,samptype,topyear,channel,sysvec)
    topiary.Loop(outFile,origevnts,samptype,topyear,channel,sysvec)



    


