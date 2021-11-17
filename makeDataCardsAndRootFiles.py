import ROOT
import glob
import os
import sys
import gecorg_test as go
import numpy as np
import pandas as pd
import configparser
import argparse

def writeDataCard(processes,rootFileName,channel):
    signame = processes["processnames"][0]
    nbkg = len(processes["processnames"])-1
    obsEvents = 0
    binname = signame+"_"+channel
    namestr = " ".join(processes["processnames"])
    rates   = [str(hist.Integral()) for hist in processes["hists"]]
    prepCardName = go.makeOutFile('Run2_2017_2018_ZllHbbMET','datacard_'+chan+'_'+signame,'.txt',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    card = open(prepCardName,"w")

    #Write the card
    card.write("imax 1\n")
    card.write("jmax {0}\n".format(nbkg))
    card.write("kmax *\n")              
    card.write("------------\n")
    card.write("shapes * * {0} $PROCESS $PROCESS_$SYSTEMATIC \n".format(rootFileName.split("/")[-1]))
    card.write("------------\n")
    card.write("bin {0} \n".format(binname))
    card.write("observation {0} \n".format(obsEvents))
    card.write("------------\n")
    card.write("bin {0} {1} {2} {3}\n".format(binname,binname,binname,binname))
    card.write("process "+namestr+"\n")
    card.write("process 0 1 2 3\n")#hardcode
    card.write("rate "+" ".join(rates)+"\n")
    card.write("------------\n")

    for syst in processes["syst"].keys():
        #vals = [x*processes["syst"][syst]["unc"] for x in processes["syst"][syst]["proc"]]
        vals = processes["syst"][syst]["proc"]
        sval = [x if x != 0 else "-" for x in vals]
        cardstr = "{0} {1} {2} {3} {4} {5}\n".format(str(syst),processes["syst"][syst]["type"],sval[0],sval[1],sval[2],sval[3])
        card.write(cardstr)

    card.write("* autoMCStats 1\n")
    card.close()

#def gatherSystematics(confsecsyst):
    
    

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    args = parser.parse_args()
    
    zptcut  = str(args.zptcut)#'150.0'
    hptcut  = str(args.hptcut)#'300.0'
    metcut  = str(args.metcut)#'200.0'
    btagwp  = str(args.btagwp)#'0.8'

    #zptcut  = '150.0'
    #hptcut  = '300.0'
    #metcut  = '200.0'
    #btagwp  = '0.8'

    chan    = 'mumu'
    sigxs   = 1.0
    rebindiv = 2

    limrangelow = 1400
    limrangehigh = 3000

    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()

    ####load in the files with the nominal distributions
    bkgs = go.backgrounds(config.get('nominal','pathnom'),zptcut,hptcut,metcut,btagwp,config.get('nominal','strnom'))
    sig  = go.signal(config.get('nominal','pathsignom'),zptcut,hptcut,metcut,btagwp,sigxs,101.27,config.get('nominal','strnom'))
    dyEst = ROOT.TFile(config.get('nominal','pathnom')+'Run2_2017_2018_dy_extraploation'+config.get('nominal','strnom')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')

    ####Prepping holders####
    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty1 = empty.Clone()
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    
    ####Getting the Estimations####
    hdat = empty.Clone()
    hdy    = dyEst.Get("extrphist").Clone()
    htt = bkgs.getAddedHist(empty1,"TT","sr","h_zp_jigm")
    hzz  = bkgs.getAddedHist(empty2,"ZZTo2L2Q","sr","h_zp_jigm")
    hwz  = bkgs.getAddedHist(empty3,"WZTo2L2Q","sr","h_zp_jigm")
    hvv  = hzz.Clone()
    hvv.Add(hwz)

    ####Rename and restucture
    htt.SetName("TT")
    hvv.SetName("VV")
    hdy.SetName("DY")
    hdat.SetName("data_obs")
    htt.Rebin(rebindiv)
    hvv.Rebin(rebindiv)
    hdat.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
    htt.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
    hvv.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
    hdy.GetXaxis().SetRangeUser(limrangelow,limrangehigh)

    ####For Each signal, make a datacard, and a root file with all systematics
    siginfo = sig.prepsigsr
    sigcolors = go.colsFromPalette(siginfo,ROOT.kCMYK)
    for s,sig in enumerate(siginfo):
        name = sig["name"]
        signame = "holder"
        if "Tune" in name:
            strippedname = name.split("_Tune")[0]
            signame = strippedname.replace("-","")
        else:
            signame = name.replace("-","")
        print("------- Looking at signal sample ",signame)

        ####Make Files
        prepRootName = go.makeOutFile('Run2_2017_2018_ZllHbbMET',chan+'_'+signame,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        prepRootFile = ROOT.TFile(prepRootName,"recreate")
            
        ####Nominal Signal
        hsigori = sig["tfile"].Get("h_zp_jigm")
        hsigori.Sumw2(ROOT.kTRUE)#Throws a warning that it is already created

        #####Debug of bin uncertainties
        #print("Looking at signal errors and uncertainties before scaling")
        #for ibin in range(hsigori.GetNbinsX()+1):
        #    print(" Bin Content: ",hsigori.GetBinContent(ibin))
        #    print(" Bin Error: ",hsigori.GetBinError(ibin))

        hsig = hsigori.Clone()
        hsig.SetName(signame)
        hsig.Scale(sig["scale"])

       # print("Looking at signal errors and uncertainties after scaling")
       # for ibin in range(hsig.GetNbinsX()+1):
       #     print(" Bin Content: ",hsig.GetBinContent(ibin))
       #     print(" Bin Error: ",hsig.GetBinError(ibin))

        hsig.Rebin(rebindiv)
        hsig.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
        prepRootFile.cd() 
        htt.Write()
        hvv.Write()
        hdy.Write()
        hdat.Write()
        hsig.Write()

        #####Gather Systematics
        systdict = {"lumi_13TeV":{"type":"lnN","unc":1.018,"proc":[0,1.018,1.018,1.018]},#GOOD LUMI
                    }

        for syst in systs:
            if syst == 'nominal':
                continue
            print("------- Looking at systematic ",syst)

            
            systbkgsup  = go.backgrounds(config.get(syst,'pathup'),zptcut,hptcut,metcut,btagwp,config.get(syst,'strup'))
            systbkgsdwn = go.backgrounds(config.get(syst,'pathdwn'),zptcut,hptcut,metcut,btagwp,config.get(syst,'strdwn'))
            systsigup   = go.signal(config.get(syst,'pathsigup'),zptcut,hptcut,metcut,btagwp,sigxs,101.27,config.get(syst,'strup'))
            systsigdwn  = go.signal(config.get(syst,'pathsigdwn'),zptcut,hptcut,metcut,btagwp,sigxs,101.27,config.get(syst,'strdwn'))

            if len(systbkgsup.bkgs["DYJetsToLL"][18]["sb"][0]) < 1:
                print("        There are no DYJets entires for this systematic.")
                print("        moving on. This systematic will not included.")
                continue
            
            appcode = config.get(syst,'applist').split(',')
            applist = [float(x) for x in appcode]
            systdict[syst] = {"type":config.get(syst,'type'),"unc":1.0,"proc":applist}
            
            if rebindiv == 2:
                dyEstup     = ROOT.TFile(config.get(syst,'pathup')+'Run2_2017_2018_dy_extraploation'+config.get(syst,'strup')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
                dyEstdwn    = ROOT.TFile(config.get(syst,'pathdwn')+'Run2_2017_2018_dy_extraploation'+config.get(syst,'strdwn')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')

            ####Make it useful
            sigup  = systsigup.prepsigsr
            sigdwn = systsigdwn.prepsigsr
            if name != sigup[s]["name"] !=  sigdwn[s]["name"]:
                print("Cannot find matching entries in systematics collections, aborting")
                break

            ####Prepping holders####
            empty4 = empty.Clone()
            empty5 = empty.Clone()
            empty6 = empty.Clone()
            empty7 = empty.Clone()
            empty8 = empty.Clone()
            empty9 = empty.Clone()

            ####Gathering the Systematic
            #Background
            httup  = systbkgsup.getAddedHist(empty4,"TT","sr","h_zp_jigm")
            httdwn = systbkgsdwn.getAddedHist(empty5,"TT","sr","h_zp_jigm")
            hzzup  = systbkgsup.getAddedHist(empty6,"ZZTo2L2Q","sr","h_zp_jigm")
            hzzdwn = systbkgsdwn.getAddedHist(empty7,"ZZTo2L2Q","sr","h_zp_jigm")
            hwzup  = systbkgsup.getAddedHist(empty8,"WZTo2L2Q","sr","h_zp_jigm")
            hwzdwn = systbkgsdwn.getAddedHist(empty9,"WZTo2L2Q","sr","h_zp_jigm")
            hvvup  = hzzup.Clone()
            hvvdwn = hzzdwn.Clone()
            hvvup.Add(hwzup)
            hvvdwn.Add(hwzdwn)
            hdyup  = dyEstup.Get("extrphist").Clone()
            hdydwn = dyEstdwn.Get("extrphist").Clone()

            #Rename and Restructure
            httup.SetName("TT_"+syst+"Up")
            httdwn.SetName("TT_"+syst+"Down")
            hvvup.SetName("VV_"+syst+"Up")
            hvvdwn.SetName("VV_"+syst+"Down")
            hdyup.SetName("DY_"+syst+"Up")
            hdydwn.SetName("DY_"+syst+"Down")
            httup.Rebin(rebindiv)
            httdwn.Rebin(rebindiv)
            hvvup.Rebin(rebindiv)
            hvvdwn.Rebin(rebindiv)
            httup.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
            httdwn.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
            hvvup.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
            hvvdwn.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
            hdyup.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
            hdydwn.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
            
            #Signal
            hsigupori = sigup[s]["tfile"].Get("h_zp_jigm")
            hsigup = hsigupori.Clone()
            hsigup.SetName(signame+"_"+syst+"Up")
            hsigup.Scale(sigup[s]["scale"])
            hsigup.Rebin(rebindiv)
            hsigup.GetXaxis().SetRangeUser(limrangelow,limrangehigh)
            
            hsigdwnori = sigdwn[s]["tfile"].Get("h_zp_jigm")
            hsigdwn = hsigdwnori.Clone()
            hsigdwn.SetName(signame+"_"+syst+"Down")
            hsigdwn.Scale(sigdwn[s]["scale"])
            hsigdwn.Rebin(rebindiv)
            hsigdwn.GetXaxis().SetRangeUser(limrangelow,limrangehigh)

            #Write the histograms
            prepRootFile.cd()
            print("------- Writing the up/dwn histograms")
            hsigup.Write()
            hsigdwn.Write()
            httup.Write()
            httdwn.Write()
            hvvup.Write()
            hvvdwn.Write()
            hdyup.Write()
            hdydwn.Write()

            #This is need to get the errors on the bins correct
            #built for nominal case
            #but not needed for Combine
            #for ibin in range(hsig.GetNbinsX()+1):
            # oribin = hsigori.GetBinContent(ibin)
            # orierr = hsigori.GetBinError(ibin)
            # newbinval = oribin*sig["scale"]
            # hsig.SetBinContent(ibin,newbinval)
            # hsig.SetBinError(ibin,newbinval**(1/2))

            ##### Make Plot of each syst for each signal sample
            tc = ROOT.TCanvas("tc",signame,600,600)
            hsig.SetLineColor(ROOT.kBlack)
            hsigup.SetLineColor(ROOT.kRed)
            hsigdwn.SetLineColor(ROOT.kBlue)
            hsigup.SetStats(0)
            hsigup.GetXaxis().SetTitle("Z' Jigsaw Mass Estmator")
            tc.cd()
            hsigup.Draw('hist')
            hsig.Draw('histsame')
            hsigdwn.Draw('histsame')
            sigplotname = go.makeOutFile(signame,syst+'updowncomp','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
            tc.SaveAs(sigplotname)


        #For writing the datacard
        procdict = {"processnames":[signame,"DY","TT","VV"],
                    "hists":[hsig,hdy,htt,hvv],
                    "method":["mc","alpha","mc","mc"],
                    "syst":systdict,
        }
        print("------- Defined the Datacard Dict")
        print("------- Writing the Datacard")
        writeDataCard(procdict,prepRootName,chan)
        print("------- About to close the root file")
        prepRootFile.Close()

