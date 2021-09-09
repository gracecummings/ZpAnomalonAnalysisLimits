import ROOT
import glob
import os
import sys
#import gecorg as go
import gecorg_test as go
import numpy as np
import pandas as pd
import configparser

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
        vals = [x*processes["syst"][syst]["unc"] for x in processes["syst"][syst]["proc"]]
        sval = [x if x != 0 else "-" for x in vals]
        cardstr = "{0} {1} {2} {3} {4} {5}\n".format(str(syst),processes["syst"][syst]["type"],sval[0],sval[1],sval[2],sval[3])
        card.write(cardstr)
    #After here, put some sort of interation through a list of systematics
    card.close()
    

if __name__=='__main__':

    #will replace with command line options
    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'
    chan    = 'mumu'
    sigxs   = 1.0

    #load in the files with the hists
    bkgs        = go.backgrounds('BkgInputs/',zptcut,hptcut,metcut,btagwp)
    systbkgsup  = go.backgrounds('BkgJECSyst/',zptcut,hptcut,metcut,btagwp,"systjecup")
    systbkgsdwn = go.backgrounds('BkgJECSyst/',zptcut,hptcut,metcut,btagwp,"systjecdwn")
    sig         = go.signal('SignalInputs/',zptcut,hptcut,metcut,btagwp,sigxs,101.27)
    systsigup   = go.signal('SignalJECSyst/',zptcut,hptcut,metcut,btagwp,sigxs,101.27,"systjecup")
    systsigdwn  = go.signal('SignalJECSyst/',zptcut,hptcut,metcut,btagwp,sigxs,101.27,"systjecdwn")
    dyEst = ROOT.TFile('BkgInputs/Run2_2017_2018_dy_extraploation_Zptcut150.0_Hptcut300.0_metcut200.0_btagwp0.8.root')


    ####Prepping holders####
    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty1 = empty.Clone()
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    empty4 = empty.Clone()
    empty5 = empty.Clone()
    empty6 = empty.Clone()
    empty7 = empty.Clone()
    empty8 = empty.Clone()
    empty9 = empty.Clone()

    ####Gathering the Systematics####
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

    #Rename and Restructure
    httup.SetName("TT_jecUp")
    httdwn.SetName("TT_jecDown")
    hvvup.SetName("VV_jecUp")
    hvvdwn.SetName("VV_jecDown")
    
    
    ####Getting the Estimations####
    hdat = empty.Clone()
    htt = bkgs.getAddedHist(empty1,"TT","sr","h_zp_jigm")
    hzz  = bkgs.getAddedHist(empty2,"ZZTo2L2Q","sr","h_zp_jigm")
    hwz  = bkgs.getAddedHist(empty3,"WZTo2L2Q","sr","h_zp_jigm")
    hvv  = hzz.Clone()
    hvv.Add(hwz)

    #Get DY from estimation
    hdy = dyEst.Get("extrphist").Clone()

    #Rename and restructure background if necessary
    htt.SetName("TT")
    hvv.SetName("VV")
    hdy.SetName("DY")
    hdat.SetName("data_obs")

    #Signal Samples
    siginfo = sig.prepsigsr
    sigJECup = systsigup.prepsigsr
    sigJECdwn = systsigdwn.prepsigsr
    for s,sig in enumerate(siginfo):
        name = sig["name"]
        signame = "holder"
        if name != sigJECup[s]["name"] !=  sigJECdwn[s]["name"]:
            print("Cannot find matching entries in JEC systematics collections, aborting")
            break
        if "Tune" in name:
            strippedname = name.split("_Tune")[0]
            signame = strippedname.replace("-","")
        else:
            signame = name.replace("-","")

        #nominal    
        hsigori = sig["tfile"].Get("h_zp_jigm")
        #hsigori.Sumw2(ROOT.kTRUE)#Throws a warning that it is already created
        hsig = hsigori.Clone()
        hsig.SetName(signame)
        hsig.Scale(sig["scale"])

        #up
        hsigupori = sigJECup[s]["tfile"].Get("h_zp_jigm")
        hsigup = hsigupori.Clone()
        hsigup.SetName(signame+"_jecUp")
        hsigup.Scale(sigJECup[s]["scale"])

        #down
        hsigdwnori = sigJECdwn[s]["tfile"].Get("h_zp_jigm")
        hsigdwn = hsigdwnori.Clone()
        hsigdwn.SetName(signame+"_jecDown")
        hsigdwn.Scale(sigJECdwn[s]["scale"])

        #This is need to get the errors on the bins correct
        #built for nominal case
        #but not needed for Combine
        #for ibin in range(hsig.GetNbinsX()+1):
        # oribin = hsigori.GetBinContent(ibin)
        # orierr = hsigori.GetBinError(ibin)
        # newbinval = oribin*sig["scale"]
        # hsig.SetBinContent(ibin,newbinval)
        # hsig.SetBinError(ibin,newbinval**(1/2))

        #For writing the datacard
        #it makes sense to have the line be the key,
        #and then to have subdicts with the channel
        procdict = {"processnames":[signame,"DY","TT","VV"],
                    "hists":[hsig,hdy,htt,hvv],
                    "method":["mc","alpha","mc","mc"],
                    "syst":{"lumi_13TeV":
                            {"type":"lnN","unc":1.018,"proc":[1,1,1,1]},
                            #"alphar_alt":
                            #{"type":"shape","unc":1,"proc":[0,1,0,0]},
                            "jec":
                            {"type":"shapeN2","unc":1,"proc":[1,1,1,1]},
                            },
                    }
    
        prepRootName = go.makeOutFile('Run2_2017_2018_ZllHbbMET',chan+'_'+signame,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        prepRootFile = ROOT.TFile(prepRootName,"recreate")
        htt.Write()
        hvv.Write()
        hdy.Write()
        hdat.Write()
        hsig.Write()
        hsigup.Write()
        hsigdwn.Write()
        httup.Write()
        httdwn.Write()
        hvvup.Write()
        hvvdwn.Write()

        prepRootFile.Close()

        writeDataCard(procdict,prepRootName,chan)

    
    
