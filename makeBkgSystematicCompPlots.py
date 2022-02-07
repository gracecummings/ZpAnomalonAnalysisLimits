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
        #vals = [x*processes["syst"][syst]["unc"] for x in processes["syst"][syst]["proc"]]
        vals = processes["syst"][syst]["proc"]
        sval = [x if x != 0 else "-" for x in vals]
        cardstr = "{0} {1} {2} {3} {4} {5}\n".format(str(syst),processes["syst"][syst]["type"],sval[0],sval[1],sval[2],sval[3])
        card.write(cardstr)
    #After here, put some sort of interation through a list of systematics
    card.close()

#def gatherSystematics(confsecsyst):
    
    

if __name__=='__main__':

    #will replace with command line options
    zptcut  = '100.0'
    hptcut  = '300.0'
    metcut  = '75.0'
    btagwp  = '0.8'
    chan    = 'mumu'
    sigxs   = 1.0
    rebindiv = 2

    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()

    ####load in the files with the nominal distributions
    bkgs = go.backgrounds(config.get('nominal','pathnom'),zptcut,hptcut,metcut,btagwp,config.get('nominal','strnom'))
    sig  = go.signal(config.get('nominal','pathsignom'),zptcut,hptcut,metcut,btagwp,sigxs,systr=config.get('nominal','strnom'))
    dyEst = ROOT.TFile(config.get('nominal','pathnom')+'/Run2_161718_dy_extraploation'+config.get('nominal','strnom')+'_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')

    ####Prepping holders####
    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty1 = empty.Clone()
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    empty41 = empty.Clone()
    empty42 = empty.Clone()
    empty43 = empty.Clone()
    empty44 = empty.Clone()
    empty45 = empty.Clone()
    
    ####Getting the Estimations####
    hdat = empty.Clone()
    hdy    = dyEst.Get("extrphist").Clone()
    htt = bkgs.getAddedHist(empty1,"TT","sr","h_zp_jigm")
    hzz  = bkgs.getAddedHist(empty2,"ZZTo2L2Q","sr","h_zp_jigm")
    hwz  = bkgs.getAddedHist(empty3,"WZTo2L2Q","sr","h_zp_jigm")
    hdytr= bkgs.getAddedHist(empty41,"DYJetsToLL","sr","h_zp_jigm")
    hvv  = hzz.Clone()
    hvv.Add(hwz)

    ####Troubleshooting
    httchecksb = bkgs.getAddedHist(empty42,"TT","sb","h_zp_jigm")
    hzzchecksb  = bkgs.getAddedHist(empty43,"ZZTo2L2Q","sb","h_zp_jigm")
    hwzchecksb  = bkgs.getAddedHist(empty44,"WZTo2L2Q","sb","h_zp_jigm")
    hdytrchecksb= bkgs.getAddedHist(empty45,"DYJetsToLL","sb","h_zp_jigm")

    #print(bkgs.bkgs["DYJetsToLL"][17]["sb"][0])
    #print(len(bkgs.bkgs["DYJetsToLL"][17]["sb"][0]))

    print("The SR yields, no cut on M_Zp")
    #print("   dyest: ",hdy.Integral())
    print("  dyhist: ",hdytr.Integral())
    print("   ttbar: ",htt.Integral())
    print("      zz: ",hzz.Integral())
    print("      wz: ",hwz.Integral())

    print("The SB yields, no cut on M_Zp")
    print("      dy: ",hdytrchecksb.Integral())
    print("   ttbar: ",httchecksb.Integral())
    print("      zz: ",hzzchecksb.Integral())
    print("      wz: ",hwzchecksb.Integral())

    ####Rename and restucture
    htt.SetName("TT")
    hvv.SetName("VV")
    hdy.SetName("DY")
    hdat.SetName("data_obs")
    htt.Rebin(rebindiv)
    hvv.Rebin(rebindiv)
    hdat.GetXaxis().SetRangeUser(1500,5000)
    htt.GetXaxis().SetRangeUser(1500,5000)
    hvv.GetXaxis().SetRangeUser(1500,5000)
    hdy.GetXaxis().SetRangeUser(1500,5000)

    for syst in systs:
        if "nominal" == syst:
            continue
        print("------- Looking at systematic ",syst)
        systbkgsup  = go.backgrounds(config.get(syst,'pathup'),zptcut,hptcut,metcut,btagwp,config.get(syst,'strup'))
        systbkgsdwn = go.backgrounds(config.get(syst,'pathdwn'),zptcut,hptcut,metcut,btagwp,config.get(syst,'strdwn'))

        if rebindiv == 2:
            dyEstup     = ROOT.TFile(config.get(syst,'pathup')+'/Run2_161718_dy_extraploation'+config.get(syst,'strup')+'_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
            dyEstdwn    = ROOT.TFile(config.get(syst,'pathdwn')+'/Run2_161718_dy_extraploation'+config.get(syst,'strdwn')+'_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')


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
        httup.GetXaxis().SetRangeUser(1500,5000)
        httdwn.GetXaxis().SetRangeUser(1500,5000)
        hvvup.GetXaxis().SetRangeUser(1500,5000)
        hvvdwn.GetXaxis().SetRangeUser(1500,5000)
        hdyup.GetXaxis().SetRangeUser(1500,5000)
        hdydwn.GetXaxis().SetRangeUser(1500,5000)
        
        
        #This is need to get the errors on the bins correct
        #built for nominal case
        #but not needed for Combine
        #for ibin in range(hsig.GetNbinsX()+1):
        # oribin = hsigori.GetBinContent(ibin)
        # orierr = hsigori.GetBinError(ibin)
        # newbinval = oribin*sig["scale"]
        # hsig.SetBinContent(ibin,newbinval)
        # hsig.SetBinError(ibin,newbinval**(1/2))
        
        #Debug for background
        #Extra plots
        dynormup = np.load(config.get(syst,'pathup')+'/Run2_161718_dynormalization_'+config.get(syst,'strup')+'_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
        dynormnom = np.load(config.get(syst,'pathnom')+'/Run2_161718_dynormalization_'+config.get(syst,'strnom')+'_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
        dynormdwn = np.load(config.get(syst,'pathdwn')+'/Run2_161718_dynormalization_'+config.get(syst,'strdwn')+'_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]

        print("DY norm up:  ",dynormup)
        print("DY norm nom: ",dynormnom)
        print("DY norm dwn: ",dynormdwn)
        emptypt = tf1.Get('h_h_pt')
        emptypt.Reset("ICESM")
        empty10 = emptypt.Clone()
        empty11 = emptypt.Clone()
        empty12 = emptypt.Clone()
        empty13 = emptypt.Clone()
        empty14 = emptypt.Clone()
        empty15 = emptypt.Clone()
        empty16 = emptypt.Clone()
        empty17 = emptypt.Clone()
        empty18 = emptypt.Clone()
        empty19 = emptypt.Clone()
        empty20 = emptypt.Clone()
        empty21 = emptypt.Clone()
        
        hdyuppt  = systbkgsup.getAddedHist(empty19,"DYJetsToLL","sr","h_h_pt")#.Scale(dynormup)
        hdyuppt.Scale(dynormup)
        hdydwnpt = systbkgsdwn.getAddedHist(empty20,"DYJetsToLL","sr","h_h_pt")#.Scale(dynormdwn)
        hdydwnpt.Scale(dynormdwn)
        httuppt  = systbkgsup.getAddedHist(empty10,"TT","sr","h_h_pt")
        httdwnpt = systbkgsdwn.getAddedHist(empty11,"TT","sr","h_h_pt")
        hzzuppt  = systbkgsup.getAddedHist(empty12,"ZZTo2L2Q","sr","h_h_pt")
        hzzdwnpt = systbkgsdwn.getAddedHist(empty13,"ZZTo2L2Q","sr","h_h_pt")
        hwzuppt  = systbkgsup.getAddedHist(empty14,"WZTo2L2Q","sr","h_h_pt")
        hwzdwnpt = systbkgsdwn.getAddedHist(empty15,"WZTo2L2Q","sr","h_h_pt")
        hvvuppt  = hzzuppt.Clone()
        hvvdwnpt = hzzdwnpt.Clone()
        hvvuppt.Add(hwzuppt)
        hvvdwnpt.Add(hwzdwnpt)
        hdypt    = bkgs.getAddedHist(empty21,"DYJetsToLL","sr","h_h_pt")#.Scale(dynormnom)
        hdypt.Scale(dynormnom)
        httpt    = bkgs.getAddedHist(empty16,"TT","sr","h_h_pt")
        hzzpt    = bkgs.getAddedHist(empty17,"ZZTo2L2Q","sr","h_h_pt")
        hwzpt    = bkgs.getAddedHist(empty18,"WZTo2L2Q","sr","h_h_pt")
        hvvpt    = hzzpt.Clone()
        hvvpt.Add(hwzpt)

    
        #Plotting parameters
        htt.SetLineColor(ROOT.kBlack)
        hvv.SetLineColor(ROOT.kBlack)
        hdy.SetLineColor(ROOT.kBlack)
        httup.SetLineColor(ROOT.kRed)
        hvvup.SetLineColor(ROOT.kRed)
        hdyup.SetLineColor(ROOT.kRed)
        httdwn.SetLineColor(ROOT.kBlue)
        hvvdwn.SetLineColor(ROOT.kBlue)
        hdydwn.SetLineColor(ROOT.kBlue)
        httpt.SetLineColor(ROOT.kBlack)
        hvvpt.SetLineColor(ROOT.kBlack)
        hdypt.SetLineColor(ROOT.kBlack)
        httuppt.SetLineColor(ROOT.kRed)
        hvvuppt.SetLineColor(ROOT.kRed)
        hdyuppt.SetLineColor(ROOT.kRed)
        httdwnpt.SetLineColor(ROOT.kBlue)
        hvvdwnpt.SetLineColor(ROOT.kBlue)
        hdydwnpt.SetLineColor(ROOT.kBlue)
        #httuppt.SetMaximum(10)
        
        dyupx = hdyup.GetXaxis()
        dyupx.SetTitle("Z' Jigsaw Mass Estimator")
        dyuptx = hdyuppt.GetXaxis()
        dyuptx.SetTitle("Higgs Candidate pT")
        ttupx = httup.GetXaxis()
        ttupx.SetTitle("Z' Jigsaw Mass Estimator")
        ttuptx = httuppt.GetXaxis()
        ttuptx.SetTitle("Higgs Candidate pT")
        vvupx = hvvup.GetXaxis()
        vvupx.SetTitle("Z' Jigsaw Mass Estimator")
        vvuptx = hvvuppt.GetXaxis()
        vvuptx.SetTitle("Higgs Candidate pT")
        
        hdyup.SetStats(0)
        httup.SetStats(0)
        hvvup.SetStats(0)
        hdyuppt.SetStats(0)
        httuppt.SetStats(0)
        hvvuppt.SetStats(0)
        
        #Relevant Info
        dystrup  = "dy up, integral: "+str(round(hdyup.Integral(),3))
        dystrnom = "dy nom, integral: "+str(round(hdy.Integral(),3))
        dystrdwn = "dy dwn, integral: "+str(round(hdydwn.Integral(),3))
        ttstrup  = "tt up, integral: "+str(round(httup.Integral(),3))
        ttstrnom = "tt nom, integral: "+str(round(htt.Integral(),3))
        ttstrdwn = "tt dwn, integral: "+str(round(httdwn.Integral(),3))
        vvstrup  = "vv up, integral: "+str(round(hvvup.Integral(),3))
        vvstrnom = "vv nom, integral: "+str(round(hvv.Integral(),3))
        vvstrdwn = "vv dwn, integral: "+str(round(hvvdwn.Integral(),3))
        
        dystruppt  = "dy up, integral: "+str(round(hdyuppt.Integral(),3))
        dystrnompt = "dy nom, integral: "+str(round(hdypt.Integral(),3))
        dystrdwnpt = "dy dwn, integral: "+str(round(hdydwnpt.Integral(),3))
        ttstruppt  = "tt up, integral: "+str(round(httuppt.Integral(),3))
        ttstrnompt = "tt nom, integral: "+str(round(httpt.Integral(),3))
        ttstrdwnpt = "tt dwn, integral: "+str(round(httdwnpt.Integral(),3))
        vvstruppt  = "vv up, integral: "+str(round(hvvuppt.Integral(),3))
        vvstrnompt = "vv nom, integral: "+str(round(hvvpt.Integral(),3))
        vvstrdwnpt = "vv dwn, integral: "+str(round(hvvdwnpt.Integral(),3))
        
        #Labels
        dyjiglab = ROOT.TPaveText(.5,.4,.8,.6,"NBNDC")
        dyjiglab.AddText(dystrup)
        dyjiglab.AddText(dystrnom)
        dyjiglab.AddText(dystrdwn)
        dyjiglab.SetFillColor(0)
        dyjiglab.GetLine(0).SetTextColor(ROOT.kRed)
        dyjiglab.GetLine(1).SetTextColor(ROOT.kBlack)
        dyjiglab.GetLine(2).SetTextColor(ROOT.kBlue)
        
        dyptlab = ROOT.TPaveText(.5,.4,.8,.6,"NBNDC")
        dyptlab.AddText(dystruppt)
        dyptlab.AddText(dystrnompt)
        dyptlab.AddText(dystrdwnpt)
        dyptlab.SetFillColor(0)
        dyptlab.GetLine(0).SetTextColor(ROOT.kRed)
        dyptlab.GetLine(1).SetTextColor(ROOT.kBlack)
        dyptlab.GetLine(2).SetTextColor(ROOT.kBlue)
        
        ttjiglab = ROOT.TPaveText(.5,.4,.8,.6,"NBNDC")
        ttjiglab.AddText(ttstrup)
        ttjiglab.AddText(ttstrnom)
        ttjiglab.AddText(ttstrdwn)
        ttjiglab.SetFillColor(0)
        ttjiglab.GetLine(0).SetTextColor(ROOT.kRed)
        ttjiglab.GetLine(1).SetTextColor(ROOT.kBlack)
        ttjiglab.GetLine(2).SetTextColor(ROOT.kBlue)
        
        ttptlab = ROOT.TPaveText(.5,.4,.8,.6,"NBNDC")
        ttptlab.AddText(ttstruppt)
        ttptlab.AddText(ttstrnompt)
        ttptlab.AddText(ttstrdwnpt)
        ttptlab.SetFillColor(0)
        ttptlab.GetLine(0).SetTextColor(ROOT.kRed)
        ttptlab.GetLine(1).SetTextColor(ROOT.kBlack)
        ttptlab.GetLine(2).SetTextColor(ROOT.kBlue)
        
        vvjiglab = ROOT.TPaveText(.5,.4,.8,.6,"NBNDC")
        vvjiglab.AddText(vvstrup)
        vvjiglab.AddText(vvstrnom)
        vvjiglab.AddText(vvstrdwn)
        vvjiglab.SetFillColor(0)
        vvjiglab.GetLine(0).SetTextColor(ROOT.kRed)
        vvjiglab.GetLine(1).SetTextColor(ROOT.kBlack)
        vvjiglab.GetLine(2).SetTextColor(ROOT.kBlue)
        
        vvptlab = ROOT.TPaveText(.5,.4,.8,.6,"NBNDC")
        vvptlab.AddText(vvstruppt)
        vvptlab.AddText(vvstrnompt)
        vvptlab.AddText(vvstrdwnpt)
        vvptlab.SetFillColor(0)
        vvptlab.GetLine(0).SetTextColor(ROOT.kRed)
        vvptlab.GetLine(1).SetTextColor(ROOT.kBlack)
        vvptlab.GetLine(2).SetTextColor(ROOT.kBlue)
        
        #TCanvases
        tcbkg  = ROOT.TCanvas("tcbkg","backgrounds",1400,800)
        pdy    = ROOT.TPad("pdy","dysr",0,0.5,0.33,1.0)
        ptt    = ROOT.TPad("ptt","ttsr",0.33,0.5,0.66,1.0)
        pvv    = ROOT.TPad("pvv","vvsr",0.66,0.5,1.0,1.0)
        ptdy   = ROOT.TPad("ptdy","dysr",0,0.0,0.33,0.5)
        pttt   = ROOT.TPad("pttt","ttsr",0.33,0.0,0.66,0.5)
        ptvv   = ROOT.TPad("ptvv","vvsr",0.66,0.0,1.0,0.5)
        
        tcbkg.cd()
        pdy.Draw()
        pdy.cd()
        hdyup.Draw('hist')
        hdy.Draw('histsame')
        hdydwn.Draw('histsame')
        dyjiglab.Draw()
        pdy.Update()
        tcbkg.cd()
        ptt.Draw()
        ptt.cd()
        httup.Draw('hist')
        htt.Draw('histsame')
        httdwn.Draw('histsame')
        ttjiglab.Draw()
        ptt.Update()
        tcbkg.cd()
        pvv.Draw()
        pvv.cd()
        hvvup.Draw('hist')
        hvv.Draw('histsame')
        hvvdwn.Draw('histsame')
        vvjiglab.Draw()
        pvv.Update()
        tcbkg.cd()
        ptdy.Draw()
        ptdy.cd()
        hdyuppt.Draw('hist')
        hdypt.Draw('histsame')
        hdydwnpt.Draw('histsame')
        dyptlab.Draw()
        ptdy.Update()
        tcbkg.cd()
        pttt.Draw()
        pttt.cd()
        httuppt.Draw('hist')
        httpt.Draw('histsame')
        httdwnpt.Draw('histsame')
        ttptlab.Draw()
        pttt.Update()
        tcbkg.cd()
        ptvv.Draw()
        ptvv.cd()
        hvvuppt.Draw('hist')
        hvvpt.Draw('histsame')
        hvvdwnpt.Draw('histsame')
        vvptlab.Draw()
        ptvv.Update()
        tcbkg.cd()
        tcbkg.Update()
        
        bkgcomp = go.makeOutFile('Run2_2017_2018_ZllHbbMET',chan+'_bkg_'+syst,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        tcbkg.SaveAs(bkgcomp)
    
        
