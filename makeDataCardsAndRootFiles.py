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
    sig         = go.signal('SignalInputs/',zptcut,hptcut,metcut,btagwp,sigxs,59.74)
    systsigup   = go.signal('SignalJECSyst/',zptcut,hptcut,metcut,btagwp,sigxs,59.74,"systjecup")
    systsigdwn  = go.signal('SignalJECSyst/',zptcut,hptcut,metcut,btagwp,sigxs,59.74,"systjecdwn")
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
    sigcolors = go.colsFromPalette(siginfo,ROOT.kCMYK)
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
                            {"type":"shapeN2","unc":1,"proc":[1,0,1,1]},
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

        #Make Plot for each signal sample
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
        sigplotname = go.makeOutFile(signame,'jecupdowncomp','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        tc.SaveAs(sigplotname)
    #Debug for background

    #Extra plots
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
    #hdyup.SetLineColor(ROOT.kRed)
    httdwn.SetLineColor(ROOT.kBlue)
    hvvdwn.SetLineColor(ROOT.kBlue)
    #hdydwn.SetLineColor(ROOT.kBlue)
    httpt.SetLineColor(ROOT.kBlack)
    hvvpt.SetLineColor(ROOT.kBlack)
    #hdypt.SetLineColor(ROOT.kBlack)
    httuppt.SetLineColor(ROOT.kRed)
    hvvuppt.SetLineColor(ROOT.kRed)
    #hdyuppt.SetLineColor(ROOT.kRed)
    httdwnpt.SetLineColor(ROOT.kBlue)
    hvvdwnpt.SetLineColor(ROOT.kBlue)
    #hdydwnpt.SetLineColor(ROOT.kBlue)
    #httuppt.SetMaximum(10)
    
    ttupx = httup.GetXaxis()
    ttupx.SetTitle("Z' Jigsaw Mass Estimator")
    ttuptx = httuppt.GetXaxis()
    ttuptx.SetTitle("Higgs Candidate pT")
    vvupx = hvvup.GetXaxis()
    vvupx.SetTitle("Z' Jigsaw Mass Estimator")
    vvuptx = hvvuppt.GetXaxis()
    vvuptx.SetTitle("Higgs Candidate pT")

    httup.SetStats(0)
    hvvup.SetStats(0)
    httuppt.SetStats(0)
    hvvuppt.SetStats(0)

    #Relevant Info
    ttstrup  = "tt up, integral: "+str(round(httup.Integral(),3))
    ttstrnom = "tt nom, integral: "+str(round(htt.Integral(),3))
    ttstrdwn = "tt dwn, integral: "+str(round(httdwn.Integral(),3))
    vvstrup  = "vv up, integral: "+str(round(hvvup.Integral(),3))
    vvstrnom = "vv nom, integral: "+str(round(hvv.Integral(),3))
    vvstrdwn = "vv dwn, integral: "+str(round(hvvdwn.Integral(),3))

    ttstruppt  = "tt up, integral: "+str(round(httuppt.Integral(),3))
    ttstrnompt = "tt nom, integral: "+str(round(httpt.Integral(),3))
    ttstrdwnpt = "tt dwn, integral: "+str(round(httdwnpt.Integral(),3))
    vvstruppt  = "vv up, integral: "+str(round(hvvuppt.Integral(),3))
    vvstrnompt = "vv nom, integral: "+str(round(hvvpt.Integral(),3))
    vvstrdwnpt = "vv dwn, integral: "+str(round(hvvdwnpt.Integral(),3))

    #Labels
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
    #hdy.Draw()
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
    #hdypt.Draw()
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

    bkgjeccomp = go.makeOutFile('Run2_2017_2018_ZllHbbMET',chan+'_bkg_jeccomp','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tcbkg.SaveAs(bkgjeccomp)
    
        
