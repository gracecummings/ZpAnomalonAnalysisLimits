import tdrstyle
import CMS_lumi
import ROOT
import glob
import os
import gecorg_test as go
#import gecorg as gog
import numpy as np
#import pandas as pd
#import configparser
import argparse

tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = "59.74 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"

def plotMzp(pad,hist,islog=False,logmin=0.1,isData=False):
    maxi = hist.GetMaximum()
    #print("Plotting maximum: ",maxi)
    mr   = round(maxi,0)
    histmax = mr+mr*0.30
    histmin = 0
    if islog:
        histmax = mr*10
        histmin = logmin

    #print("Maximum used to draw (i hope): ",histmax)
    #hist.SetMaximum(histmax)
    #hist.SetMinimum(histmin)
    hist.SetMarkerStyle(8)
    hist.SetMarkerSize(0.5)
    if isData:
        hist.SetMarkerColor(ROOT.kBlack)
        hist.SetLineColor(ROOT.kBlack)
        drawopts = "SAMEE2"
    else:
        hist.SetMarkerColor(ROOT.kBlue)
        drawopts = "E1"
    xax = hist.GetXaxis()
    yax = hist.GetYaxis()
    xax.SetTitle("M_{Z'}")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 200 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    yax.SetLimits(histmin,histmax)
    
    hist.Draw(drawopts)

def setPlotParamsMC(plot,title,color):
    plot.SetFillColor(color)
    plot.SetLineColor(color)
    plot.GetXaxis().SetTitle(title)
    plot.GetYaxis().SetTitle("Events / {0}".format(plot.GetBinWidth(1)))

def setPlotParamsData(plot,title):
    plot.SetMarkerStyle(8)
    plot.SetMarkerSize(0.5)
    plot.SetMarkerColor(ROOT.kBlack)
    plot.SetLineColor(ROOT.kBlack)
    plot.GetXaxis().SetTitle(title)
    plot.GetYaxis().SetTitle("Events / {0}".format(plot.GetBinWidth(1)))

def makeRatioHist(hnumerator,hdenominator,title):
    hist = hnumerator.Clone()
    hist.Divide(hnumerator,hdenominator)
    hist.SetMarkerSize(0.5)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.GetYaxis().SetRangeUser(-0.5,2.5)
    hist.GetYaxis().SetTitle(title)
    hist.GetYaxis().SetTitleSize(0.15)
    hist.GetYaxis().SetTitleOffset(0.3)
    hist.GetYaxis().SetLabelSize(0.12)
    hist.GetYaxis().SetLabelOffset(0.017)
    hist.GetYaxis().SetNdivisions(503)
    hist.GetXaxis().SetLabelSize(0.10)
    hist.GetXaxis().SetLabelOffset(0.017)
    return hist

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-dir","--directory", type=str,help = "date folder with output")
    args = parser.parse_args()
    
    #pathbkg    = args.directory#This should be to mumu MC
    #pathdata   = args.directory#This should be to emu

    pathbkg  = 'pfMET_nominal'#This should be to mumu MC
    pathdata = 'pfMET_emu'#This should be to emu
    zptcut  = args.zptcut#'150.0'
    hptcut  = args.hptcut#'300.0'
    metcut  = args.metcut#'200.0'
    btagwp  = args.btagwp#'0.8'

    systr = 'systnominal_btagnom'

    emuscale = np.load(pathdata+'/Run2_2017_2018_ttemunormalization_systnominal_btagnom_Zptcut100.0_Hptcut0.0_metcut0.0_btagwp0.0.npy')[0]

    print("Scaling the emu yield by ",emuscale)

    bkgs  = go.backgrounds(pathbkg,zptcut,hptcut,metcut,btagwp,systr)
    data  = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,systr)

    #print(bkgs.bkgs)
    #print(data.data)

    tf1 = ROOT.TFile(bkgs.a18ttsb[0])
    emptymzp = tf1.Get('h_zp_jigm').Clone()
    emptymzp.Reset("ICESM")
    emptysdm = tf1.Get('h_h_sd').Clone()
    emptysdm.Reset("ICESM")

    emptymzp1 = emptymzp.Clone()
    emptymzp2 = emptymzp.Clone()
    emptymzp3 = emptymzp.Clone()
    emptymzp4 = emptymzp.Clone()

    emptysdm1 = emptysdm.Clone()
    emptysdm2 = emptysdm.Clone()

    #SB would be the input to the alpha method. SR would be the SR estimation
    hsbtt = bkgs.getAddedHist(emptymzp1,"TT","sb","h_zp_jigm",[18])
    hsbdat = data.getAddedHist(emptymzp2,"sb","h_zp_jigm",[18])
    hsrtt = bkgs.getAddedHist(emptymzp3,"TT","sr","h_zp_jigm",[18])
    hsrdat = data.getAddedHist(emptymzp4,"sr","h_zp_jigm",[18])#This is the emu channel,so we are good
    hsbtt.Rebin(2)
    hsbdat.Rebin(2)
    hsrtt.Rebin(2)
    hsrdat.Rebin(2)
    
    #TR is the dataset for the background template for the alpha method norm
    htrtt = bkgs.getAddedHist(emptysdm1,"TT","tr","h_h_sd",[18])
    htrdat = data.getAddedHist(emptysdm2,"tr","h_h_sd",[18])

    #Scale'em
    hsbdat.Scale(emuscale)
    hsrdat.Scale(emuscale)
    htrdat.Scale(emuscale)

    #make'em purty
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    setPlotParamsMC(hsbtt,"Z' Jigsaw Mass Estimator M_{Zp}",bkgcols[1])
    setPlotParamsMC(hsrtt,"Z' Jigsaw Mass Estimator M_{Zp}",bkgcols[1])
    setPlotParamsMC(htrtt,"Higgs Candidate soft drop mass",bkgcols[1])
    setPlotParamsData(hsrdat,"Z' Jigsaw Mass Estimator M_{Zp}")
    setPlotParamsData(hsbdat,"Z' Jigsaw Mass Estimator M_{Zp}")
    setPlotParamsData(htrdat,"Higgs Candidate soft drop mass")

    #make ratios
    hdivttsb = makeRatioHist(hsbdat,hsbtt,"emudata/mc")
    hdivttsr = makeRatioHist(hsrdat,hsrtt,"emudata/mc")
    hdivtttr = makeRatioHist(htrdat,htrtt,"emudata/mc")
    
    tc = ROOT.TCanvas("tc","tt channel overlays",1500,600)
    histpad1 = ROOT.TPad("histpad1","pad1",0,.2,.33,1)
    ratpad1 = ROOT.TPad("ratpad1","ratio1",0,0,.33,.2)
    histpad2 = ROOT.TPad("histpad2","pad2",.33,.2,.67,1)
    ratpad2  = ROOT.TPad("ratpad2","ratio2",.33,0,.67,.2)
    histpad3 = ROOT.TPad("histpad3","pad3",.67,.2,1,1)
    ratpad3  = ROOT.TPad("ratpad3","ratio3",.67,0,1,.2)
    ratlinemzp = ROOT.TLine(hdivttsb.GetBinLowEdge(1),1,hdivttsb.GetBinWidth(1)*hdivttsb.GetNbinsX(),1)
    ratlinemsd = ROOT.TLine(hdivtttr.GetBinLowEdge(1),1,hdivtttr.GetBinWidth(1)*hdivtttr.GetNbinsX(),1)
    ttleg  = ROOT.TLegend(0.55,0.60,0.9,0.8)
    ttleg.AddEntry(hsbtt,"TT, $\mu\mu$","f" )
    ttleg.AddEntry(hsbdat,"data, $e\mu$ scaled","ep" )
    ttleg.SetBorderSize(0)

    descriptxt0 = ROOT.TPaveText(0.60,0.40,0.94,0.55,"NBNDC")
    descriptxt0.AddText("Preselection")
    descriptxt0.AddText("Soft Drop Mass Total Region")
    descriptxt0.SetFillColor(0)


    descriptxt1 = ROOT.TPaveText(0.55,0.3,0.90,0.45,"NBNDC")
    descriptxt1.AddText("Preselection")
    descriptxt1.AddText("Soft Drop Mass Side Band")
    descriptxt1.SetFillColor(0)

    descriptxt2 = ROOT.TPaveText(0.5,0.3,0.90,0.45,"NBNDC")
    descriptxt2.AddText("Preselection")
    descriptxt2.AddText("Soft Drop Mass Signal Region")
    descriptxt2.SetFillColor(0)

    tc.cd()
    histpad1.Draw()
    histpad1.cd()
    htrtt.Draw("HIST")
    CMS_lumi.CMS_lumi(histpad1,4,13)
    htrdat.Draw("SAME")
    ttleg.Draw()
    descriptxt0.Draw()
    histpad1.Update()
    tc.cd()
    ratpad1.Draw()
    ratpad1.cd()
    hdivtttr.Draw()
    ratlinemsd.Draw()
    tc.cd()
    histpad2.Draw()
    histpad2.cd()
    hsbtt.Draw("HIST")
    CMS_lumi.CMS_lumi(histpad2,4,13)
    hsbdat.Draw("SAME")
    ttleg.Draw()
    descriptxt1.Draw()
    histpad2.Update()
    tc.cd()
    ratpad2.Draw()
    ratpad2.cd()
    hdivttsb.Draw()
    ratlinemzp.Draw()
    tc.cd()
    histpad3.Draw()
    histpad3.cd()
    hsrtt.Draw("HIST")
    CMS_lumi.CMS_lumi(histpad3,4,13)
    hsrdat.Draw("SAME")
    ttleg.Draw()
    descriptxt2.Draw()
    histpad3.Update()
    tc.cd()
    ratpad3.Draw()
    ratpad3.cd()
    hdivttsr.Draw()
    ratlinemzp.Draw()
    
    extrpplots = go.makeOutFile('Run2_2018','emu_scaled_set_on_mumu_mc_'+systr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(extrpplots)
