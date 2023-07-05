import tdrstyle
import CMS_lumi
import ROOT
import glob
import os
import gecorg_test as go
#import gecorg as gog
import numpy as np
import pandas as pd
import configparser
import argparse


tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = "137.6 fb^{-1}"
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

def getErrorOnIntegral(hist):
    sumsqrs = 0
    for i in range(hist.GetNbinsX()+1):
        err = hist.GetBinError(i)
        sumsqrs += err*err
    errint = sumsqrs**(1/2)
    return errint

def getErrorOnDivison(vallist,errlist):
    valerrs = zip(vallist,errlist)
    relunsqr = [(pair[1]/pair[0])*(pair[1]/pair[0]) for pair in valerrs]#relative error squared
    sums = sum(relunsqr)
    sqrtsq = sums**(1/2)
    uncondiv = vallist[1]/vallist[0]*sqrtsq
    return uncondiv
    

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-dir","--directory", type=str,help = "date folder with output")
    parser.add_argument("-v","--validationregion", type=bool,help = "is this a validation region?")
    args = parser.parse_args()


    pathbkg    = args.directory#'pfMETNominal/'
    pathdata   = args.directory#'pfMETNominal/'
    zptcut  = args.zptcut#'150.0'
    hptcut  = args.hptcut#'300.0'
    metcut  = args.metcut#'200.0'
    btagwp  = args.btagwp#'0.8'

    systr = 'systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco'

    ttmcemu  = go.backgrounds('emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref',zptcut,hptcut,metcut,btagwp,systr)
    ttmcmumu = go.backgrounds('mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref',zptcut,hptcut,metcut,btagwp,systr)
    #ttmcemu  = go.backgrounds('analysis_output_ZpAnomalon/2023-03-01',zptcut,hptcut,metcut,btagwp,systr)#prefire
    #ttmcmumu = go.backgrounds('analysis_output_ZpAnomalon/2023-03-01/muonchannel',zptcut,hptcut,metcut,btagwp,systr)#prefire
    #dataemu  = go.run2('emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref',zptcut,hptcut,metcut,btagwp,'systnominal_btagnom_muidnom')

    print(ttmcmumu.bkgs['TT'][18]['tr'])
    print(ttmcemu.bkgs['TT'][18]['tr'])
    #print(dataemu.data)
    
    tf1 = ROOT.TFile(ttmcemu.a18tttr[0])
    emptypt = tf1.Get('h_z_pt').Clone()
    emptyphi = tf1.Get('h_z_phiw').Clone()
    emptyeta = tf1.Get('h_z_eta').Clone()
    emptypt.Reset("ICESM")
    emptyphi.Reset("ICESM")
    emptyeta.Reset("ICESM")
    empty1 = emptypt.Clone()
    empty2 = emptyphi.Clone()
    empty3 = emptyeta.Clone()
    empty4 = emptypt.Clone()
    empty5 = emptyphi.Clone()
    empty6 = emptyeta.Clone()
    empty7 = emptypt.Clone()
    empty8 = emptypt.Clone()
    empty9 = emptypt.Clone()
    empty10 = emptypt.Clone()

    print("Getting the emu channel stuffs")
    httemupt = ttmcemu.getAddedHist(empty1,"TT","tr","h_z_pt")
    httemuphi = ttmcemu.getAddedHist(empty2,"TT","tr","h_z_phiw")
    httemueta = ttmcemu.getAddedHist(empty3,"TT","tr","h_z_eta")
    httemupt.Rebin(2)
    httemuphi.Rebin(2)
    httemueta.Rebin(2)

    print("################Getting the mumu channel stuffs")
    httmumupt = ttmcmumu.getAddedHist(empty4,"TT","tr","h_z_pt")
    httmumuphi = ttmcmumu.getAddedHist(empty5,"TT","tr","h_z_phiw")
    httmumueta = ttmcmumu.getAddedHist(empty6,"TT","tr","h_z_eta")
    httmumupt.Rebin(2)
    httmumuphi.Rebin(2)
    httmumueta.Rebin(2)

    httemupt.SetLineColor(ROOT.kRed)
    httemuphi.SetLineColor(ROOT.kRed)
    httemueta.SetLineColor(ROOT.kRed)

    httmumupt.SetLineColor(ROOT.kBlue)
    httmumuphi.SetLineColor(ROOT.kBlue)
    httmumueta.SetLineColor(ROOT.kBlue)
    httmumuphi.SetMaximum(15)
    httmumuphi.SetMinimum(0)
    httmumueta.SetMaximum(20)

    httmumupt.GetXaxis().SetTitle("Z candidate p_{T}")
    httmumuphi.GetXaxis().SetTitle("Z candidate $\phi$")
    httmumueta.GetXaxis().SetTitle("Z candidate $\eta$")
    httmumupt.GetYaxis().SetTitle("Events / {0} GeV".format(httmumupt.GetBinWidth(1)))
    httmumuphi.GetYaxis().SetTitle("Events")
    httmumueta.GetYaxis().SetTitle("Events")

    httptdiv = httmumupt.Clone()
    httphidiv = httmumuphi.Clone()
    httetadiv = httmumueta.Clone()

    httptdiv.Divide(httmumupt,httemupt)
    httphidiv.Divide(httmumuphi,httemuphi)
    httetadiv.Divide(httmumueta,httemueta)

    httptdiv.SetMarkerStyle(8)
    httptdiv.SetMarkerSize(0.5)
    httptdiv.SetMarkerColor(ROOT.kBlack)
    httptdiv.GetYaxis().SetRangeUser(-.5,2.5)
    httptdiv.GetYaxis().SetTitle("mumu/emu")
    httptdiv.GetYaxis().SetTitleSize(0.15)
    httptdiv.GetYaxis().SetTitleOffset(0.3)
    httptdiv.GetYaxis().SetLabelSize(0.12)
    httptdiv.GetYaxis().SetLabelOffset(0.017)
    httptdiv.GetYaxis().SetNdivisions(503)
    httptdiv.GetXaxis().SetLabelSize(0.10)
    httptdiv.GetXaxis().SetLabelOffset(0.017)

    httphidiv.SetMarkerStyle(8)
    httphidiv.SetMarkerSize(0.5)
    httphidiv.SetMarkerColor(ROOT.kBlack)
    httphidiv.GetYaxis().SetRangeUser(-.5,2.5)
    httphidiv.GetYaxis().SetTitle("mumu/emu")
    httphidiv.GetYaxis().SetTitleSize(0.15)
    httphidiv.GetYaxis().SetTitleOffset(0.3)
    httphidiv.GetYaxis().SetLabelSize(0.12)
    httphidiv.GetYaxis().SetLabelOffset(0.017)
    httphidiv.GetYaxis().SetNdivisions(503)
    httphidiv.GetXaxis().SetLabelSize(0.10)
    httphidiv.GetXaxis().SetLabelOffset(0.017)

    httetadiv.SetMarkerStyle(8)
    httetadiv.SetMarkerSize(0.5)
    httetadiv.SetMarkerColor(ROOT.kBlack)
    httetadiv.GetYaxis().SetRangeUser(-.5,2.5)
    httetadiv.GetYaxis().SetTitle("mumu/emu")
    httetadiv.GetYaxis().SetTitleSize(0.15)
    httetadiv.GetYaxis().SetTitleOffset(0.3)
    httetadiv.GetYaxis().SetLabelSize(0.12)
    httetadiv.GetYaxis().SetLabelOffset(0.017)
    httetadiv.GetYaxis().SetNdivisions(503)
    httetadiv.GetXaxis().SetLabelSize(0.10)
    httetadiv.GetXaxis().SetLabelOffset(0.017)

    intttmumupt = httmumupt.Integral()
    intttemupt = httemupt.Integral()
    intttmumuphi = httmumuphi.Integral()
    intttemuphi = httemuphi.Integral()
    intttmumueta = httmumueta.Integral()
    intttemueta = httemueta.Integral()

    print("mumu pt: ",intttmumupt)
    print("emu pt: ",intttemupt)
    print("mumu phi: ",intttmumuphi)
    print("emu phi: ",intttemuphi)
    print("mumu eta: ",intttmumueta)
    print("emu eta: ",intttemueta)
    
    print("scale from ZpT cut: ",intttmumuphi/intttemuphi)

    erremu  = getErrorOnIntegral(httemuphi)
    errmumu = getErrorOnIntegral(httmumuphi)

    yields = [intttmumuphi,intttemuphi]
    uncs   = [errmumu,erremu]

    unconscale = getErrorOnDivison(yields,uncs)
    scaleFromMC = intttmumuphi/intttemuphi

    print("scale from ZpT cut: {0} +- {1}".format(scaleFromMC,unconscale))

    httmumupt.GetYaxis().SetRangeUser(0,70)
    httmumuphi.GetYaxis().SetRangeUser(0,45)
    httmumueta.GetYaxis().SetRangeUser(0,150)

    #Save the normalization
    scalefilename = go.makeOutFile('Run2_2016_2017_2018','ttemunormalization_'+systr,'.npy',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    scalefile     = open(scalefilename,'wb')
    np.save(scalefile,np.array([scaleFromMC,unconscale]))
    scalefile.close()

    #scale derivation plots
    tc = ROOT.TCanvas("tc","tt channel overlays",1500,600)
    histpad1 = ROOT.TPad("histpad1","emuag",0,.2,.33,1)
    ratpad1 = ROOT.TPad("ratpad1","ratio",0,0,.33,.2)
    histpad2 = ROOT.TPad("histpad2","emuagsr",.33,.2,.67,1)
    ratpad2  = ROOT.TPad("ratpad2","ratiosr",.33,0,.67,.2)
    histpad3 = ROOT.TPad("histpad3","emuagsr",.67,.2,1,1)
    ratpad3  = ROOT.TPad("ratpad3","ratiosr",.67,0,1,.2)
    ratlinept = ROOT.TLine(httptdiv.GetBinLowEdge(1),1,httptdiv.GetBinWidth(1)*httptdiv.GetNbinsX(),1)
    ratlinephi = ROOT.TLine(httphidiv.GetBinLowEdge(1),1,httphidiv.GetBinWidth(1)*httphidiv.GetNbinsX(),1)
    ratlineeta = ROOT.TLine(httetadiv.GetBinLowEdge(1),1,httetadiv.GetBinWidth(1)*httetadiv.GetNbinsX(),1)
    ttleg  = ROOT.TLegend(0.55,0.65,0.9,0.8)
    ttleg.AddEntry(httemupt,"TT , $e\mu$","l")
    ttleg.AddEntry(httmumupt,"TT , $\mu\mu$","l")
    ttleg.SetBorderSize(0)
    normcomp = ROOT.TPaveText(0.45,0.3,0.90,0.55,"NBNDC")
    normcomp.AddText("mumu/emu = {0} +- {1}".format(round((intttmumuphi/intttemuphi),3),round(unconscale,3)))
    normcomp.SetFillColor(0)

    tc.cd()
    histpad1.Draw()
    histpad1.cd()
    httmumupt.Draw('hist')
    CMS_lumi.CMS_lumi(histpad1,4,13)
    httemupt.Draw('histsame')
    ttleg.Draw()
    normcomp.Draw()
    histpad1.Update()
    tc.cd()
    ratpad1.Draw()
    ratpad1.cd()
    httptdiv.Draw()
    ratlinept.Draw()
    tc.cd()
    histpad2.Draw()
    histpad2.cd()
    httmumuphi.Draw('hist')
    CMS_lumi.CMS_lumi(histpad2,4,13)
    httemuphi.Draw('histsame')
    ttleg.Draw()
    histpad2.Update()
    tc.cd()
    ratpad2.Draw()
    ratpad2.cd()
    httphidiv.Draw()
    ratlinephi.Draw()
    tc.cd()
    histpad3.Draw()
    histpad3.cd()
    httmumueta.Draw('hist')
    CMS_lumi.CMS_lumi(histpad3,4,13)
    httemueta.Draw('histsame')
    ttleg.Draw()
    histpad3.Update()
    tc.cd()
    ratpad3.Draw()
    ratpad3.cd()
    httetadiv.Draw()
    ratlineeta.Draw()

    ttover = go.makeOutFile('Run2_161718','ttbar_emu_mumu_mc_comp'+systr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(ttover)

