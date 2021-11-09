import argparse
import tdrstyle
import CMS_lumi
import ROOT
import glob
import os
import gecorg_test as go
import numpy as np
import pandas as pd
import configparser

tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = "101.27 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"

def setLogAxis(pad,islog):
    if islog:
        pad.SetLogy()

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

def getMaxXY(hist):
    binmax = hist.GetMaximumBin()
    xmax   = hist.GetBinCenter(binmax)
    hmax   = hist.GetMaximum()
    return(int(binmax),int(xmax),float(hmax))

def getMinXY(hist,minthresh):
    laststatsbin = hist.FindLastBinAbove(minthresh)
    lastbincent = hist.GetBinCenter(hist.GetNbinsX())
    if laststatsbin > 0:
        lastbincent = hist.GetBinCenter(laststatsbin)
    else:
        laststatsbin = hist.GetNbinsX()
    return(int(laststatsbin),int(lastbincent))

    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    #parser.add_argument("-dir","--directory", type=str,help = "date folder with output")
    #probably should add systematic stuff as a input
    parser.add_argument("-v","--validationregion", type=bool,help = "is this a validation region?")
    args = parser.parse_args()

    #Get systematic info
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()
    
    #Starting parameters
    zptcut  = str(args.zptcut)#'150.0'
    hptcut  = str(args.hptcut)#'300.0'
    metcut  = str(args.metcut)#'200.0'
    btagwp  = str(args.btagwp)#'0.8'
    dynorm  = 1.
    validation = False
    rstr = "signalblind"
    rebindiv = 2
    systname = 'nominal'
    systclass = 'pathnom'
    syststr = 'strnom'
    sigmabars = 1#1 sigma bands

    #Get the samples
    pathbkg    = config.get(systname,systclass)
    pathdata   = 'pfMETNominal_CorrCounting_Optimization/'
    datastr = 'systnominal_btagnom'
    
    if validation:
        rstr = "validationblind"
        dynorm = np.load(path+'/Run2_2017_2018_dynormalization_validationblind_Zptcut150.0_Hptcut300.0_metcut200.0_btagwp0.8.npy')[0]
        bkgs = go.validation('analysis_output_ZpAnomalon/2021-07-14',zptcut,hptcut,metcut,btagwp)
        data = go.validation('analysis_output_ZpAnomalon/2021-07-14',zptcut,hptcut,metcut,btagwp)
        sbstring = "30 < m_{hcand,SD} <= 55"
        srstring = "55 < m_{hcand,SD} < 70"
        srregionstring = "Validation Region"
    else:
        dynorm = np.load(pathbkg+'Run2_2017_2018_dynormalization_'+config.get(systname,syststr)+'_signalblind_Zptcut'+zptcut+'_Hptcut'+hptcut+'_metcut'+metcut+'_btagwp'+btagwp+'.npy')[0]
        bkgs = go.backgrounds(pathbkg,zptcut,hptcut,metcut,btagwp,config.get(systname,syststr))
        data = go.run2(pathdata,zptcut,hptcut,metcut,btagwp)
        sbstring = "30 < m_{hcand,SD} < 70"
        srstring = "110 <= m_{hcand,SD} < 150"
        srregionstring = "Signal Region"
        
    print("Using the DY normalization factor: ",dynorm)
    
    #Get the histograms needed
    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    empty4 = empty.Clone()
    empty5 = empty.Clone()
    empty6 = empty.Clone()
    empty7 = empty.Clone()
    empty8 = empty.Clone()
    empty9 = empty.Clone()
    empty10 = empty.Clone()

    hsbdy = bkgs.getAddedHist(empty,"DYJetsToLL","sb","h_zp_jigm")
    hsrdy = bkgs.getAddedHist(empty2,"DYJetsToLL","sr","h_zp_jigm")
    hsbtt = bkgs.getAddedHist(empty3,"TT","sb","h_zp_jigm")
    hsbzz = bkgs.getAddedHist(empty4,"ZZTo2L2Q","sb","h_zp_jigm")
    hsbwz = bkgs.getAddedHist(empty5,"WZTo2L2Q","sb","h_zp_jigm")
    hsbvv = hsbzz.Clone()
    hsbvv.Add(hsbwz)

    if not validation:
        hdatsb = data.getAddedHist(empty9,"sb","h_zp_jigm")
    if validation:
        hdatsb = data.getAddedHistData(empty9,"sb","h_zp_jigm")
        hdatvr = data.getAddedHistData(empty10,"vr","h_zp_jigm")
    hdatsbsub = hdatsb.Clone()

    #rebin
    hsbdy.Rebin(rebindiv)
    hsrdy.Rebin(rebindiv)
    hsbtt.Rebin(rebindiv)
    hsbvv.Rebin(rebindiv)
    hdatsbsub.Rebin(rebindiv)
    hdatsb.Rebin(rebindiv)

    #Apply the normalization
    hsbdy.Scale(dynorm)
    hsrdy.Scale(dynorm)
    

    vvsbmaxbin,vvsbxmax,vvsbymax = getMaxXY(hsbvv)
    ttsbmaxbin,ttsbxmax,ttsbymax = getMaxXY(hsbtt)
    dysbmaxbin,dysbxmax,dysbymax = getMaxXY(hsbdy)
    dysrmaxbin,dysrxmax,dysrymax = getMaxXY(hsrdy)
    datsbmaxbin,datsbxmax,datsbymax = getMaxXY(hdatsb)

    maxvals = [datsbxmax,dysrxmax,dysbxmax,ttsbxmax,vvsbxmax]

    vvsbminbin,vvsbxmin = getMinXY(hsbvv,0.1)
    ttsbminbin,ttsbxmin = getMinXY(hsbtt,0.1)
    dysbminbin,dysbxmin = getMinXY(hsbdy,0.1)
    dysrminbin,dysrxmin = getMinXY(hsrdy,0.1)
    datsbminbin,datsbxmin = getMinXY(hdatsb,0.1)

    #print(" x of the vv max: ",vvsbxmax)
    #print(" x of the tt max: ",ttsbxmax)
    #print(" x of the dy max: ",dysbxmax)
    #print(" x of the dy max: ",dysrxmax)
    #print(" x of dat sb max: ",datsbxmax)
    #print("\n")
    #print("\n")
    #print("x of the vv sb threshold min: ",vvsbxmin)
    #print("x of the tt sb threshold min: ",ttsbxmin)
    #print("x of the dy sb threshold min: ",dysbxmin)
    #print("x of the dy sr threshold min: ",dysrxmin)
    #print("x of data sb threshold min:   ",datsbxmin)
    #print("\n")
    #print("\n")
    

    ROOT.gSystem.CompileMacro("../ZpAnomalonAnalysisUproot/cfunctions/alphafits.C","kfc")
    ROOT.gSystem.Load("../ZpAnomalonAnalysisUproot/cfunctions/alphafits_C")

    #Do fits, make shit
    print("================= doing dy sb fit =================")
    sbfit = ROOT.expFit(hsbdy,"sbl","QR0+",dysbxmax,dysbxmin)#5000)
    uncbands = ROOT.expFitErrBands(hsbdy,"Qsbl","R0+",sigmabars,dysbxmax,dysbxmin)#5000)
    print("================= doing dy sr fit ==================")
    srdyunc = ROOT.expFitErrBands(hsrdy,"sbl","QR0+",sigmabars,dysrxmax,dysrxmin)#4000)
    srfit = ROOT.expFit(hsrdy,"srl","QR0+",dysrxmax,dysrxmin)#4000)#be aware, diff range from err and sb
    print("================== doing tt sb fit ==================")
    sbttfit = ROOT.expFit(hsbtt,"sbl","QR0+",ttsbxmax,ttsbxmin)#3000)
    sbttunc = ROOT.expFitErrBands(hsbtt,"sbl","QR0+",sigmabars,ttsbxmax,ttsbxmin)#3000)
    print("================== doing VV sb fit ==================")
    sbvvfit = ROOT.expFit(hsbvv,"sbl","QR0+",vvsbxmax,vvsbxmin)#3000)
    sbvvunc = ROOT.expFitErrBands(hsbvv,"sbl","QR0+",sigmabars,vvsbxmax,vvsbxmin)#3000)
    print("================= doing data sb fit =================")
    sbdatfit = ROOT.expFit(hdatsb,"sbl","QR0+",datsbxmax,datsbxmin)#5000)
    sbdatunc = ROOT.expFitErrBands(hdatsb,"sbl","QR0+",sigmabars,datsbxmax,datsbxmin)#5000)
    print("================= doing alpha ratio fit =============")
    alpha = ROOT.alphaRatioMakerExp(hsbdy,hsrdy)


    #Set some paremeter for plotting
    uncbands.SetStats(ROOT.kFALSE)
    uncbands.SetFillColorAlpha(2,0.35)
    uncbands.SetMarkerStyle(8)
    uncbands.SetMarkerSize(0)
    sbfit.SetMarkerStyle(8)
    sbfit.SetMarkerSize(0)
    srdyunc.SetFillColor(2)
    srdyunc.SetMarkerSize(0)
    sbttunc.SetFillColor(2)
    sbttunc.SetMarkerSize(0)
    sbvvunc.SetFillColor(2)
    sbvvunc.SetMarkerSize(0)
    sbdatunc.SetMarkerStyle(8)
    sbdatunc.SetMarkerSize(.5)
    sbdatunc.SetMarkerColor(ROOT.kBlack)
    sbdatunc.SetLineColor(ROOT.kBlack)
    hdatsb.SetLineColor(ROOT.kBlack)
    #hsbdy.GetXaxis().SetRangeUser(1500,5000)
    #hsrdy.GetXaxis().SetRangeUser(1500,5000)
    #hsbtt.GetXaxis().SetRangeUser(1500,5000)
    #hsrtt.GetXaxis().SetRangeUser(1500,5000)
    #hsrvv.GetXaxis().SetRangeUser(1500,5000)
    #hsbvv.GetXaxis().SetRangeUser(1500,5000)
    #hdatsbsub.GetXaxis().SetRangeUser(1500,5000)
    #hdatsb.GetXaxis().SetRangeUser(1500,5000)

    #Let's get colorful
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)

    #Bkg hists stacks
    hsbkg = ROOT.THStack("hsbkg","")
    hsbdyc = hsbdy.Clone()
    hsbdyc.SetFillColor(bkgcols[0])
    hsbdyc.SetLineColor(bkgcols[0])
    hsbttc = hsbtt.Clone()
    hsbttc.SetFillColor(bkgcols[1])
    hsbttc.SetLineColor(bkgcols[1])
    hsbvvc = hsbvv.Clone() 
    hsbvvc.SetFillColor(bkgcols[2])
    hsbvvc.SetLineColor(bkgcols[2])
    hsbkg.Add(hsbvvc)
    hsbkg.Add(hsbttc)
    hsbkg.Add(hsbdyc)

    #Bkg Fit stacks
    hsbkgfit = ROOT.THStack("hsbkgfits","")
    sbdyfitc = uncbands.Clone()
    sbttfitc = sbttunc.Clone()
    sbvvfitc = sbvvunc.Clone()
    sbdyfitc.SetFillColor(bkgcols[0])
    sbdyfitc.SetLineColor(bkgcols[0])
    sbttfitc.SetFillColor(bkgcols[1])
    sbttfitc.SetLineColor(bkgcols[1])
    sbvvfitc.SetFillColor(bkgcols[2])
    sbvvfitc.SetLineColor(bkgcols[2])
    hsbkgfit.Add(sbvvfitc)
    hsbkgfit.Add(sbttfitc)
    hsbkgfit.Add(sbdyfitc)

    #Build an added MC hist histogram for errors
    hbkg = hsbdy.Clone()
    hbkg.Add(hsbtt)
    hbkg.Add(hsbvv)


    #Build an added MC fit hist
    hfits = sbdyfitc.Clone()
    hfits.Add(sbttfitc)
    hfits.Add(sbvvfitc)
    
    #Build the data/stackbkg ratio
    hdivnom = hdatsb.Clone()
    hdivnom.Divide(hdatsb,hbkg)
    hdivnom.SetMarkerStyle(8)
    hdivnom.SetMarkerSize(.5)
    hdivnom.GetYaxis().SetRangeUser(0,2)
    hdivnom.GetYaxis().SetTitle("data/MC")
    hdivnom.GetYaxis().SetTitleSize(0.15)
    hdivnom.GetYaxis().SetTitleOffset(0.3)
    hdivnom.GetYaxis().SetLabelSize(0.12)
    hdivnom.GetYaxis().SetLabelOffset(0.017)
    hdivnom.GetYaxis().SetNdivisions(503)
    hdivnom.GetXaxis().SetLabelSize(0.10)
    hdivnom.GetXaxis().SetLabelOffset(0.017)

    hdivfit = hdatsb.Clone()
    hdivfit.Divide(hfits)
    hdivfit.SetMarkerStyle(8)
    hdivfit.SetMarkerSize(.5)
    hdivfit.GetYaxis().SetRangeUser(0,2)
    hdivfit.GetYaxis().SetTitle("data/MC")
    hdivfit.GetYaxis().SetTitleSize(0.15)
    hdivfit.GetYaxis().SetTitleOffset(0.3)
    hdivfit.GetYaxis().SetLabelSize(0.12)
    hdivfit.GetYaxis().SetLabelOffset(0.017)
    hdivfit.GetYaxis().SetNdivisions(503)
    hdivfit.GetXaxis().SetLabelSize(0.10)
    hdivfit.GetXaxis().SetLabelOffset(0.017)

    #Make the first TCanvas. This is the shape fits
    tc = ROOT.TCanvas("tc","shapes",1100,800)
    p11 = ROOT.TPad("p11","dysr",0,0,0.33,.5)
    p12 = ROOT.TPad("p12","ttsr",0.33,0,0.66,.5)
    p13 = ROOT.TPad("p13","vvsr",0.66,0,1.0,.5)
    p21 = ROOT.TPad("p21","dysb",0,.5,0.33,1.0)
    p22 = ROOT.TPad("p22","ttsb",0.33,.5,0.66,1.0)
    p23 = ROOT.TPad("p23","vvsb",0.66,.5,1.0,1.0)

    #Define the legends
    l21 = ROOT.TLegend(0.55,0.65,0.9,0.8)
    l11 = ROOT.TLegend(0.55,0.65,0.9,0.8)
    l22 = ROOT.TLegend(0.55,0.65,0.9,0.8)
    l23 = ROOT.TLegend(0.55,0.65,0.9,0.8)
    l111 = ROOT.TLegend(0.55,0.4,0.9,0.8)

    #Define the labels
    label = ROOT.TPaveText(.5,.4,.9,.5,"NBNDC")
    label.AddText(sbstring)
    label.AddText("150 < m_{hcand,SD}")
    label.SetFillColor(0)
    label2 = ROOT.TPaveText(.5,.4,.9,.5,"NBNDC")
    label2.AddText(srstring)#higgs mass
    label2.AddText(srregionstring)#higgs mass
    label2.SetFillColor(0)

    #Define a line
    ratline = ROOT.TLine(hdivnom.GetBinLowEdge(1),1,hdivnom.GetBinWidth(1)*hdivnom.GetNbinsX(),1)
    
    #Draw
    histmax = 15.
    histstackmax = 60.
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptStat(0)
    islog = True
    
    tc.cd()
    p21.Draw()
    p21.cd()
    plotMzp(p21,hsbdy)
    CMS_lumi.CMS_lumi(p21,4,13)
    p21.Update()
    uncbands.Draw("e4same")#err bands are 2 sigma bands
    sbfit.Draw("SAME")
    hsbdy.Draw("SAME")

    l21.AddEntry(hsbdy,"DY Jets SB MC","ep")
    l21.AddEntry(sbfit,"2 Param Exp fit","l")
    l21.AddEntry(uncbands,str(sigmabars)+" $\sigma$ uncertainty","f")
    l21.SetBorderSize(0)
    l21.Draw()
    label.Draw()
    p21.Update()
     
    tc.cd()
    p11.Draw()
    p11.cd()
    plotMzp(p11,hsrdy)
    CMS_lumi.CMS_lumi(p11,4,13)
    srdyunc.Draw("e4,same,c")
    srfit.Draw("same")
    hsrdy.Draw("SAME")
    p11.Update()

    l11.AddEntry(hsrdy,"DY Jets SR MC","ep")
    l11.AddEntry(srfit,"2 Param Exp fit","l")
    l11.AddEntry(srdyunc,str(sigmabars)+" $\sigma$ uncertainty\)","f")
    l11.SetBorderSize(0)
    l11.Draw()
    

    label2.Draw()
    p11.Update()

    tc.cd()
    p22.Draw()
    p22.cd()
    #setLogAxis(p22,islog)
    #plotMzp(p22,hsbtt,islog)
    plotMzp(p22,hsbtt)


    CMS_lumi.CMS_lumi(p22,4,13)
    sbttunc.Draw("e3,same,c")
    sbttfit.Draw("SAME")
    hsbtt.Draw("SAME")
    p22.Update()

    l22.AddEntry(hsbtt,"ttbar SB MC","ep")
    l22.AddEntry(sbttfit,"2 Param Exp fit","l")
    l22.AddEntry(sbttunc,str(sigmabars)+" $\sigma$ uncertainty","f")
    l22.SetBorderSize(0)
    l22.Draw()
    p22.Update()


    tc.cd()
    p12.Draw()
    p12.cd()


    l111.AddEntry(hdatsb,"Data SB - full","ep")
    l111.AddEntry(sbdatfit,"2 Param Exp fit","l")
    l111.AddEntry(sbdatunc,str(sigmabars)+" $\sigma$ uncertainty","f")
    l111.SetBorderSize(0)
    l111.AddEntry(hsbdyc,"DYJetsToLL","f")
    l111.AddEntry(hsbttc,"TT","f")
    l111.AddEntry(hsbvvc,"VV","f")
    hsbkg.SetMaximum(60.)
    hsbkg.SetMinimum(0.)

    hsbkg.Draw("HIST")
    xax = hsbkg.GetXaxis()
    yax = hsbkg.GetYaxis()
    xax.SetTitle("M_{Z'}")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 45 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    plotMzp(p12,hdatsb,isData=True)
    sbdatunc.Draw("e3,same,c")
    sbdatfit.Draw("SAME")
    hdatsb.Draw("SAME,E1")
    CMS_lumi.CMS_lumi(p12,4,13)
    p12.Update()
    l111.Draw()
    
    tc.cd()
    p23.Draw()
    p23.cd()
    setLogAxis(p23,islog)
    plotMzp(p23,hsbvv,islog,0.001)
    CMS_lumi.CMS_lumi(p23,4,13)
    sbvvunc.Draw("e4,same,c")
    sbvvfit.Draw("SAME")
    hsbvv.Draw("SAME")

    p23.Update()

    l23.AddEntry(hsbvv,"VV SB MC","ep")
    l23.AddEntry(sbvvfit,"2 Param Exp fit","l")
    l23.AddEntry(sbvvunc,"1 $\sigma$ uncertainty","f")
    l23.SetBorderSize(0)
    l23.Draw()
    p23.Update()

    tc.cd()
    p13.Draw()
    p13.cd()
    alpha.GetYaxis().SetTitle("alpha(M_{Z'})")
    alpha.GetXaxis().SetTitle("M_{Z'}")
    alpha.Draw()


    #print("\n")
    #print("=========== printing some values =============")
    #print("    Value of DY SB Fit in 1500 bin center bin: ",uncbands.GetBinContent(8))
    #print("    Value of tt SB Fit in 1500 bin center bin: ",sbttunc.GetBinContent(8))
    #print("    Value of vv SB Fit in 1500 bin center bin: ",sbvvunc.GetBinContent(8))
    #print("    Value of data SB in 1500 bin center bin  : ",hdatsb.GetBinContent(8))
    #print("\n")
    #print("    Value of DY SB Fit in next bin: ",uncbands.GetBinContent(9))
    #print("    Value of tt SB Fit in next bin: ",sbttunc.GetBinContent(9))
    #print("    Value of vv SB Fit in next bin: ",sbvvunc.GetBinContent(9))
    #print("    Value of data SB in next bin  : ",hdatsb.GetBinContent(9))
    #print("\n")
    

    
    figshapes = go.makeOutFile('Run2_2017_2018','alpha_shapes_'+config.get(systname,syststr)+'_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(figshapes)

    #Stacked Plots Canvas
    tc1 = ROOT.TCanvas("tc1","stacks",1800,800)
    pdstacknom = ROOT.TPad("pdstacknom","normal stack",0,.25,.33,1)
    pdstackfit = ROOT.TPad("pdstackfits","fit stack",.33,.25,.67,1)
    pdstackdat = ROOT.TPad("pdstackdat","data v. data fit",.67,.25,1,1)
    pdrationom = ROOT.TPad("pdrationom","mc stack and data ratio",0,0,0.33,.25)
    pdratiofit = ROOT.TPad("pdratiofit","mc fits and data ratio",0.33,0,0.67,.25)
    pdratiodat = ROOT.TPad("pdratiodat","data and data fit ratio",0.67,0,1,.25)

    tc1.cd()
    pdstacknom.Draw()
    pdstacknom.cd()
    hsbkg.SetMaximum(histstackmax)
    hsbkg.Draw("HIST")
    xax = hsbkg.GetXaxis()
    yax = hsbkg.GetYaxis()
    xax.SetTitle("M_{Z'}")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 200 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    CMS_lumi.CMS_lumi(pdstacknom,4,13)
    hdatsb.Draw("same,e1")

    tc1.cd()
    pdrationom.Draw()
    pdrationom.cd()
    hdivnom.Draw()
    ratline.Draw()

    tc1.cd()
    pdstackfit.Draw()
    pdstackfit.cd()
    hsbkgfit.SetMaximum(histstackmax)
    hsbkgfit.Draw("HIST,C")
    xax = hsbkgfit.GetXaxis()
    yax = hsbkgfit.GetYaxis()
    xax.SetTitle("M_{Z'}")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 200 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    CMS_lumi.CMS_lumi(pdstackfit,4,13)
    hdatsb.Draw("same,e1")

    tc1.cd()
    pdratiofit.Draw()
    pdratiofit.cd()
    hdivfit.Draw()
    ratline.Draw()

    tc1.cd()
    pdstackdat.Draw()
    pdstackdat.cd()
    setLogAxis(pdstackdat,True)
    #sbdatfit.Draw("c,e")
    sbdatunc.SetMinimum(0.1)
    sbdatunc.SetFillColor(ROOT.kRed)
    sbdatunc.Draw("E4")
    sbdatfit.Draw("same,L")
    hdatsb.Draw("same,e")
    CMS_lumi.CMS_lumi(pdstackdat,4,13)

    
    
    stackshapes = go.makeOutFile('Run2_2017_2018','alpha_stacks_'+config.get(systname,syststr)+'_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc1.SaveAs(stackshapes)
    
    #Subtracted Background canvas
    tc2 = ROOT.TCanvas("tc","ratio",1800,800)
    pd112 = ROOT.TPad("pd112","datsb",0,.25,.33,1)
    pd122 = ROOT.TPad("pd122","alpha",.33,.25,.67,1)
    pd132 = ROOT.TPad("pd132","alphaextra",.67,.25,1,1)
    pd221 = ROOT.TPad("pd212","sbratio",0,0,0.33,.25)
    pd232 = ROOT.TPad("pd232","etrapratio",.67,0,1,.25)
    hdatsbsub.Add(sbttunc,-1)
    hdatsbsub.Add(sbvvunc,-1)
    sbdatuncsub = sbdatunc.Clone()
    sbdatuncsub.SetFillColor(0)
    sbdatuncsub.Add(sbttunc,-1)
    sbdatuncsub.Add(sbvvunc,-1)

    #print("Integral of data sub sideband ",hdatsbsub.Integral())
    #print("Integral of data sideband     ",hdatsb.Integral())
    #print("The dy sb value in bin 8         : ",hsbdy.GetBinContent(8))
    #print("The subtacted data value in bin 8: ",sbdatuncsub.GetBinContent(8))
    #print("The dy sb value in bin 9         : ",hsbdy.GetBinContent(9))
    #print("The subtacted data value in bin 9: ",sbdatuncsub.GetBinContent(9))
    #print("The dy sb value in bin 10         : ",hsbdy.GetBinContent(10))
    #print("The subtacted data value in bin 10: ",sbdatuncsub.GetBinContent(10))
    
    tc2.cd()
    pd112.Draw()
    pd112.cd()
    CMS_lumi.CMS_lumi(pd112,4,13)
    hsbdy.SetFillColor(bkgcols[0])
    hsbdy.GetXaxis().SetRangeUser(1500,5000)
    hsbdy.Draw("hist")
    sbdatuncsub.GetXaxis().SetRangeUser(1500,5000)
    sbdatuncsub.Draw("hist,same,c")
    

    #####old stuff#####
    #print("=================doing sb data fit==============")
    #sbdatsubfit = ROOT.expFit(hdatsbsub,"sbl","R0+",1500,2500)
    #sbdatsubunc = ROOT.expFitErrBands(hdatsbsub,"sbl","QR0+",2,1500,2500)
    #sbdatsubunc.SetFillColor(2)
    #sbdatsubunc.SetMarkerSize(0)
    #sbfit.SetFillColor(bkgcols[0])
    #sbfit.SetLineColor(bkgcols[0])
    #sbfit.SetFillStyle(1001)
    #hsbdy.SetFillColor(bkgcols[0])
    #hsbdy.SetLineColor(bkgcols[0])
    #plotMzp(pd112,hdatsbsub,isData=True)
    #sbfit.Draw("SAMEC")
    #CMS_lumi.CMS_lumi(pd112,4,13)
    #hsbdy.Draw("HISTSAME")
    #sbdatsubunc.Draw("e3,same,c")
    #sbdatsubfit.Draw("SAME")
    #hdatsbsub.GetXaxis().SetRangeUser(1500,5000)
    #hdatsbsub.Draw("SAME")

    #lstack = ROOT.TLegend(0.40,0.6,0.93,0.8`)
    #lstack.AddEntry(hdatsbsub,"Data SB, VV tt subtracted","ep")
    #lstack.AddEntry(sbdatsubfit,"2 Param Exp fit","l")
    #lstack.AddEntry(sbdatsubunc,"2 $\sigma$ uncertainty","f")
    #lstack.AddEntry(sbfit,"DY MC Fit","f")
    #lstack.AddEntry(hsbdy,"DY SB MC","f")
    #lstack.SetBorderSize(0)
    #lstack.Draw()
    tc2.cd()

    pd221.Draw()
    pd221.cd()
    hsbdiv = hdatsbsub.Clone()
    hsbdiv.Divide(hdatsbsub,hsbdy)
    hsbdiv.Draw()
    tc2.cd()

    tc2.cd()
    pd122.Draw()
    pd122.cd()
    alpha.GetYaxis().SetTitle("alpha(M_{Z'})")
    alpha.GetXaxis().SetTitle("M_{Z'}")
    alpha.Draw()

    tc2.cd()
    pd132.Draw()
    pd132.cd()
    hsrdy.SetFillColor(bkgcols[0])
    hsrdy.SetLineColor(bkgcols[0])
    hdatsbsub1 = hdatsbsub.Clone()#check errors
    print("=========doing sb extrapolation fit==================")
    extrap  = ROOT.alphaExtrapolation(hsbdy,hsrdy,hdatsbsub1)
    extrphist = ROOT.alphaExtrapolationHist(hsbdy,hsrdy,hdatsbsub1,1)
    hsrdy.GetXaxis().SetRangeUser(1500,5000)
    hsrdy.GetYaxis().SetRangeUser(0,20)
    hsrdy.Draw("HIST")
    extrap.Draw("SAME")
    CMS_lumi.CMS_lumi(pd132,4,13)

    lstack1 = ROOT.TLegend(0.30,0.6,0.93,0.8)
    lstack1.AddEntry(extrap,"DY Predict, alpha*(data SB fit)","p")
    lstack1.AddEntry(hsrdy,"DY MC SR","f")
    if validation:
        hdatvr.SetMarkerStyle(8)
        hdatvr.SetMarkerSize(0.5)
        hdatvr.Draw("SAME")
        lstack1.AddEntry(hdatvr,"VR Data, no subtraction","ep")
    lstack1.SetBorderSize(0)
    lstack1.Draw()

    datavis = go.makeOutFile('Run2_2017_2018','alpha_sub_tester_'+config.get(systname,syststr)+'_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc2.SaveAs(datavis)

    rootOutName = go.makeOutFile('Run2_2017_2018','dy_extraploation'+config.get(systname,syststr),'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    rootFile = ROOT.TFile(rootOutName,"recreate")
    extrphist.Write()
    rootFile.Close()
