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
CMS_lumi.lumi_13TeV = "101.27 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"

def makeTPadOfFitGOF(fit,normfits = False):
    lims = [0.6,0.3,.9,0.45]
    if normfits:
        lims = [0.17,0.8,0.52,0.9]
    chi2 = fit.GetChisquare()
    ndof = fit.GetNDF()
    fitlabel = ROOT.TPaveText(lims[0],lims[1],lims[2],lims[3],"NBNDC")
    #print("chi 2 ",chi2)
    #print("\Chi^2 = {0} , ndof = {1}".format(round(chi2,2),ndof))
    fitlabel.AddText("\Chi^2 = {0} , ndof = {1}".format(round(chi2,2),ndof))
    fitlabel.AddText("\Chi^2 / ndof = {0}".format(round(chi2/ndof,4)))
    fitlabel.SetFillColor(0)
    return chi2,ndof,fitlabel


def setLogAxis(pad,islog):
    if islog:
        pad.SetLogy()

def plotMsd(pad,hist,islog=False,logmin=0.1,isData=False):
    maxi = hist.GetMaximum()
    mr   = round(maxi,0)
    histmax = mr+mr*0.30
    histmin = 0
    if islog:
        histmax = mr*10
        histmin = logmin
    hist.SetMaximum(histmax)
    hist.SetMinimum(histmin)
    hist.SetMarkerStyle(8)
    hist.SetMarkerSize(0.5)
    if isData:
        hist.SetMarkerColor(ROOT.kBlack)
        drawopts = "SAMEE2"
    else:
        hist.SetMarkerColor(ROOT.kBlue)
        drawopts = "E1"
    xax = hist.GetXaxis()
    yax = hist.GetYaxis()
    xax.SetTitle("Higgs Candidate Soft Drop Mass")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 5 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    
    hist.Draw(drawopts)

def makeRatios(hsbkg,hsdat):
    hsumb = hsbkg.GetStack().Last()
    binlist = np.zeros(hsumb.GetNbinsX()+1)
    fill    = np.zeros(hsumb.GetNbinsX()+1)
    ratiolist = np.zeros(hsumb.GetNbinsX()+1)
    rerrlist = np.zeros(hsumb.GetNbinsX()+1)

    for ibin in range(hsumb.GetNbinsX()+1):#CHECK
        bincen = hsumb.GetBinCenter(ibin)
        bkgmc  = hsumb.GetBinContent(ibin)
        data   = hsdat.GetBinContent(ibin)
        binlist[ibin] = bincen
        if ibin != 0:
            ratiolist[ibin] = -1
            datunc = hsdat.GetBinError(ibin)
            datval = hsdat.GetBinContent(ibin)
            bkgval = hsumb.GetBinContent(ibin)
            bkgerr = hsumb.GetBinError(ibin)
            if bkgmc != 0 and data != 0:
                ratiolist[ibin] = datval/bkgval
                #rerrlist[ibin] = datval/bkgval*sqrt((datunc/data)**2+(bkguncs[hname][ibin-1]/bkgmc)**2)
            if bkgmc == 0:
                ratiolist[ibin] = -1
                rerrlist[ibin] = 0
            else:
                ratiolist[ibin] = -1
                rerrlist[ibin] = 0
        
        #remove underflow bin#hopefuly can get rid of this
        ratiolist = np.delete(ratiolist,0)
        binlist   = np.delete(binlist,0)
        rerrlist  = np.delete(rerrlist,0)


ROOT.gSystem.CompileMacro("../ZpAnomalonAnalysisUproot/cfunctions/alphafits.C","kfc")
ROOT.gSystem.Load("../ZpAnomalonAnalysisUproot/cfunctions/alphafits_C")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-dir","--directory", type=str,help = "date folder with output")
    parser.add_argument("-v","--validationregion", type=bool,help = "is this a validation region?")
    args = parser.parse_args()


    #will replace with command line options
    #pathbkg    = 'BkgInputsNominalJECBtagSyst/'
    #pathdata   = 'DataInputsNominalJECBtagSyst/'
    pathbkg    = args.directory#'pfMETNominal/'
    pathdata   = args.directory#'pfMETNominal/'
    zptcut  = args.zptcut#'150.0'
    hptcut  = args.hptcut#'300.0'
    metcut  = args.metcut#'200.0'
    btagwp  = args.btagwp#'0.8'

    #ranges in question
    lsb = [30,70]
    hsb = [150,400]
    br  = [70,150]
    sr  = [110,150]
    vr  = [55,70]
    
    validation = False
    rstr = "signalblind"
    if validation:
        rstr = "validationblind"

    systr = 'systnominal_btagnom'

    bkgs  = go.backgrounds(pathbkg,zptcut,hptcut,metcut,btagwp,systr)
    data  = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,systr)

    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_h_sd')
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    empty4 = empty.Clone()
    empty5 = empty.Clone()
    empty6 = empty.Clone()
    empty7 = empty.Clone()
    empty8 = empty.Clone()
    empty9 = empty.Clone()

    #Gather basics histograms
    hdatsb = data.getAddedHist(empty9,"sb","h_h_sd")
    htrdy  = bkgs.getAddedHist(empty2,"DYJetsToLL","tr","h_h_sd")
    htrtt  = bkgs.getAddedHist(empty6,"TT","tr","h_h_sd")
    htrzz  = bkgs.getAddedHist(empty7,"ZZTo2L2Q","tr","h_h_sd")
    htrwz  = bkgs.getAddedHist(empty8,"WZTo2L2Q","tr","h_h_sd")
    htrvv  = htrzz.Clone()
    htrvv.Add(htrwz)

    #Background hist styles
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    htrdy.SetFillColor(bkgcols[0])
    htrdy.SetLineColor(bkgcols[0])
    htrtt.SetFillColor(bkgcols[1])
    htrtt.SetLineColor(bkgcols[1])
    htrwz.SetFillColor(bkgcols[2])
    htrwz.SetLineColor(bkgcols[2])
    htrzz.SetFillColor(bkgcols[3])
    htrzz.SetLineColor(bkgcols[3])

    plotmax = 35.
    linemax = plotmax*.8
    labelmin = plotmax*.75

    #stack for no scaling
    hsbkg = ROOT.THStack("hsbkg","")
    hsbkg.Add(htrzz)
    hsbkg.Add(htrwz)
    hsbkg.Add(htrtt)
    hsbkg.Add(htrdy)
    hsbkg.SetMaximum(plotmax)
    hsbkg.SetMinimum(0.0)

    #makes some fits
    dyfit = ROOT.poly5Fit(htrdy,"dyl","QR0+",30,250)
    #ttfit = ROOT.gaus2Fit(htrtt,"ttl","QR0+",30,400)
    ttfit = ROOT.gaus2Fit2(htrtt,"ttl","QR0+",30,400)
    vvfit = ROOT.gausPoly1Fit(htrvv,"vvl","QR0+",30,250,90,5)
    normfits = ROOT.totalFit(hsbkg.GetStack().Last(),htrdy,htrtt,htrvv,hdatsb,"R0+",validation)
    bkgfit = normfits[0]#fit to un-normalized MC
    sbdatfit = normfits[1]#Fit with blinded region, bad for visuals
    totnormfit = normfits[2]#Uses params of fit above to draw full line
    lsbdatfit = normfits[3]#low sideband of fit above
    hsbdatfit = normfits[4]#high sideband of fit 2 above
    dynormprefit = bkgfit.GetParameters()[17] 
    dynormpostfit = totnormfit.GetParameters()[17]

    print("MC only Normalization:  ",dynormprefit)
    print("Data sideband Normalization: ",dynormpostfit)

    #Get Some fit info
    dychi2f,dyndoff,dyfitpavetext = makeTPadOfFitGOF(dyfit)
    ttchi2f,ttndoff,ttfitpavetext = makeTPadOfFitGOF(ttfit)
    vvchi2f,vvndoff,vvfitpavetext = makeTPadOfFitGOF(vvfit)
    mcnormchi2,mcnormdof,mcnormtpavetext = makeTPadOfFitGOF(bkgfit,normfits=True)
    normchi2,normdof,normtpavetext = makeTPadOfFitGOF(totnormfit,normfits=True)
    
    #Save the normalization
    normfilename = go.makeOutFile('Run2_2017_2018','dynormalization_'+systr+'_'+rstr,'.npy',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    normfile     = open(normfilename,'wb')
    np.save(normfile,np.array([dynormpostfit]))
    normfile.close()
    
    #Define some fit styles
    bkgfit.SetLineStyle(6)
    totnormfit.SetLineColor(ROOT.kBlue)
    totnormfit.SetLineStyle(9)

    #Stack for DY Scaling
    hsbkgnorm = ROOT.THStack("hsbkgnorm","")
    htrdyclone = htrdy.Clone()
    htrdyclone.SetFillColor(bkgcols[0])
    htrdyclone.SetLineColor(bkgcols[0])
    htrdyclone.Scale(dynormpostfit)
    hsbkgnorm.Add(htrzz)
    hsbkgnorm.Add(htrwz)
    hsbkgnorm.Add(htrtt)
    hsbkgnorm.Add(htrdyclone)
    hsbkgnorm.SetMaximum(plotmax)
    hsbkgnorm.SetMinimum(0.0)

    #Ratio plots
    #Needs backgrounds stored as an added hist
    #Stacks do not track uncertainties properly
    hbkgaddtionscaled = htrzz.Clone()
    hbkgaddtionscaled.Add(htrwz)
    hbkgaddtionscaled.Add(htrtt)
    hbkgaddtionscaled.Add(htrdyclone)
    hdivscaled = hdatsb.Clone()
    hdivscaled.Divide(hdatsb,hbkgaddtionscaled)
    hbkgaddtion = htrzz.Clone()
    hbkgaddtion.Add(htrwz)
    hbkgaddtion.Add(htrtt)
    hbkgaddtion.Add(htrdy)
    hdiv = hdatsb.Clone()
    hdiv.Divide(hdatsb,hbkgaddtion)

    #labels and legends
    dyleg  = ROOT.TLegend(0.55,0.65,0.9,0.8)
    dyleg.AddEntry(htrdy,"DY","ep")
    dyleg.AddEntry(dyfit,"5th deg poly fit","l")
    dyleg.SetBorderSize(0)
    ttleg  = ROOT.TLegend(0.55,0.65,0.9,0.8)
    ttleg.AddEntry(htrtt,"TT","ep")
    ttleg.AddEntry(ttfit,"2 Gaussian fit","l")
    ttleg.SetBorderSize(0)
    vvleg  = ROOT.TLegend(0.55,0.65,0.9,0.8)
    vvleg.AddEntry(htrvv,"VV","ep")
    vvleg.AddEntry(vvfit,"GausPol1 Fit","l")
    vvleg.SetBorderSize(0)
    stackleg = ROOT.TLegend(0.55,0.40,0.93,0.8)
    stackleg.AddEntry(htrdy,"DYJetsToLL","f")
    stackleg.AddEntry(htrtt,"TT","f")
    stackleg.AddEntry(htrwz,"WZTo2L2Q","f")
    stackleg.AddEntry(htrzz,"ZZTo2L2Q","f")
    stackleg.AddEntry(bkgfit,"Bkg MC fit Prefit - Likelihood","l")
    stackleg.SetBorderSize(0)
    normcomplabel = ROOT.TPaveText(0.17,0.8,0.52,0.9,"NBNDC")
    normcomplabel.AddText("Prefit DY Norm:  {0}".format(round(dynormprefit,4)))
    normcomplabel.AddText("Postfit DY Norm: {0}".format(round(dynormpostfit,4)))
    normcomplabel.SetFillColor(0)
    unnormlabel = ROOT.TPaveText(0.6,0.3,0.93,0.4,"NBNDC")
    unnormlabel.AddText("DY MC w/out SB data norm")
    unnormlabel.SetFillColor(0)
    normlabel = ROOT.TPaveText(0.6,0.3,0.93,0.4,"NBNDC")
    normlabel.AddText("DY MC with SB data norm")
    normlabel.SetFillColor(0)
    brl = ROOT.TLine(br[0],0,br[0],linemax)
    brh = ROOT.TLine(br[1],0,br[1],linemax)
    vrl = ROOT.TLine(vr[0],0,vr[0],linemax)
    srl = ROOT.TLine(sr[0],0,sr[0],linemax)
    srlabel = ROOT.TPaveText(125,labelmin,135,linemax,"NB")
    srlabel.SetFillColor(0)
    srlabel.AddText("SR")
    vrlabel = ROOT.TPaveText(57.5,labelmin,67.5,linemax,"NB")
    vrlabel.SetFillColor(0)
    vrlabel.AddText("VR")
    zrlabel = ROOT.TPaveText(85,labelmin,95,linemax,"NB")
    zrlabel.SetFillColor(0)
    zrlabel.AddText("ZR")
    
    #make some output
    tc = ROOT.TCanvas("tc","shapes",1100,400)
    p11 = ROOT.TPad("p11","dysr",0,0,0.33,1.0)
    p12 = ROOT.TPad("p12","ttsr",0.33,0,0.66,1.0)
    p13 = ROOT.TPad("p13","vvsr",0.66,0,1.0,1.0)


    tc.cd()
    p11.Draw()
    p11.cd()
    plotMsd(p11,htrdy)
    CMS_lumi.CMS_lumi(p11,4,13)
    dyfit.Draw('same')
    dyleg.Draw()
    dyfitpavetext.Draw()
    p11.Update()
    
    tc.cd()
    p12.Draw()
    p12.cd()
    plotMsd(p12,htrtt)
    CMS_lumi.CMS_lumi(p12,4,13)
    ttfit.Draw("same")
    ttleg.Draw()
    ttfitpavetext.Draw()
    p12.Update()

    tc.cd()
    p13.Draw()
    p13.cd()
    plotMsd(p13,htrvv)
    CMS_lumi.CMS_lumi(p13,4,13)
    vvfit.Draw("same")
    vvleg.Draw()
    vvfitpavetext.Draw()
    p13.Update()
    tc.cd()
    
    normshapes = go.makeOutFile('Run2_2017_2018','norm_shapes_'+systr+'_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(normshapes)

    tc1 = ROOT.TCanvas("tc1","stacked",1500,800)
    pd11 = ROOT.TPad("pd11","bkgonly",0,.25,0.5,1.0)
    pd12 = ROOT.TPad("pd12","datfit",0.5,0.25,1.0,1.0)
    pd21 = ROOT.TPad("pd21","bkgonlyratio",0,0,0.5,.25)
    pd22 = ROOT.TPad("pd22","datfitratio",0.5,0,1.0,.25)
    ratline = ROOT.TLine(hdivscaled.GetBinLowEdge(1),1,hdivscaled.GetBinWidth(1)*hdivscaled.GetNbinsX(),1)

    pd11.Draw()
    pd11.cd()
    hsbkg.Draw('HIST')
    xax = hsbkg.GetXaxis()
    yax = hsbkg.GetYaxis()
    xax.SetTitle("Higgs Candidate Soft Drop Mass")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 5 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    CMS_lumi.CMS_lumi(pd11,4,13)
    bkgfit.Draw('SAME')
    totnormfit.Draw("SAME")
    lsbdatfit.Draw("SAME")
    hsbdatfit.Draw("SAME")
    hdatsb.SetMarkerStyle(8)
    hdatsb.SetMarkerSize(0.5)
    hdatsb.SetMarkerColor(ROOT.kBlack)
    hdatsb.Draw("SAME")
    stackleg.AddEntry(sbdatfit,"Fit to Data SB","l")
    stackleg.AddEntry(totnormfit,"Data SB Extrap - Likelihood","l")
    stackleg.AddEntry(sbdatfit,"Data SB","ep")
    stackleg.Draw()
    unnormlabel.Draw()
    mcnormtpavetext.Draw()
    brl.Draw()
    brh.Draw()
    srl.Draw()
    srlabel.Draw()
    zrlabel.Draw()
    if validation:
        vrl.Draw()
        vrlabel.Draw()
    pd11.Update()
    tc1.cd()
    tc1.Update()

    pd21.Draw()
    pd21.cd()
    hdiv.SetMarkerStyle(8)
    hdiv.SetMarkerSize(0.5)
    hdiv.SetMarkerColor(ROOT.kBlack)
    hdiv.GetYaxis().SetRangeUser(0,2)
    hdiv.GetYaxis().SetTitle("data/bkg")
    hdiv.GetYaxis().SetTitleSize(0.15)
    hdiv.GetYaxis().SetTitleOffset(0.3)
    hdiv.GetYaxis().SetLabelSize(0.12)
    hdiv.GetYaxis().SetLabelOffset(0.017)
    hdiv.GetYaxis().SetNdivisions(503)
    hdiv.GetXaxis().SetLabelSize(0.10)
    hdiv.GetXaxis().SetLabelOffset(0.017)

    hdiv.Draw()
    ratline.Draw()
    pd21.Update()
    tc1.cd()
    tc1.Update()


    pd12.Draw()
    pd12.cd()
    hsbkgnorm.Draw('HIST')
    xax = hsbkgnorm.GetXaxis()
    yax = hsbkgnorm.GetYaxis()
    xax.SetTitle("Higgs Candidate Soft Drop Mass")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 5 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.015)
    CMS_lumi.CMS_lumi(pd12,4,13)
    totnormfit.Draw("SAME")
    lsbdatfit.Draw("SAME")
    hsbdatfit.Draw("SAME")
    hdatsb.SetMarkerStyle(8)
    hdatsb.SetMarkerSize(0.5)
    hdatsb.SetMarkerColor(ROOT.kBlack)
    hdatsb.Draw("SAME")
    bkgfit.Draw("SAME")
    stackleg.Draw()
    normcomplabel.Draw()
    normlabel.Draw()
    brl.Draw()
    brh.Draw()
    srl.Draw()
    srlabel.Draw()
    zrlabel.Draw()
    if validation:
        vrl.Draw()    
        vrlabel.Draw()
    p12.Update()
    tc1.cd()
    tc1.Update()

    pd22.Draw()
    pd22.cd()
    hdivscaled.SetMarkerStyle(8)
    hdivscaled.SetMarkerSize(0.5)
    hdivscaled.SetMarkerColor(ROOT.kBlack)
    hdivscaled.GetYaxis().SetRangeUser(0,2)
    hdivscaled.GetYaxis().SetTitle("data/bkg")
    hdivscaled.GetYaxis().SetTitleSize(0.15)
    hdivscaled.GetYaxis().SetTitleOffset(0.3)
    hdivscaled.GetYaxis().SetLabelSize(0.12)
    hdivscaled.GetYaxis().SetLabelOffset(0.017)
    hdivscaled.GetYaxis().SetNdivisions(503)
    hdivscaled.GetXaxis().SetLabelSize(0.10)
    hdivscaled.GetXaxis().SetLabelOffset(0.017)
    #hdivscaled.GetXaxis().SetNdivisions(503)

    hdivscaled.Draw()
    ratline.Draw()
    pd22.Update()
    tc1.cd()
    tc1.Update()


    stackedfit = go.makeOutFile('Run2_2017_2018','norm_stackfit_'+systr+'_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc1.SaveAs(stackedfit)

    
    #ttbarhist = go.makeOutFile('Run2_2017_2018','ttbar_hist','.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))

    #rootfile = ROOT.TFile(ttbarhist,"recreate")
    #htrtt.Write()
    #rootfile.Close()
