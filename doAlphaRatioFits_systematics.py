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

def plotExtrapolations(srdyhist,nomextrap,shiftedfits,dycolor,addname,systr):
    #define the basic plotting styles
    srdyhist.SetFillColor(dycolor)
    srdyhist.SetLineColor(dycolor)
    srdyhist.GetYaxis().SetRangeUser(0,7)
    srdyhist.GetXaxis().SetRangeUser(1500,3000)
    nomextrap.SetLineColor(ROOT.kBlack)#nominal subtractions
    tcsub = ROOT.TCanvas("tcsub","tcsub",800,800)
    leg = ROOT.TLegend(0.6,0.45,0.88,0.80)
    leg.SetBorderSize(0)
    tcsub.Draw()
    tcsub.cd()
    srdyhist.Draw("hist")
    nomextrap.Draw("same")
    leg.AddEntry(srdyhist,"SR DY","f")
    leg.AddEntry(nomextrap,"nominal extrapolation","l")
    shiftcols = go.colsFromPalette(shiftedfits,ROOT.kCMYK)
    for i,key in enumerate(shiftedfits.keys()):
        fit = shiftedfits[key]
        fit.SetLineColor(shiftcols[i])
        fit.Draw("same")
        leg.AddEntry(fit,key,"l")
    leg.Draw()
    nomextrap.Draw("same")
    tcsub.Update()
    #outfile = go.makeOutFile('Run2_'+yearstr,"datasidebandsubtractionmethodcomp",'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    outfile = go.makeOutFile('Run2_'+yearstr,addname+'_shifted_shapes_'+systr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tcsub.SaveAs(outfile)


def getTheShiftedExtrapolations(noma,nomsub,shifta,shiftsub,direc):
    shiftedExtraps = {}
    for key in shiftsub.keys():
        name = "extrap_"+key
        fit = shiftsub[key]
        shiftedExtraps[name] = ROOT.alphaExtrapolationFromFits(fit,noma,name)
    for key in shifta.keys():
        name = "extrap_"+key
        fit = shifta[key]
        shiftedExtraps[name] = ROOT.alphaExtrapolationFromFits(nomsub,fit,name)
    return shiftedExtraps

def plotAlphaExtrapCompMeth(originalalpha,newalpha,addname,systr):
    #define the basic plotting styles
    originalalpha.SetLineColor(ROOT.kBlack)
    newalpha.SetLineColor(ROOT.kRed)#nominal subtractions
    tcsub = ROOT.TCanvas("tcsub","tcsub",800,800)
    leg = ROOT.TLegend(0.6,0.45,0.88,0.80)
    leg.SetBorderSize(0)
    tcsub.Draw()
    tcsub.cd()
    originalalpha.Draw()
    newalpha.Draw("same")
    leg.AddEntry(originalalpha,"original method","l")
    leg.AddEntry(newalpha,"new way","l")
    leg.Draw()
    tcsub.Update()
    outfile = go.makeOutFile('Run2_'+yearstr,"alphaextrapmethodcomp",'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tcsub.SaveAs(outfile)


def plotSubtractedDistributionCompMeth(sbdyhist,sbdatfitfit,originalfromhist,dycolor,addname,systr):
    #define the basic plotting styles
    sbdyhist.SetFillColor(dycolor)
    sbdyhist.SetLineColor(dycolor)
    sbdyhist.GetYaxis().SetRangeUser(0,100)
    sbdyhist.GetXaxis().SetRangeUser(1500,3000)
    sbdatfitfit.SetLineColor(ROOT.kRed)#nominal subtractions
    tcsub = ROOT.TCanvas("tcsub","tcsub",800,800)
    leg = ROOT.TLegend(0.6,0.45,0.88,0.80)
    leg.SetBorderSize(0)
    tcsub.Draw()
    tcsub.cd()
    sbdyhist.Draw("hist")
    originalfromhist.Draw("hist,same,c")
    sbdatfitfit.Draw("same")
    leg.AddEntry(sbdyhist,"SB DY","f")
    leg.AddEntry(originalfromhist,"hist sub way","l")
    leg.AddEntry(sbdatfitfit,"nominal subtraction","l")
    leg.Draw()
    sbdatfitfit.Draw("same")
    tcsub.Update()
    outfile = go.makeOutFile('Run2_'+yearstr,"datasidebandsubtractionmethodcomp",'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tcsub.SaveAs(outfile)

def plotSubtractedDistributions(sbdyhist,sbdatfitfit,shiftedfits,dycolor,addname,systr):
    #define the basic plotting styles
    sbdyhist.SetFillColor(dycolor)
    sbdyhist.SetLineColor(dycolor)
    sbdyhist.GetYaxis().SetRangeUser(0,100)
    sbdyhist.GetXaxis().SetRangeUser(1500,3000)
    sbdatfitfit.SetLineColor(ROOT.kBlack)#nominal subtractions
    tcsub = ROOT.TCanvas("tcsub","tcsub",800,800)
    leg = ROOT.TLegend(0.6,0.45,0.88,0.80)
    leg.SetBorderSize(0)
    tcsub.Draw()
    tcsub.cd()
    sbdyhist.Draw("hist")
    sbdatfitfit.Draw("same")
    leg.AddEntry(sbdyhist,"SB DY","f")
    leg.AddEntry(sbdatfitfit,"nominal subtraction","l")
    shiftcols = go.colsFromPalette(shiftedfits,ROOT.kCMYK)
    for i,key in enumerate(shiftedfits.keys()):
        fit = shiftedfits[key]
        fit.SetLineColor(shiftcols[i])
        fit.Draw("same")
        leg.AddEntry(fit,key,"l")
    leg.Draw()
    sbdatfitfit.Draw("same")
    tcsub.Update()
    #outfile = go.makeOutFile('Run2_'+yearstr,"datasidebandsubtractionmethodcomp",'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    outfile = go.makeOutFile('Run2_'+yearstr,addname+'_shifted_shapes_'+systr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tcsub.SaveAs(outfile)

def plotShiftedAlphas(nomalpha,salphas,systr,addname=''):
    tcshift = ROOT.TCanvas("tcashift","tcashift",800,800)
    leg = ROOT.TLegend(0.2,0.75,0.6,0.9)
    leg.SetBorderSize(0)
    tcshift.cd()
    nomalpha.SetLineColor(ROOT.kBlack)
    nomalpha.Draw()
    leg.AddEntry(nomalpha,"nominal alpha","l")
    shiftcols = go.colsFromPalette(salphas,ROOT.kCMYK)
    for i,key in enumerate(salphas.keys()):
        fit = salphas[key]
        fit.SetLineColor(shiftcols[i])
        fit.Draw("same")
        leg.AddEntry(fit,key,"l")
    nomalpha.Draw("same")
    leg.Draw()
    tcshift.Update()
    shiftfitsfile = go.makeOutFile('Run2_'+yearstr,addname+'_shifted_shapes_'+systr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tcshift.SaveAs(shiftfitsfile)
    
    
def plotShiftedFits(nomfit,hist,shiftedfits,name,direc,plotmax,systr,errhist=None,addname=""):
    tcshift = ROOT.TCanvas("tcshift","tcshift",800,800)
    hist.SetLineColor(ROOT.kBlack)
    hist.SetLineWidth(1)
    leg = ROOT.TLegend(0.6,0.45,0.88,0.80)
    #leg = ROOT.TLegend(0.2,0.2,0.5,0.6)
    leg.SetBorderSize(0)
    tcshift.cd()
    #tcshift.SetLogy()
    nomfit.Draw()
    nomfit.GetYaxis().SetRangeUser(0,plotmax)
    if errhist:
        errhist.SetFillColor(ROOT.kGreen-6)
        errhist.SetMarkerSize(0)
        errhist.Draw("sameCE3")
        leg.AddEntry(errhist,"nominal fit envelope","f")
    shiftcols = go.colsFromPalette(shiftedfits,ROOT.kCMYK)
    for i,fit in enumerate(shiftedfits):
        fit.SetLineColor(shiftcols[i])
        #print("While plotting a fit this is the param ",fit.GetParameter(0))
        #print("While plotting a fit this is the param ",fit.GetParameter(1))
        leg.AddEntry(fit,"par[{0}] {1}".format(i,direc),"l")
        fit.Draw("same")
        #fit.Draw()
    nomfit.Draw("same")
    leg.AddEntry(nomfit,"nominal","l")
    hist.Draw("sameE1")
    leg.AddEntry(hist,"fitted hist","pe")
    leg.Draw()
    tcshift.Update()
    shiftfitsfile = go.makeOutFile('Run2_'+yearstr,name+'_'+addname+'_shifted_shapes'+direc+'_'+systr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tcshift.SaveAs(shiftfitsfile)
    
def getShiftedSubtractions(sbdatfitfit,nomdat,nomtt,nomvv,shiftdata,shifttt,shiftvv,direc):
    newsubdatasbs = {}
    listoffits = []
    for i,fit in enumerate(shiftdata):
        name = "subdatasb_datafit_par"+str(i)+direc
        ssub = ROOT.subtractionFromFits(fit,nomtt,nomvv,name)
        newsubdatasbs[name] = ssub
    for i,fit in enumerate(shifttt):
        name = "subdatasb_ttfit_par"+str(i)+direc
        newsubdatasbs[name] = ROOT.subtractionFromFits(nomdat,fit,nomvv,name)
    for i,fit in enumerate(shiftvv):
        name = "subdatasb_vvfit_par"+str(i)+direc
        ssubvv = ROOT.subtractionFromFits(nomdat,nomtt,fit,name)
        newsubdatasbs[name] = ssubvv
    return newsubdatasbs

def getShiftedAlphaRatios(sbnom,srnom,sbshifts,srshifts,direc):
    newalphas = {}
    #do shifted sideband fits
    for i,fit in enumerate(sbshifts):
        name = "alpha_DY_sb_par"+str(i)+direc
        salpha = ROOT.alphaRatioMakerExpExternalParams(fit,srnom,name)
        newalphas[name] = salpha
    for i,fit in enumerate(srshifts):
        name = "alpha_DY_sr_par"+str(i)+direc
        salpha = ROOT.alphaRatioMakerExpExternalParams(sbnom,fit,name)
        newalphas[name] = salpha
    return newalphas
    

def doExpFitShifts(nomfit,parnum,name,lowr,highr,shiftedparamsin = ROOT.TVector()):
    print("doing the shifted fits!")
    fitpars = []
    fiterrs = []
    shiftedfits = []
    for i in range(parnum):
        parvecedit = ROOT.TVector(parnum)
        for j in range(parnum):
            parvecedit[j] = shiftedparamsin[i][j]
        fit = ROOT.expFitSetParsAndErrs(name+"par"+str(i),parvecedit,lowr,highr)
        shiftedfits.append(fit)
    #parvecedit = ROOT.TVector(parnum)
    #parsysdict = {"2":"up","0":"down"}
    #for par in range(parnum):
    #    parvecedit[par] = nomfit.GetParameter(par)
    #    fitpars.append(nomfit.GetParameter(par))
    #    fiterrs.append(nomfit.GetParError(par))
    #for i in range(parnum):
    #    #print("Looking to shift par ",i)
    #    #print("The shifted pars")
    #    parvecedit[i] = shiftedparamsin[i]
    #    parvecedit.Print()
    #    fitup = ROOT.expFitSetParsAndErrs(name+"par"+str(i),parvecedit,lowr,highr)
    #    shiftedfits.append(fitup)
    #    parvecedit[i] = fitpars[i]#reset the fit params

    return fitpars,fiterrs,shiftedfits

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
    #binmax = hist.GetMaximumBin()
    #xmax   = hist.GetBinCenter(binmax)
    #hmax   = hist.GetMaximum()
    binmax = hist.GetMaximumBin()+1#Bin next to max
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

def castFitIntoHistogram(empty,dicfits):
    histlist = []
    for key in dicfits.keys():
        if "up" in key:
            name = key.replace("up","Up")
        if "dwn" in key:
            name = key.replace("dwn","Down")
        hnew = empty.Clone()
        hnew.SetName(name)
        hnew.SetTitle(name)
        fit = dicfits[key]
        for b in range(hnew.GetNbinsX()+1):
            hnew.SetBinContent(b,fit.Eval(hnew.GetBinCenter(b)))
            hnew.SetBinError(b,0.0)
        histlist.append(hnew)
    return histlist
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-dir","--directory", type=str,help = "date folder with output")
    parser.add_argument("-sn","--systname", type=str,help = "systematic name")
    parser.add_argument("-sd","--systdirection", type=str,help = "systematic direction,: up, dwn, nom")
    #probably should add systematic stuff as a input
    parser.add_argument("-v","--validationregion", type=bool,help = "is this a validation region?")
    parser.add_argument("-normsyst","--normsyst", type=str,help = "normalization systematics?")
    parser.add_argument("-normsystval","--normval", type=float,help = "normalization value used")
    args = parser.parse_args()

    #Years
    years   = [16,17,18]
    yearstr = go.yearFormatter(years)

    #Get systematic info
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()

    #Plotting sytle
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = go.lumiFormatter(years)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Simulation Preliminary"

    #Starting parameters
    zptcut  = str(args.zptcut)#'150.0'
    hptcut  = str(args.hptcut)#'300.0'
    metcut  = str(args.metcut)#'200.0'
    btagwp  = str(args.btagwp)#'0.8'
    dynorm  = 1.
    validation = args.validationregion
    rstr = "signalblind"
    rebindiv = 2
    systname = args.systname
    systclass = 'path'+args.systdirection
    syststr = 'str'+args.systdirection
    sigmabars = 1#1 sigma bands

    #Get the samples
    if not args.directory:
        pathbkg    = config.get(systname,systclass)
        pathdata   = config.get(systname,systclass)
    #datastr = 'systnominal_btagnom'

    else:
        pathbkg = args.directory
        pathdata = args.directory
    
    if validation:
        print("You are trying to do a validation of the alpha method run")
        rstr = "validationblind"
        dynorm = np.load('mumu_2022-03-31_ProperREOIDSF_AlphaMethodExtrap//Run2_161718_dynormalization_systnominal_kfnom_btagnom_muidnom_elidnom_elreconom_signalblind_Zptcut100.0_Hptcut300.0_metcut0.0_btagwp0.8.npy')[0]
        bkgs = go.validation('mumu_2022-03-31_ProperREOIDSF_validationOfAlphaMethodExtrap',zptcut,hptcut,metcut,btagwp,'systnominal_kfnom_btagnom_muidnom_elidnom_elreconom')
        print(bkgs)
        data = go.validation('mumu_2022-03-31_ProperREOIDSF_validationOfAlphaMethodExtrap',zptcut,hptcut,metcut,btagwp,syststr)
        sbstring = "30 < m_{hcand,SD} <= 70, met > 0"
        srstring = "110 < m_{hcand,SD} < 150 met < 75."
        srregionstring = "Validation Region"
    else:
        if len(years) == 2:
            dynorm = np.load(pathbkg+'/Run2_2017_2018_dynormalization_'+config.get(systname,syststr)+'_signalblind_Zptcut'+zptcut+'_Hptcut'+hptcut+'_metcut'+metcut+'_btagwp'+btagwp+'.npy')[0]
        elif len(years) == 3 and not args.normsyst:
            #dynorm = np.load(pathbkg+'/Run2_161718_dynormalization_'+config.get(systname,syststr)+'_signalblind_Zptcut'+zptcut+'_Hptcut'+hptcut+'_metcut'+metcut+'_btagwp'+btagwp+'.npy')[0]
            #print("WARNING: Hard coded in the normalization")
            dynormname = ''
            dynorm = np.load(pathbkg+'/Run2_161718_dynormalization_alphatest_systnominal_kfnom_btagnom_muidnom_elidnom_elreconom_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
            #dynorm = np.load('analysis_output_ZpAnomalon/Run2_161718_dynormalization_systnominal_kfnom_btagnom_muidnom_elidnom_elreconom_signalblind_Zptcut100.0_Hptcut300.0_metcut0.0_btagwp0.8.npy')[0]
            #dynorm = np.load
        elif args.normsyst:
            dynormname = args.normsyst
            dynorm = args.normval
        bkgs = go.backgrounds(pathbkg,zptcut,hptcut,metcut,btagwp,config.get(systname,syststr))
        #data = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,config.get(systname,syststr).replace("_elidnom_elreconom","").replace("_kfnom",""))
        #data = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,config.get(systname,syststr))
        print("WARNING: Hardcdoed data gathering - make sure it if the correct path")
        data = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,"alphatest_systnominal_btagnom_muidnom")
        #data = go.run2("analysis_output_ZpAnomalon/2022-05-17/",zptcut,hptcut,metcut,btagwp,"systnominal_btagnom_muidnom")
        sbstring = "30 < m_{hcand,SD} < 70"
        srstring = "110 <= m_{hcand,SD} < 150"
        srregionstring = "Signal Region"

    #print(bkgs.bkgs["DYJetsToLL"][18]["sb"][1])
    #print(bkgs.bkgs["DYJetsToLL"][18]["sr"][1])
    #print(data.data[18]["sb"][1])
    #print(data.data[18]["vr"][1])
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
    empty22 = empty.Clone()
    empty55 = empty.Clone()
    empty44 = empty.Clone()
              
    if not validation:
        hdatsb = data.getAddedHist(empty9,"sb","h_zp_jigm")
    if validation:
        hdatsb = data.getAddedHistData(empty9,"sb","h_zp_jigm")
        hdatvr = data.getAddedHistData(empty10,"vr","h_zp_jigm")
        hdatvr.Rebin(rebindiv)

    hdatsbsub = hdatsb.Clone()    
    hsbdy = bkgs.getAddedHist(empty,"DYJetsToLL","sb","h_zp_jigm")
    hsrdy = bkgs.getAddedHist(empty2,"DYJetsToLL","sr","h_zp_jigm")
    hsrtt = bkgs.getAddedHist(empty22,"TT","sr","h_zp_jigm")
    hsbtt = bkgs.getAddedHist(empty3,"TT","sb","h_zp_jigm")
    hsbzz = bkgs.getAddedHist(empty4,"ZZTo2L2Q","sb","h_zp_jigm")
    hsbwz = bkgs.getAddedHist(empty5,"WZTo2L2Q","sb","h_zp_jigm")
    hsrzz = bkgs.getAddedHist(empty44,"ZZTo2L2Q","sr","h_zp_jigm")
    hsrwz = bkgs.getAddedHist(empty55,"WZTo2L2Q","sr","h_zp_jigm")
    hsbvv = hsbzz.Clone()
    hsbvv.Add(hsbwz)
    hsrvv = hsrzz.Clone()
    hsrvv.Add(hsrwz)
    

    #rebin
    hsbdy.Rebin(rebindiv)
    hsrdy.Rebin(rebindiv)
    hsbtt.Rebin(rebindiv)
    hsbvv.Rebin(rebindiv)
    hdatsbsub.Rebin(rebindiv)
    hdatsb.Rebin(rebindiv)
    hsrtt.Rebin(rebindiv)
    hsrvv.Rebin(rebindiv)

    #Apply the normalization
    hsbdy.Scale(dynorm)
    hsrdy.Scale(dynorm)

    vvsbmaxbin,vvsbxmax,vvsbymax = getMaxXY(hsbvv)
    ttsbmaxbin,ttsbxmax,ttsbymax = getMaxXY(hsbtt)
    dysbmaxbin,dysbxmax,dysbymax = getMaxXY(hsbdy)
    dysrmaxbin,dysrxmax,dysrymax = getMaxXY(hsrdy)
    datsbmaxbin,datsbxmax,datsbymax = getMaxXY(hdatsb)
    vvsbxmax = 1500#manually to align with the others

    maxvals = [datsbxmax,dysrxmax,dysbxmax,ttsbxmax,vvsbxmax]

    vvsbminbin,vvsbxmin = getMinXY(hsbvv,0.1)#good lower limit
    ttsbminbin,ttsbxmin = getMinXY(hsbtt,0.1)
    dysbminbin,dysbxmin = getMinXY(hsbdy,0.1)#fit OK, changing this does not help
    dysrminbin,dysrxmin = getMinXY(hsrdy,0.1)
    datsbminbin,datsbxmin = getMinXY(hdatsb,0.1)#good lower limit

    print("tt sb limits of fit: ",ttsbxmax,ttsbxmin)
    print("dy sb limits of fit: ",dysbxmax,dysbxmin)
    print("dat sb limits of fit: ",datsbxmax,datsbxmin)
    print("vv sb limits of fit: ",vvsbxmax,vvsbxmin)
    
    ROOT.gSystem.CompileMacro("../ZpAnomalonAnalysisUproot/cfunctions/alphafits.C","kfc")
    ROOT.gSystem.Load("../ZpAnomalonAnalysisUproot/cfunctions/alphafits_C")

    #Do bkg shape fits
    print("================= doing dy sb fit =================")
    sbfit = ROOT.expFit(hsbdy,"dysbl","QRE0+",dysbxmax,dysbxmin)#5000)
    uncbands = ROOT.expFitErrBands(hsbdy,"Qsbl","QER0+",sigmabars,dysbxmax,dysbxmin)#5000)
    dysbfitdecorrparamsup = ROOT.expFitDecorrParamsShiftedUp(hsbdy.Clone(),"dysbshiftfinder","ER0+",dysbxmax,dysbxmin)
    dysbparsdecoup,dysbfiterrdecoup,dysbfitsdecoup = doExpFitShifts(sbfit,2,"dysbup",dysbxmax,dysbxmin,dysbfitdecorrparamsup)
    dysbfitdecorrparamsdn = ROOT.expFitDecorrParamsShiftedDown(hsbdy.Clone(),"dysbshiftfinder","ER0+",dysbxmax,dysbxmin)
    dysbparsdecodn,dysbfiterrdecodn,dysbfitsdecodn = doExpFitShifts(sbfit,2,"dysbdwn",dysbxmax,dysbxmin,dysbfitdecorrparamsdn)
    print("================= doing dy sr fit ==================")
    srdyunc = ROOT.expFitErrBands(hsrdy,"sbl","QER0+",sigmabars,dysrxmax,dysrxmin)#4000)
    srfit = ROOT.expFit(hsrdy,"dysrl","EQR0+",dysrxmax,dysrxmin)#4000)
    dysrfitdecorrparamsup = ROOT.expFitDecorrParamsShiftedUp(hsrdy.Clone(),"dysrshiftfinder","ER0+",dysrxmax,dysrxmin)
    dysrparsdecoup,dysrfiterrdecoup,dysrfitsdecoup = doExpFitShifts(srfit,2,"dysrup",dysrxmax,dysrxmin,dysrfitdecorrparamsup)
    dysrfitdecorrparamsdn = ROOT.expFitDecorrParamsShiftedDown(hsrdy.Clone(),"dysrshiftfinder","ER0+",dysrxmax,dysrxmin)
    dysrparsdecodn,dysrfiterrdecodn,dysrfitsdecodn = doExpFitShifts(srfit,2,"dysrdwn",dysrxmax,dysrxmin,dysrfitdecorrparamsdn)
    print("================== doing tt sb fit ==================")
    sbttfit = ROOT.expFit(hsbtt,"ttsbl","ER0+",ttsbxmax,ttsbxmin)#3000)
    sbttunc = ROOT.expFitErrBands(hsbtt,"sbl","QR0+",sigmabars,ttsbxmax,ttsbxmin)#3000)
    ttsbfitdecorrparamsup = ROOT.expFitDecorrParamsShiftedUp(hsbtt.Clone(),"ttsbshiftfinder","ER0+",ttsbxmax,ttsbxmin)
    ttsbparsdecoup,ttsbfiterrdecoup,ttsbfitsdecoup = doExpFitShifts(sbttfit,2,"ttsbup",ttsbxmax,ttsbxmin,ttsbfitdecorrparamsup)
    ttsbfitdecorrparamsdn = ROOT.expFitDecorrParamsShiftedDown(hsbtt.Clone(),"ttsbshiftfinder","ER0+",ttsbxmax,ttsbxmin)
    ttsbparsdecodn,ttsbfiterrdecodn,ttsbfitsdecodn = doExpFitShifts(sbttfit,2,"ttsbdwn",ttsbxmax,ttsbxmin,ttsbfitdecorrparamsdn)
    print("================== doing VV sb fit ==================")
    sbvvfit = ROOT.expFit(hsbvv,"vvsbl","QER0+",vvsbxmax,vvsbxmin)#3000)
    sbvvunc = ROOT.expFitErrBands(hsbvv,"sbl","QER0+",sigmabars,vvsbxmax,vvsbxmin)#3000)
    vvsbfitdecorrparamsup = ROOT.expFitDecorrParamsShiftedUp(hsbvv.Clone(),"vvsbshiftfinder","ER0+",vvsbxmax,vvsbxmin)
    vvsbparsdecoup,vvsbfiterrdecoup,vvsbfitsdecoup = doExpFitShifts(sbvvfit,2,"vvsbup",vvsbxmax,vvsbxmin,vvsbfitdecorrparamsup)
    vvsbfitdecorrparamsdn = ROOT.expFitDecorrParamsShiftedDown(hsbvv.Clone(),"vvsbshiftfinder","ER0+",vvsbxmax,vvsbxmin)
    vvsbparsdecodn,vvsbfiterrdecodn,vvsbfitsdecodn = doExpFitShifts(sbvvfit,2,"vvsbdwn",vvsbxmax,vvsbxmin,vvsbfitdecorrparamsdn)
    print("================= doing data sb fit =================")
    sbdatfit = ROOT.expFit(hdatsb,"datsbl","ER0+",datsbxmax,datsbxmin)#5000)
    sbdatunc = ROOT.expFitErrBands(hdatsb,"sbl","ER0+",sigmabars,datsbxmax,datsbxmin)#5000)
    datsbfitdecorrparamsup = ROOT.expFitDecorrParamsShiftedUp(hdatsb.Clone(),"datsbshiftfinder","ER0+",datsbxmax,datsbxmin)
    #print("The vector of shifted parameters returned to python")
    #datsbfitdecorrparamsup.Print()
    datsbparsdecoup,datsbfiterrdecoup,datsbfitsdecoup = doExpFitShifts(sbdatfit,2,"datsbup",datsbxmax,datsbxmin,datsbfitdecorrparamsup)
    datsbfitdecorrparamsdn = ROOT.expFitDecorrParamsShiftedDown(hdatsb.Clone(),"datsbshiftfinder","ER0+",datsbxmax,datsbxmin)
    datsbparsdecodn,datsbfiterrdecodn,datsbfitsdecodn = doExpFitShifts(sbdatfit,2,"datsbdwn",datsbxmax,datsbxmin,datsbfitdecorrparamsdn)
    print("================= doing alpha ratio fit =============")
    alpha = ROOT.alphaRatioMakerExp(hsbdy,hsrdy,dysbxmax,dysbxmin,dysrxmax,dysrxmin)
    alphaups = getShiftedAlphaRatios(sbfit,srfit,dysbfitsdecoup,dysrfitsdecoup,"up")
    alphadns = getShiftedAlphaRatios(sbfit,srfit,dysbfitsdecodn,dysrfitsdecodn,"dwn")

    #Get Fit info
    dysbchi2,dysbndof,dysbfitgofpad = makeTPadOfFitGOF(sbfit)
    dysrchi2,dysrndof,dysrfitgofpad = makeTPadOfFitGOF(srfit)
    ttsbchi2,ttsbndof,ttsbfitgofpad = makeTPadOfFitGOF(sbttfit)
    vvsbchi2,vvsbndof,vvsbfitgofpad = makeTPadOfFitGOF(sbvvfit)
    datsbchi2,datsbndof,datsbfitgofpad = makeTPadOfFitGOF(sbdatfit)
    
    #Do the subtractions
    hdatsbsub.Add(sbttunc,-1)
    hdatsbsub.Add(sbvvunc,-1)
    sbdatuncsub = sbdatunc.Clone()
    sbdatuncsub.SetFillColor(0)
    sbdatuncsub.Add(sbttunc,-1)
    sbdatuncsub.Add(sbvvunc,-1)
    subdatafit = ROOT.subtractionFromFits(sbdatfit,sbttfit,sbvvfit,"subtracteddatasbfit")
    shiftedsubsup = getShiftedSubtractions(subdatafit,sbdatfit,sbttfit,sbvvfit,datsbfitsdecoup,ttsbfitsdecoup,vvsbfitsdecoup,"up")
    shiftedsubsdn = getShiftedSubtractions(subdatafit,sbdatfit,sbttfit,sbvvfit,datsbfitsdecodn,ttsbfitsdecodn,vvsbfitsdecodn,"dwn")
    print(shiftedsubsdn)
    print("=========doing sb extrapolation fit==================")
    extrap  = ROOT.alphaExtrapolation(hsbdy,hsrdy,sbdatuncsub,dysbxmax,dysbxmin,dysrxmax,dysrxmin,datsbxmax,datsbxmin)
    extrphist = ROOT.alphaExtrapolationHist(hsbdy,hsrdy,sbdatuncsub,1,dysbxmax,dysbxmin,dysrxmax,dysrxmin,datsbxmax,datsbxmin)
    #alpharver = ROOT.alphaExtrapolationFromFits(subdatafit,alpha,"alphaextrapnohistinput")
    shiftedextrapsup = getTheShiftedExtrapolations(alpha,subdatafit,alphaups,shiftedsubsup,"up")
    shiftedextrapsdn = getTheShiftedExtrapolations(alpha,subdatafit,alphadns,shiftedsubsdn,"dwn")
    sextrpuphists = castFitIntoHistogram(empty,shiftedextrapsup)
    sextrpdnhists = castFitIntoHistogram(empty,shiftedextrapsdn)

    extrphnoerrs = extrphist.Clone()
    for i in range(extrphnoerrs.GetNbinsX()+1):
        extrphnoerrs.SetBinError(i,0.0)
    extrphnoerrs.SetName("extrphistnoerrs")

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

    hdivdat = hdatsb.Clone()
    hdivdat.Divide(hdatsb,sbdatunc)
    hdivdat.SetMarkerStyle(8)
    hdivdat.SetMarkerSize(.5)
    hdivdat.GetYaxis().SetRangeUser(0,2)
    hdivdat.GetYaxis().SetTitle("data/dataFit")
    hdivdat.GetYaxis().SetTitleSize(0.15)
    hdivdat.GetYaxis().SetTitleOffset(0.3)
    hdivdat.GetYaxis().SetLabelSize(0.12)
    hdivdat.GetYaxis().SetLabelOffset(0.017)
    hdivdat.GetYaxis().SetNdivisions(503)
    hdivdat.GetXaxis().SetLabelSize(0.10)
    hdivdat.GetXaxis().SetLabelOffset(0.017)

    hdivdatdysb = sbdatuncsub.Clone()
    hdivdatdysb.Divide(sbdatuncsub,hsbdy)
    hdivdatdysb.SetMarkerStyle(8)
    hdivdatdysb.SetMarkerSize(.5)
    hdivdatdysb.GetYaxis().SetRangeUser(0,2)
    hdivdatdysb.GetYaxis().SetTitle("dataFit/DY")
    hdivdatdysb.GetYaxis().SetTitleSize(0.15)
    hdivdatdysb.GetYaxis().SetTitleOffset(0.3)
    hdivdatdysb.GetYaxis().SetLabelSize(0.12)
    hdivdatdysb.GetYaxis().SetLabelOffset(0.017)
    hdivdatdysb.GetYaxis().SetNdivisions(503)
    hdivdatdysb.GetXaxis().SetLabelSize(0.10)
    hdivdatdysb.GetXaxis().SetLabelOffset(0.017)

    hdivextrapdysr = extrphist.Clone()
    hdivextrapdysr.Divide(extrphist,hsrdy)
    hdivextrapdysr.SetMarkerStyle(8)
    hdivextrapdysr.SetMarkerSize(.5)
    hdivextrapdysr.GetYaxis().SetRangeUser(0,2)
    hdivextrapdysr.GetYaxis().SetTitle("extrap/DY")
    hdivextrapdysr.GetYaxis().SetTitleSize(0.15)
    hdivextrapdysr.GetYaxis().SetTitleOffset(0.3)
    hdivextrapdysr.GetYaxis().SetLabelSize(0.12)
    hdivextrapdysr.GetYaxis().SetLabelOffset(0.017)
    hdivextrapdysr.GetYaxis().SetNdivisions(503)
    hdivextrapdysr.GetXaxis().SetLabelSize(0.10)
    hdivextrapdysr.GetXaxis().SetLabelOffset(0.017)

    if validation:
        hdivextrapvrsr = hdatvr.Clone()
        hdivextrapvrsr.Divide(hdatvr,extrphist)
        hdivextrapvrsr.SetMarkerStyle(8)
        hdivextrapvrsr.SetMarkerSize(.5)
        hdivextrapvrsr.GetYaxis().SetRangeUser(-5,5)
        hdivextrapvrsr.GetYaxis().SetTitle("data/extrap")
        hdivextrapvrsr.GetYaxis().SetTitleSize(0.15)
        hdivextrapvrsr.GetYaxis().SetTitleOffset(0.3)
        hdivextrapvrsr.GetYaxis().SetLabelSize(0.12)
        hdivextrapvrsr.GetYaxis().SetLabelOffset(0.017)
        hdivextrapvrsr.GetYaxis().SetNdivisions(503)
        hdivextrapvrsr.GetXaxis().SetLabelSize(0.10)
        hdivextrapvrsr.GetXaxis().SetLabelOffset(0.017)

        hdivextrapvr = hdatvr.Clone()
        hvrbkgclone = hsrtt.Clone()
        hvrbkgclone.Add(hsrdy)
        hvrbkgclone.Add(hsrvv)
        hdivextrapvr.Divide(hdatvr,hvrbkgclone)
        hdivextrapvr.SetMarkerStyle(8)
        hdivextrapvr.SetMarkerSize(.5)
        hdivextrapvr.GetYaxis().SetRangeUser(-5,5)
        hdivextrapvr.GetYaxis().SetTitle("data/MCbkgs")
        hdivextrapvr.GetYaxis().SetTitleSize(0.15)
        hdivextrapvr.GetYaxis().SetTitleOffset(0.3)
        hdivextrapvr.GetYaxis().SetLabelSize(0.12)
        hdivextrapvr.GetYaxis().SetLabelOffset(0.017)
        hdivextrapvr.GetYaxis().SetNdivisions(503)
        hdivextrapvr.GetXaxis().SetLabelSize(0.10)
        hdivextrapvr.GetXaxis().SetLabelOffset(0.017)
    
    #Make the first TCanvas. This is the shape fits
    tc = ROOT.TCanvas("tc","shapes",1200,800)
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
    l111 = ROOT.TLegend(0.55,0.45,0.9,0.8)
    ldysbcomps = ROOT.TLegend(0.55,0.4,0.9,0.8)
    ldysbcomps.SetBorderSize(0)
    ldysbcompf = ROOT.TLegend(0.55,0.4,0.9,0.8)
    ldysbcompf.SetBorderSize(0)
    ldatfitcomp = ROOT.TLegend(0.55,0.4,0.9,0.8)
    ldatfitcomp.SetBorderSize(0)
    

    #Define the labels
    label = ROOT.TPaveText(.35,.5,.9,.6,"NBNDC")
    label.AddText(sbstring)
    label.AddText("150 < m_{hcand,SD}")
    label.SetFillColor(0)
    label2 = ROOT.TPaveText(.35,.5,.9,.6,"NBNDC")
    label2.AddText(srstring)#higgs mass
    label2.AddText(srregionstring)#higgs mass
    label2.SetFillColor(0)
    labelsbs = ROOT.TPaveText(.5,.35,.9,.35,"NBNDC")
    labelsbs.AddText("Sideband")
    labelsbs.AddText("Stacked MC")
    labelsbs.SetFillColor(0)
    labelsbf = ROOT.TPaveText(.5,.25,.9,.35,"NBNDC")
    labelsbf.AddText("Sideband")
    labelsbf.AddText("Stacked SB MC Fits")
    labelsbf.SetFillColor(0)

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
    dysbfitgofpad.Draw()
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
    dysrfitgofpad.Draw()
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
    ttsbfitgofpad.Draw()
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
    hsbkg.SetMaximum(150.)
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
    yax.SetLabelOffset(0.025)
    plotMzp(p12,hdatsb,isData=True)
    sbdatunc.Draw("e3,same,c")
    sbdatfit.Draw("SAME")
    hdatsb.GetYaxis().SetTitleOffset(0.4)
    hdatsb.Draw("SAME,E1")
    CMS_lumi.CMS_lumi(p12,4,13)
    p12.Update()
    datsbfitgofpad.Draw()
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
    vvsbfitgofpad.Draw()
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
    

    
    figshapes = go.makeOutFile('Run2_'+yearstr,'alpha_shapes_'+config.get(systname,syststr)+'_'+rstr+'_'+dynormname,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
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
    #hsbkg.SetMaximum(histstackmax)
    hsbkg.SetMaximum(150.0)
    hsbkg.Draw("HIST")
    xax = hsbkg.GetXaxis()
    yax = hsbkg.GetYaxis()
    xax.SetTitle("M_{Z'}")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 200 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.025)
    CMS_lumi.CMS_lumi(pdstacknom,4,13)
    hdatsb.Draw("same,e1")
    ldysbcomps.AddEntry(hsbdyc,"DYJetsToLL","f")
    ldysbcomps.AddEntry(hsbttc,"TT","f")
    ldysbcomps.AddEntry(hsbvvc,"VV","f")
    ldysbcomps.AddEntry(hdatsb,"Data SB","ep")
    ldysbcomps.Draw()
    labelsbs.Draw()

    tc1.cd()
    pdrationom.Draw()
    pdrationom.cd()
    hdivnom.Draw()
    ratline.Draw()

    tc1.cd()
    pdstackfit.Draw()
    pdstackfit.cd()
    #hsbkgfit.SetMaximum(histstackmax)
    hsbkgfit.SetMaximum(150)
    hsbkgfit.Draw("HIST,C")
    xax = hsbkgfit.GetXaxis()
    yax = hsbkgfit.GetYaxis()
    xax.SetTitle("M_{Z'}")
    xax.SetTitleSize(0.05)
    xax.SetLabelSize(0.035)
    yax.SetTitle("Events / 200 GeV")
    yax.SetTitleSize(0.05)
    yax.SetLabelSize(0.04)
    yax.SetLabelOffset(0.025)
    CMS_lumi.CMS_lumi(pdstackfit,4,13)
    hdatsb.Draw("same,e1")
    ldysbcompf.AddEntry(sbdyfitc,"DYJetsToLL - fit","f")
    ldysbcompf.AddEntry(sbttfitc,"TT - fit","f")
    ldysbcompf.AddEntry(sbvvfitc,"VV - fit","f")
    ldysbcompf.AddEntry(hdatsb,"Data SB","ep")
    ldysbcompf.Draw()
    labelsbf.Draw()

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
    ldatfitcomp.AddEntry(hdatsb,"Data SB","ep")
    ldatfitcomp.AddEntry(sbdatfit,"Fit to Data SB","l")
    ldatfitcomp.AddEntry(sbdatunc,"1 $\sigma$ uncertainty","f")
    ldatfitcomp.Draw()

    tc1.cd()
    pdratiodat.Draw()
    pdratiodat.cd()
    hdivdat.Draw()
    ratline.Draw()
    
    
    stackshapes = go.makeOutFile('Run2_'+yearstr,'alpha_stacks_'+config.get(systname,syststr)+'_'+rstr+'_'+dynormname,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc1.SaveAs(stackshapes)
    
    #Subtracted Background canvas
    tc2 = ROOT.TCanvas("tc","ratio",1800,800)
    pd112 = ROOT.TPad("pd112","datsb",0,.25,.33,1)
    pd122 = ROOT.TPad("pd122","alpha",.33,.25,.67,1)
    pd132 = ROOT.TPad("pd132","alphaextra",.67,.25,1,1)
    pd221 = ROOT.TPad("pd212","sbratio",0,0,0.33,.25)
    pd232 = ROOT.TPad("pd232","etrapratio",.67,0,1,.25)

    tc2.cd()
    pd112.Draw()
    pd112.cd()
    hsbdy.SetFillColor(bkgcols[0])
    hsbdy.SetLineColor(bkgcols[0])
    hsbdy.GetYaxis().SetRangeUser(0,100)
    #hsbdy.GetXaxis().SetRangeUser(1500,5000)
    hsbdy.Draw("hist")
    #sbdatuncsub.GetXaxis().SetRangeUser(1500,5000)
    sbdatuncsub.Draw("hist,same,c")
    CMS_lumi.CMS_lumi(pd112,4,13)
    lsubcomp = ROOT.TLegend(0.30,0.6,0.93,0.8)
    lsubcomp.SetBorderSize(0)
    lsubcomp.AddEntry(sbdatuncsub,"SB Data fit w/ tt+vv sub","l")
    lsubcomp.AddEntry(hsbdy,"DY MC SB","f")
    lsubcomp.Draw()

    tc2.cd()

    pd221.Draw()
    pd221.cd()
    hdivdatdysb.Draw()
    ratline.Draw()
    tc2.cd()

    if not validation:
        pd232.Draw()
        pd232.cd()
        hdivextrapdysr.Draw()
        ratline.Draw()
        tc2.cd()

    if validation:
        pd232.Draw()
        pd232.cd()
        hdivextrapdysr.Draw()
        #hdivextrapvr.SetMarkerColor(ROOT.kRed)
        #hdivextrapvr.Draw("SAME")
        #hdivextrapvrsr.Draw()
        ratline.Draw()
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
    hsrdy.GetYaxis().SetRangeUser(0,10)

    #for b,bb in
    #hsrdy.GetXaxis().SetRangeUser(1500,5000)
    #hsrdy.GetYaxis().SetRangeUser(0,20)
    hsrdy.Draw("HIST")
    extrap.Draw("SAME")
    CMS_lumi.CMS_lumi(pd132,4,13)

    lstack1 = ROOT.TLegend(0.30,0.6,0.93,0.8)
    lstack1.AddEntry(extrap,"DY Predict, alpha*(data SB fit)","p")
    lstack1.AddEntry(hsrdy,"DY MC SR","f")
    if validation:
        hsrbkg = ROOT.THStack("hsrbkg","")
        hsrbkg.Add(hsrdy)
        hsrtt.SetFillColor(bkgcols[1])
        hsrtt.SetLineColor(bkgcols[1])
        hsrvv.SetFillColor(bkgcols[2])
        hsrvv.SetLineColor(bkgcols[2])
        hsrbkg.Add(hsrtt)
        hsrbkg.Add(hsrvv)
        hsrbkg.Draw("HIST")
        extrap.Draw("SAME")
        hdatvr.SetMarkerStyle(8)
        hdatvr.SetMarkerSize(0.5)
        hdatvr.Draw("SAME")
        lstack1.AddEntry(hsrtt,"TT MC SR","f")
        lstack1.AddEntry(hsrvv,"VV MC SR","f")
        lstack1.AddEntry(hdatvr,"VR Data, no subtraction","ep")
        CMS_lumi.CMS_lumi(pd132,4,13)
    lstack1.SetBorderSize(0)
    lstack1.Draw()

    datavis = go.makeOutFile('Run2_'+yearstr,'alpha_sub_tester_'+config.get(systname,syststr)+'_'+rstr+'_'+dynormname,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc2.SaveAs(datavis)

    plotShiftedFits(sbfit,hsbdy,dysbfitsdecoup,"DYSB","up",90,config.get(systname,syststr),uncbands,"decorrelatedShifts")
    plotShiftedFits(sbfit,hsbdy,dysbfitsdecodn,"DYSB","dwn",90,config.get(systname,syststr),uncbands,"decorrelatedShifts")
    plotShiftedFits(srfit,hsrdy,dysrfitsdecoup,"DYSR","up",5,config.get(systname,syststr),srdyunc,"decorrelatedShifts")
    plotShiftedFits(srfit,hsrdy,dysrfitsdecodn,"DYSR","dwn",5,config.get(systname,syststr),srdyunc,"decorrelatedShifts")
    plotShiftedFits(sbttfit,hsbtt,ttsbfitsdecoup,"TTSB","up",55,config.get(systname,syststr),sbttunc,"decorrelatedShifts")
    plotShiftedFits(sbttfit,hsbtt,ttsbfitsdecodn,"TTSB","dwn",55,config.get(systname,syststr),sbttunc,"decorrelatedShifts")
    plotShiftedFits(sbvvfit,hsbvv,vvsbfitsdecoup,"VVSB","up",5,config.get(systname,syststr),sbvvunc,"decorrelatedShifts")
    plotShiftedFits(sbvvfit,hsbvv,vvsbfitsdecodn,"VVSB","dwn",5,config.get(systname,syststr),sbvvunc,"decorrelatedShifts")
    plotShiftedFits(sbdatfit,hdatsb,datsbfitsdecoup,"dataSB","up",150,config.get(systname,syststr),sbdatunc,"decorrelatedShifts")
    plotShiftedFits(sbdatfit,hdatsb,datsbfitsdecodn,"dataSB","dwn",150,config.get(systname,syststr),sbdatunc,"decorrelatedShifts")
    plotShiftedAlphas(alpha,alphaups,config.get(systname,syststr),"alphaups")
    plotShiftedAlphas(alpha,alphadns,config.get(systname,syststr),"alphadwns")
    plotSubtractedDistributions(hsbdy,subdatafit,shiftedsubsup,bkgcols[0],"subtractedDataDistsUp",config.get(systname,syststr))
    plotSubtractedDistributions(hsbdy,subdatafit,shiftedsubsdn,bkgcols[0],"subtractedDataDistsDwn",config.get(systname,syststr))
    plotExtrapolations(hsrdy,extrap,shiftedextrapsup,bkgcols[0],"shiftedExtrapsUp",config.get(systname,syststr))
    plotExtrapolations(hsrdy,extrap,shiftedextrapsdn,bkgcols[0],"shiftedExtrapsDwn",config.get(systname,syststr))
    #plotSubtractedDistributionCompMeth(hsbdy,subdatafit,sbdatuncsub,bkgcols[0],"method comparison",config.get(systname,syststr))
    #plotAlphaExtrapCompMeth(extrap,alpharver,"methodcomp",config.get(systname,syststr))

    rootOutName = go.makeOutFile('Run2_'+yearstr,'dy_extraploation'+config.get(systname,syststr)+'_'+dynormname,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    rootFile = ROOT.TFile(rootOutName,"recreate")
    extrphist.Write()
    extrphnoerrs.Write()
    savehists = [hsbdy,hsrdy,hsbtt,hsbvv,hdatsb,hdatsbsub]
    savehistnames = ["DYSB","DYSR","TT","VV","data","datasubtact"]
    for h,hist in enumerate(savehists):
        hist.SetName("h_zp_mjig_"+savehistnames[h])
        hist.Write()
    for h,hist in enumerate(sextrpuphists):
        hist.Write()
        sextrpdnhists[h].Write()
    sbfit.Write()
    srfit.Write()
    sbttfit.Write()
    sbvvfit.Write()
    sbdatfit.Write()
    sbdatuncsub.Write()
    alpha.Write()
    rootFile.Close()
