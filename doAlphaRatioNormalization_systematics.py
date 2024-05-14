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
import pickle
import argparse

#tdrstyle.setTDRStyle()
#CMS_lumi.lumi_13TeV = "101.27 fb^{-1}"
#CMS_lumi.writeExtraText = 1
#CMS_lumi.extraText = "Simulation Preliminary"

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

def makeTPadOfIntegrals(hist,fit,fitlowr,fithighr):
    histint = hist.Integral(hist.FindBin(fitlowr),hist.FindBin(fithighr))
    fitint = fit.Integral(fitlowr,fithighr)/hist.GetBinWidth(1)
    intlabel = ROOT.TPaveText(.2,.8,.45,.9,"NBNDC")
    intlabel.AddText("MC Int: "+str(round(histint,2)))
    intlabel.AddText("Fit Int    : "+str(round(fitint,2)))
    intlabel.SetFillColor(0)
    return intlabel

def setLogAxis(pad,islog):
    if islog:
        pad.SetLogy()

def plotMsd(pad,hist,histmax,islog=False,logmin=0.1,isData=False):
    #maxi = hist.GetMaximum()
    #`mr   = round(maxi,0)
    #histmax = mr+mr*0.30
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


def doTemplateFitShifts(nomfit,multi,parnum,name,shiftedparamsin = ROOT.TVector()):
    fitpars = []
    fiterrs = []
    shiftedfits = []
    parvecedit = ROOT.TVector(parnum)
    parsysdict = {"2":"up","0":"down"}
    for par in range(parnum):
        parvecedit[par] = nomfit.GetParameter(par)
        fitpars.append(nomfit.GetParameter(par))
        fiterrs.append(nomfit.GetParError(par))
    for i in range(parnum):
        par = fitpars[i]
        if shiftedparamsin.Norm1() > 0:
            parvecedit[i] = shiftedparamsin[i]
        else:
            parup =par+multi*fiterrs[i]
            parvecedit[i] = parup#change the value of that param's index to shifted param
        if "TT" in name:
            fitup = ROOT.gaus2Fit2SetParsAndErrs("par"+str(i)+parsysdict[str(multi+1)],parvecedit,30,400)
            #here is where we feed the fits to the othersx
        if "VV" in name:
            fitup = ROOT.gausPoly1FitSetParsAndErrs("par"+str(i)+parsysdict[str(multi+1)],parvecedit,30,250)
        if "DY" in name:
            #print(parvecedit[0],parvecedit[1],parvecedit[2],parvecedit[3],parvecedit[4],parvecedit[5])
            fitup = ROOT.poly5mod5FitSetParsAndErrs("par"+str(i)+parsysdict[str(multi+1)],parvecedit,30,225)
        shiftedfits.append(fitup)
        parvecedit[i] = fitpars[i]#reset the fit params

    return fitpars,fiterrs,shiftedfits

def plotShiftedFits(nomfit,hist,shiftedfits,name,plotmax,errhist=None,addname=""):
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
        leg.AddEntry(fit,"par[{0}] {1}".format(i,parsysdict[str(multi+1)]),"l")
        fit.Draw("same")
    nomfit.Draw("same")
    leg.AddEntry(nomfit,"nominal","l")
    hist.Draw("sameE1")
    leg.AddEntry(hist,"fitted hist","pe")
    leg.Draw()
    tcshift.Update()
    shiftfitsfile = go.makeOutFile('Run2_'+yearstr,name+'_'+addname+'_shifted_shapes'+parsysdict[str(multi+1)]+'_'+systr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tcshift.SaveAs(shiftfitsfile)

def plotEnvelopeFits(nomfit,hist,shiftedfits,name,plotmax,errhist=None,addname=""):
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
    shiftcols = [ROOT.kRed,ROOT.kBlue]
    names = ["high","low"]
    for i,fit in enumerate(shiftedfits):
        fit.SetLineColor(shiftcols[i])
        leg.AddEntry(fit,"envoplefit {0}".format(names[i]),"l")
        fit.Draw("same")
    nomfit.Draw("same")
    leg.AddEntry(nomfit,"nominal","l")
    hist.Draw("sameE1")
    leg.AddEntry(hist,"fitted hist","pe")
    leg.Draw()
    tcshift.Update()
    shiftfitsfile = go.makeOutFile('Run2_'+yearstr,name+'_'+addname+'_shifted_shapes'+parsysdict[str(multi+1)]+'_'+systr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tcshift.SaveAs(shiftfitsfile)

def shiftedFitNormalization(listoffits,name,hsbkg,htrtt,htrdy,htrvv,hdatsb):
    fitresults = []
    shiftednorms = {}
    paramcounter = 0
    for fit in listoffits:
        if "DY" in name:
            shiftres = ROOT.totalFitDYTemplateVaried(hsbkg.GetStack().Last(),fit,htrtt.Clone(),htrvv.Clone(),hdatsb.Clone(),"R0+")
        if "TT" in name:
            shiftres = ROOT.totalFitTTTemplateVaried(hsbkg.GetStack().Last(),fit,htrdy.Clone(),htrvv.Clone(),hdatsb.Clone(),"R0+")
        if "VV" in name:
            shiftres = ROOT.totalFitVVTemplateVaried(hsbkg.GetStack().Last(),fit,htrtt.Clone(),htrdy.Clone(),hdatsb.Clone(),"R0+")
        #print("The shifted normalization: ",shiftres[2].GetParameters()[17]),
        fitresults.append(shiftres)
        shiftednorms[name+"_par"+str(paramcounter)] = shiftres[2].GetParameters()[17]
        paramcounter +=1
    return fitresults,shiftednorms

def envelopeFitNormalization(listoffits,name,hsbkg,htrtt,htrdy,htrvv,hdatsb):
    fitresults = []
    shiftednorms = {}
    paramcounter = 0
    for fit in listoffits:
        if "DY" in name:
            shiftres = ROOT.totalFitDYEnvelope(hsbkg.GetStack().Last(),fit,htrtt.Clone(),htrvv.Clone(),hdatsb.Clone(),"R0+")
        if "TT" in name:
            shiftres = ROOT.totalFitTTEnvelope(hsbkg.GetStack().Last(),fit,htrdy.Clone(),htrvv.Clone(),hdatsb.Clone(),"R0+")
        if "VV" in name:
            shiftres = ROOT.totalFitVVEvelope(hsbkg.GetStack().Last(),fit,htrtt.Clone(),htrdy.Clone(),hdatsb.Clone(),"R0+")
        #print("The shifted normalization: ",shiftres[2].GetParameters()[17]),
        fitresults.append(shiftres)
        shiftednorms[name+"_envelope"+str(paramcounter)] = shiftres[2].GetParameters()[11]
        paramcounter +=1
    return fitresults,shiftednorms

def makeNormPickleFile(name,normdict,envnormdict,nomnorm,multi):
    #print(envnormdict)
    normdict[name+"_envup"] = envnormdict[name+"_envelope0"]
    normdict[name+"_envdwn"] = envnormdict[name+"_envelope1"]
    normdict[name+"_nominal"] = nomnorm
    uncf = go.makeOutFile('Run2_'+yearstr,name+'shiftednormsforuncs_'+str(multi)+'_dynormalization_'+systr+'_'+rstr,'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    f = open(uncf,'wb')
    pickle.dump(normdict,f)
    f.close()


ROOT.gROOT.SetBatch()
ROOT.gSystem.CompileMacro("../ZpAnomalonAnalysisUproot/cfunctions/alphafits.C","kfc")
ROOT.gSystem.Load("../ZpAnomalonAnalysisUproot/cfunctions/alphafits_C")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-dir","--directory", type=str,help = "date folder with output")
    #parser.add_argument("-s","--syst", type=str,help = "systematic string")
    parser.add_argument("-sn","--systname", type=str,help = "systematic name ['nominal', 'jec', 'btag', 'unclmet', 'muonid', 'muontrig', 'pdfscale', 'qcdscale', 'prefire', 'jer', 'pileup1', 'pileup2']")
    parser.add_argument("-sd","--systdirection", type=str,help = "systematic direction,: up, dwn, nom")
    parser.add_argument("-v","--validationregion", type=bool,help = "is this a validation region?")
    args = parser.parse_args()

    #Get systematic info
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()
    systname = args.systname
    systclass = 'path'+args.systdirection
    syststr = 'str'+args.systdirection


    #Get the samples
    if not args.directory:
        pathbkg    = config.get(systname,systclass)
        pathdata   = config.get(systname,'pathdata')
    else:
        print("manually putting in a path, hope you know what you're doing")
        pathbkg = args.directory
        pathdata   = config.get(systname,'pathdata')
        
    
    #Command Line options for Selections
    zptcut  = args.zptcut#'150.0'
    hptcut  = args.hptcut#'300.0'
    metcut  = args.metcut#'200.0'
    btagwp  = args.btagwp#'0.8'
    years   = [16,17,18]
    yearstr = go.yearFormatter(years)

    #Plotting Style
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = go.lumiFormatter(years)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Simulation Preliminary"


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

    #systr = 'systnominal_btagnom_muidnom'

    #linesforjeccheck
    #bkgs  = go.backgrounds('analysis_output_ZpAnomalon/2022-06-24',zptcut,hptcut,metcut,btagwp,"alphatest_systnominal_kfnom_btagnom_muiddwn_elidnom_elreconom")
    #data  = go.run2('analysis_output_ZpAnomalon/2022-06-24',zptcut,hptcut,metcut,btagwp,"alphatest_systunclup_btagnom_muidnom")#data uncl up means nothing
    #print(bkgs.bkgs)

    #nominal paths
    #bkgs  = go.backgrounds(pathbkg,zptcut,hptcut,metcut,btagwp,systr)
    #data  = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,'alphat_systnominal_btagnom_muidnom')

    bkgs = go.backgrounds(pathbkg,zptcut,hptcut,metcut,btagwp,"alphat_"+config.get(systname,syststr))
    data = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,"alphat_"+config.get(systname,'strdata'))

    #data  = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,"systnominal_kfnom_btagnom_muidnom_elidnom_elreconom")
    #data  = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,"alphatest_systnominal_btagnom_muidnom")
    #data  = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,systr)

    #print(bkgs.bkgs["DYJetsToLL"][16]["tr"])
    #print(data.data)

    tf1 = ROOT.TFile.Open(bkgs.f17dyjetsb[0])
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
    #htrvv  = htrzz.Clone()
    #htrvv.Add(htrwz)
    print("original zz: ",htrzz.Integral())
    print("original wz: ",htrwz.Integral())

    #xs test
    print("doing xs stuff, things are replaced")
    htrzzxs  = bkgs.getAddedHistXSErr(empty.Clone(),"ZZTo2L2Q","tr","h_h_sd",1)#1 is up, -1 is down
    htrwzxs  = bkgs.getAddedHistXSErr(empty.Clone(),"WZTo2L2Q","tr","h_h_sd",1)
    print("uncertainty zz: ",htrzzxs.Integral())
    print("uncertainty wz: ",htrwzxs.Integral())
    htrvv  = htrzzxs.Clone()
    htrvv.Add(htrwzxs)
    #


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

    plotmax = 125 #550 for btag region .#35 for nominal # 125 for met region
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

    #DY fit debug!
    #dyeyeballbin = htrdy.FindBin(250)
    #dybins = htrdy.GetNbinsX()
    #dybinwidth = htrdy.GetBinWidth(2)
    #hdyshift = ROOT.TH1D("hdyshift","hdyshift",int(dybins),-1*htrdy.GetBinLowEdge(dyeyeballbin),-1*htrdy.GetBinLowEdge(dyeyeballbin)+dybinwidth*dybins)
    #for i in range(dybins+1):
     #   hdyshift.SetBinContent(i,htrdy.GetBinContent(i))
     #   hdyshift.SetBinError(i,htrdy.GetBinError(i))
        #print("(lowedge,bincontent) Good Hist, New hist: ({0},{1}) ({2},{3})".format(htrdy.GetBinLowEdge(i),htrdy.GetBinContent(i),hdyshift.GetBinLowEdge(i),hdyshift.GetBinContent(i)))

    #makes some fits
    print("======Doing the DY total region fit======")
    dyfitrange = 225
    dyfit = ROOT.poly5mod5Fit(htrdy.Clone(),"dylprint","ER0+",30,dyfitrange)
    dyfithist = ROOT.poly5mod5FitErrBands(htrdy.Clone(),"dylprint","QER0+",30,dyfitrange)
    print("======Doing the TT total region fit======")
    ttfit = ROOT.gaus2Fit2(htrtt,"ttl","SQER0+",30,400)
    ttfithist = ROOT.gaus2Fit2ErrBands(htrtt,"ttl","SER0+",30,400)
    print("======Doing the VV total region fit======")
    vvfit = ROOT.gausPoly1Fit(htrvv,"vvl","SER0+",30,250,90,5)
    print("======Doing the Total region fit======")
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


    #Do the systematic shifts, param by param
    multi = 1#-1 for down shifts
    parsysdict = {"2":"up","0":"down"}
    #DY
    #dydecorrshiftedparams                  = ROOT.poly5mod5FitDecorrParamsShifted(htrdy.Clone(),"dyl",multi,"QER0+",30,dyfitrange)#getting shifted parameters after decorrelation
    #dyfitpars,dyfiterrs,shiftedupdyfits    = doTemplateFitShifts(dyfit,multi,6,"DY")#doing the individual shifts without decorr
    #dyfitparsdeco,dyfiterrdeco,sdyfitsdeco = doTemplateFitShifts(dyfit,multi,6,"DY",dydecorrshiftedparams)#doing individual fits with decorrel
    #dyenvelopefits                         = ROOT.poly5mod5FitErrFunctionsGraphs(htrdy.Clone(),"dylforuncs","QER0+",30,dyfitrange)#finding the envelopes
    #dyshiftedbkgfits,dyshiftnorms          = shiftedFitNormalization(sdyfitsdeco,"DY",hsbkg,htrtt.Clone(),htrdy.Clone(),htrvv.Clone(),hdatsb.Clone())#re-deriving the norm with the decor params
    #dyenvbkgfits,dyenvnorms                = envelopeFitNormalization(dyenvelopefits,"DY",hsbkg,htrtt.Clone(),htrdy.Clone(),htrvv.Clone(),hdatsb.Clone())
    #makeNormPickleFile("DY",dyshiftnorms,dyenvnorms,dynormpostfit,multi)
    #TT
    #ttfitpars,ttfiterrs,shiftedupttfits    = doTemplateFitShifts(ttfit,multi,6,"TT",)#6 num of fit pars, non-decorrelated fits
    #ttdecorrshiftedparams                  = ROOT.gaus2Fit2DecorrParamsShifted(htrtt.Clone(),"ttforshiftl",multi,"ER0+")
    #ttfitparsdeco,ttfiterrdeco,sttfitsdeco = doTemplateFitShifts(ttfit,multi,6,"TT",ttdecorrshiftedparams)#doing individual fits with decorrel
    #ttenvelopefits                         = ROOT.gaus2Fit2ErrFunctionsGraphs(htrtt.Clone(),"ttlforuncs","QER0+")#finding the envelopes
    #ttshiftedbkgfits,ttshiftnorms          = shiftedFitNormalization(sttfitsdeco,"TT",hsbkg,htrtt.Clone(),htrdy.Clone(),htrvv.Clone(),hdatsb.Clone())#re-deriving the norm with the decor params
    #ttenvbkgfits,ttenvnorms                = envelopeFitNormalization(ttenvelopefits,"TT",hsbkg,htrtt.Clone(),htrdy.Clone(),htrvv.Clone(),hdatsb.Clone())
    #makeNormPickleFile("TT",ttshiftnorms,ttenvnorms,dynormpostfit,multi)
    #VV
    #vvfitpars,vvfiterrs,shiftedupvvfits = doTemplateFitShifts(vvfit,multi,5,"VV")#5 num of fit pars, non-decorrelated fits

    #Get Some fit info
    dyintlab = makeTPadOfIntegrals(htrdy,dyfit,30,250)
    ttintlab = makeTPadOfIntegrals(htrtt,ttfit,30,400)
    vvintlab = makeTPadOfIntegrals(htrvv,vvfit,30,250)
    
    dychi2f,dyndoff,dyfitpavetext = makeTPadOfFitGOF(dyfit)
    ttchi2f,ttndoff,ttfitpavetext = makeTPadOfFitGOF(ttfit)
    vvchi2f,vvndoff,vvfitpavetext = makeTPadOfFitGOF(vvfit)
    mcnormchi2,mcnormdof,mcnormtpavetext = makeTPadOfFitGOF(bkgfit,normfits=True)
    #normchi2,normdof,normtpavetext = makeTPadOfFitGOF(totnormfit,normfits=True)
    normchi2,normdof,normtpavetext = makeTPadOfFitGOF(sbdatfit,normfits=True)
    
    #Save the normalization
    normfilename = go.makeOutFile('Run2_'+yearstr,'dynormalization_'+config.get(systname,syststr)+'_'+rstr,'.npy',str(zptcut),str(hptcut),str(metcut),str(btagwp))
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
    plotMsd(p11,htrdy,100)#100
    #plotMsd(p11,hdyshift,100)#100
    CMS_lumi.CMS_lumi(p11,4,13)
    dyfit.Draw('same')
    dyleg.Draw()
    dyfitpavetext.Draw()
    dyintlab.Draw()
    p11.Update()
    
    tc.cd()
    p12.Draw()
    p12.cd()
    plotMsd(p12,htrtt,45)#45
    CMS_lumi.CMS_lumi(p12,4,13)
    ttfit.Draw("same")
    ttleg.Draw()
    ttfitpavetext.Draw()
    ttintlab.Draw()
    p12.Update()

    tc.cd()
    p13.Draw()
    p13.cd()
    plotMsd(p13,htrvv,25)#25
    CMS_lumi.CMS_lumi(p13,4,13)
    vvfit.Draw("same")
    vvleg.Draw()
    vvfitpavetext.Draw()
    vvintlab.Draw()
    p13.Update()
    tc.cd()
    
    normshapes = go.makeOutFile('Run2_'+yearstr,'norm_shapes_'+config.get(systname,syststr)+'_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
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


    stackedfit = go.makeOutFile('Run2_'+yearstr,'norm_stackfit_'+config.get(systname,syststr)+'_'+rstr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc1.SaveAs(stackedfit)


    #Plot the systematic shifted hists
    #plotShiftedFits(dyfit,htrdy,shiftedupdyfits,"DY",100,dyfithist)
    #plotShiftedFits(dyfit,htrdy,sdyfitsdeco,"DY",100,dyfithist,"decor")
    #plotShiftedFits(ttfit,htrtt,shiftedupttfits,"TT",35,ttfithist)
    #plotShiftedFits(ttfit,htrtt,sttfitsdeco,"TT",40,ttfithist,"decor")
    #plotShiftedFits(vvfit,htrvv,shiftedupvvfits,"VV",25)
    #plotEnvelopeFits(dyfit,htrdy,dyenvelopefits,"DY",100,dyfithist,"envelopefits")
    #plotEnvelopeFits(ttfit,htrtt,ttenvelopefits,"TT",40,ttfithist,"envelopefits")
    
    debughists = go.makeOutFile('Run2_'+yearstr,'norm_debug_shapes_'+config.get(systname,syststr),'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))

    rootfile = ROOT.TFile(debughists,"recreate")
    savehists = [htrtt,htrdy,htrvv,hdatsb]
    histnames = ["TT","DY","VV","data"]
    for h,hist in enumerate(savehists):
        hist.SetName("h_h_sd_"+histnames[h])
        hist.Write()
    dyfit.Write()
    vvfit.Write()
    ttfit.Write()
    bkgfit.Write()#fit to un-normalized MC
    sbdatfit.Write()#Fit with blinded region, bad for visuals
    totnormfit.Write()#Uses params of fit above to draw full line
    lsbdatfit.Write()#low sideband of fit above
    hsbdatfit.Write()#high sideband of fit 2 above

    rootfile.Close()
