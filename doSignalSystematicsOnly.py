import ROOT
import glob
import os
import sys
import gecorg_test as go
import numpy as np
import pandas as pd
import configparser
import argparse

import tdrstyle
import CMS_lumi

tdrstyle.setTDRStyle()
#CMS_lumi.lumi_13TeV = "137.6 fb^{-1}"
CMS_lumi.writeExtraText = 1

def getDeviatedOverNominal(hist,histnom):
    interest = hist.Integral()
    intnom = histnom.Integral()
    devonom = interest/intnom
    return devonom


def makeDifferenceTPave(histup,hist,histdwn):
    updiv  = getDeviatedOverNominal(histup,hist)
    dwndiv = getDeviatedOverNominal(histdwn,hist)
    updiv = "up-deviated/nominal: "+str(round(updiv,3))
    dwndiv = "dwn-deviated/nominal: "+str(round(dwndiv,3))
    lab = ROOT.TPaveText(.55,.55,.9,.65,"NBNDC")
    lab.AddText(updiv)
    lab.AddText(dwndiv)
    lab.SetFillColor(0)
    return lab

def newNameAndStructure(hist,name,rebindiv,limrangelow,limrangehigh):
    hist.Rebin(rebindiv)
    nbins = hist.GetNbinsX()
    binw  = hist.GetBinWidth(1)
    newbins = [limrangelow+x*binw for x in range(int((limrangehigh-limrangelow)/binw))]
    nh = ROOT.TH1F(name,name,len(newbins),limrangelow,limrangehigh)
    for b,le in enumerate(newbins):
        bnum = hist.FindBin(le)
        bincontent = hist.GetBinContent(bnum)
        binerror   = hist.GetBinError(bnum)
        nh.SetBinContent(b+1,bincontent)
        nh.SetBinError(b+1,binerror)

    return nh

def makeBinLowEdges(hist,lastnormbin):
    #this takes the histogram ou want to rebin, and you just give it the last normal bin
    nbins = hist.GetNbinsX()
    binw  = hist.GetBinWidth(1)
    iniedges = [hist.GetBinLowEdge(i) for i in range(nbins+2)]#+2 for the actual highedge
    if lastnormbin not in iniedges:
        print("desired bin edge not possible, try again")
    newedges = [edg for edg in iniedges if (edg <= lastnormbin)]
    newedges.append(iniedges[-1])
    newedges = newedges[1:]
    return np.array(newedges)

def getRelevantPlots(signalclass,signalname,xs,region,hname,rebinval,limrangelow,limrangehigh,empty):
    test18 = signalclass.getAddedHist(signalname,xs,empty.Clone(),region,hname,[18])
    test17 = signalclass.getAddedHist(signalname,xs,empty.Clone(),region,hname,[17])
    test16 = signalclass.getAddedHist(signalname,xs,empty.Clone(),region,hname,[16])
    testclasssum = signalclass.getAddedHist(signalname,xs,empty.Clone(),region,hname)

    test18 = newNameAndStructure(test18,signalname+'2018',rebinval,limrangelow,limrangehigh)
    test17 = newNameAndStructure(test17,signalname+'2017',rebinval,limrangelow,limrangehigh)
    test16 = newNameAndStructure(test16,signalname+'2016',rebinval,limrangelow,limrangehigh)
    testclasssum = newNameAndStructure(testclasssum,signalname+'Run2',rebinval,limrangelow,limrangehigh)

    newbinedges = makeBinLowEdges(test18,2800)
    test18 = test18.Rebin(len(newbinedges)-1,signalname+'2018',newbinedges)
    test17 = test17.Rebin(len(newbinedges)-1,signalname+'2017',newbinedges)
    test16 = test16.Rebin(len(newbinedges)-1,signalname+'2016',newbinedges)
    testclasssum = testclasssum.Rebin(len(newbinedges)-1,signalname+'Run2',newbinedges)
    
    return [test16,test17,test18],testclasssum


def applyStatsUncToSignal(hist,errseries,scale):
    #alpha = 1.- 0.682689492
    for ibin in range(hist.GetNbinsX()+1):
        if ibin == 0:
            continue
        #elif hist.GetBinContent(ibin) > 0:
        binerr = errseries[ibin-1]
        hist.SetBinError(ibin,binerr)
        #else:
        #    binerrbase = ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin))+1,1)-hist.GetBinContent(ibin)
        #    binerr = binerrbase*scale
        #    hist.SetBinError(ibin,binerr)

            
    return hist


def makeSignalInfoDict(sigclass, region,sigxs):
    sigs = sigclass.getPreppedSig(region,sigxs)
    sigdict = {}
    for sig in sigs:
        sigdict[sig["name"]] = sig
    return sigdict


def makeRatioReadable(divadds,ratiostring):
    divadds.GetXaxis().SetTitleSize(0.11)
    divadds.GetXaxis().SetTitleOffset(0.65)
    divadds.GetXaxis().SetLabelSize(0.075)
    divadds.GetYaxis().SetTitle(ratiostring)
    divadds.GetYaxis().SetTitleSize(0.11)
    divadds.GetYaxis().SetTitleOffset(.45)
    divadds.GetYaxis().SetLabelSize(0.08)
    divadds.GetYaxis().SetLabelOffset(0.02)
    divadds.GetYaxis().SetNdivisions(503)
    divadds.SetMinimum(0.)
    divadds.SetMarkerStyle(8)
    divadds.SetMaximum(2.)


def getDeviatedOverNominal(hist,histnom):
    interest = hist.Integral()
    intnom = histnom.Integral()
    devonom = interest/intnom
    return devonom


def plotSignalSystematics(hup,hdwn,hnom,name):
    tc = ROOT.TCanvas("tc","tc",700,600)
    p1 = ROOT.TPad("p1","plot",0,0.3,1.0,1.0)
    p2 = ROOT.TPad("p2","ratio",0,0,1.0,0.3)
    leg = ROOT.TLegend(0.2,0.7,0.55,0.85)
    lab = makeDifferenceTPave(hup,hnom,hdwn)
    CMS_lumi.lumi_13TeV = "137 fb^{-1}"
    #drawing style
    #ROOT.gStyle.SetErrorX(0)
    plotmax = max([hup.GetMaximum(),hdwn.GetMaximum(),hnom.GetMaximum()])*10
    hists = [hup,hdwn,hnom]
    for h in hists:
        h.SetMaximum(plotmax)
    hup.SetLineColor(ROOT.kRed)
    hdwn.SetLineColor(ROOT.kBlue)
    hnom.SetLineColor(ROOT.kBlack)
    leg.SetBorderSize(0)
    leg.AddEntry(hup,name+" up","l")
    leg.AddEntry(hnom,name+" nominal","l")
    leg.AddEntry(hdwn,name+" down","l")
    hdivup = hup.Clone()
    hdivup.Divide(hup,hnom)
    hdivdn = hdwn.Clone()
    hdivdn.Divide(hdwn,hnom)
    makeRatioReadable(hdivup,"deviated/nom")
    makeRatioReadable(hdivdn,"deviated/nom")
    hdivup.SetLineColor(ROOT.kRed)
    hdivdn.SetLineColor(ROOT.kBlue)
    hdivup.SetMarkerColor(ROOT.kRed)
    hdivdn.SetMarkerColor(ROOT.kBlue)
    hdivup.SetMarkerSize(1)
    hdivdn.SetMarkerSize(.8)


    #Draw
    tc.Draw()
    tc.cd()
    p1.Draw()
    p1.cd()
    p1.SetLogy()
    hnom.Draw("hist")
    hup.Draw("histsame")
    hdwn.Draw("histsame")
    leg.Draw()
    lab.Draw()
    CMS_lumi.CMS_lumi(p1,4,13)
    p1.Update()
    tc.cd()
    p2.Draw()
    p2.cd()
    hdivup.Draw()
    hdivdn.Draw("same")
    tc.Update()
    tc.SaveAs(go.makeOutFile('Run2_161718_ZllHbbMET_sigalSystematics',name,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp)))

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

    chan    = 'mumu'
    years   = [16,17,18]
    yearstr = go.yearFormatter(years)
    sigxs   = 1.0
    rebindiv = 2

    limrangelow = 1800
    limrangehigh = 10000

    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()

    #signal       = go.signal(config.get('nominal','pathsignom'),zptcut,hptcut,metcut,btagwp,sigxs,[16,17,18],config.get('nominal','strnom'))

    signal_run2  = go.signal_run2(config.get('nominal','pathsignom'),zptcut,hptcut,metcut,btagwp,sigxs,[16,17,18],config.get('nominal','strnom'))

    #get a histogram
    tftest = ROOT.TFile(signal_run2.sig18sr[0])
    empty = tftest.Get('h_zp_jigm')
    empty.Reset("ICESM")
    emp1 = empty.Clone()
    emp2 = empty.Clone()
    emp3 = empty.Clone()
    emp4 = empty.Clone()


    #old way to get signal
    #if validating, all below must be moved into following for loop
    #sigxs = 1.0
    #siginfo = signal.getPreppedSig('sr',sigxs)
    #signom = makeSignalInfoDict(signal,'sr',sigxs)
    #sigcolors = go.colsFromPalette(siginfo,ROOT.kCMYK)
    #nomsigs = [s["name"] for s in siginfo]
    #for s in nomsigs:
    #    name = signom[s]["name"]
    #     signame = "holder"
    #    if "Tune" in name:
    #        strippedname = name.split("_Tune")[0]
    #        signame = strippedname.replace("-","")
    #    else:
    #        signame = name.replace("-","")
    #    print("------- Looking at signal sample ",signame)

    #    if "Zp5500ND1800NS200" not in signame:
    #        continue

    #Updated Feb 2023 way to get signal, all three years

    #Get the signals
    signalsamps = signal_run2.sigsbyname
    for signalname in signalsamps.keys():
        if "Zp5500-ND1800-NS200" not in signalname:
            continue
        syst = 'systmatic'
    
        test18 = signal_run2.getAddedHist(signalname,1.0,emp1,"sr","h_zp_jigm",[18])
        test17 = signal_run2.getAddedHist(signalname,1.0,emp2,"sr","h_zp_jigm",[17])
        test16 = signal_run2.getAddedHist(signalname,1.0,emp3,"sr","h_zp_jigm",[16])
        testclasssum = signal_run2.getAddedHist(signalname,1.0,emp4,"sr","h_zp_jigm")

        test18 = newNameAndStructure(test18,signalname+'2018',2,limrangelow,limrangehigh)
        test17 = newNameAndStructure(test17,signalname+'2017',2,limrangelow,limrangehigh)
        test16 = newNameAndStructure(test16,signalname+'2016',2,limrangelow,limrangehigh)
        testclasssum = newNameAndStructure(testclasssum,signalname+'Run2',2,limrangelow,limrangehigh)

        newbinedges = makeBinLowEdges(test18,2800)
        test18 = test18.Rebin(len(newbinedges)-1,signalname+'2018',newbinedges)
        test17 = test17.Rebin(len(newbinedges)-1,signalname+'2017',newbinedges)
        test16 = test16.Rebin(len(newbinedges)-1,signalname+'2016',newbinedges)
        testclasssum = testclasssum.Rebin(len(newbinedges)-1,signalname+'Run2',newbinedges)

        #Generic Plotting
        test16.GetXaxis().SetTitle("RJR Zprime Estimator")
        test16.GetYaxis().SetTitle("Events")
        test18.GetXaxis().SetTitle("RJR Zprime Estimator")
        test18.GetYaxis().SetTitle("Events")
        test17.GetXaxis().SetTitle("RJR Zprime Estimator")
        test17.GetYaxis().SetTitle("Events")
        testclasssum.GetXaxis().SetTitle("RJR Zprime Estimator")
        testclasssum.GetYaxis().SetTitle("Events")

        year_strings = ["2016, 36.2 fb^{-1}","2017, 41.5 fb^{1}","2018 59.8 fb^{-1}"]
        year_plots = [test16,test17,test18]


        for syst in systs[1:]:
            print("Checking systematics for: ",syst)
            if "prefire" not in syst:
                print("    have not run full run 2 uncertainties. this is a hard coded message.")
                continue

            CMS_lumi.extraText = "Preliminary, "+syst

            systsigup   = go.signal_run2(config.get(syst,'pathsigup'),zptcut,hptcut,metcut,btagwp,sigxs,[16,17,18],config.get(syst,'strup'))
            systsigdwn  = go.signal_run2(config.get(syst,'pathsigdwn'),zptcut,hptcut,metcut,btagwp,sigxs,[16,17,18],config.get(syst,'strdwn'))

            yearup_plots,upsum = getRelevantPlots(systsigup,signalname,1.0,"sr","h_zp_jigm",2,limrangelow,limrangehigh,empty)
            yeardn_plots,dnsum = getRelevantPlots(systsigdwn,signalname,1.0,"sr","h_zp_jigm",2,limrangelow,limrangehigh,empty)

            #TCanvases
            tcycomp  = ROOT.TCanvas("tcycomp","signals by year",1200,500)
            p16      = ROOT.TPad("p16","2016 sig syst",0.0,0.25,0.33,1.0)
            p17      = ROOT.TPad("p17","2017 sig syst",0.33,0.25,0.67,1.0)
            p18      = ROOT.TPad("p18","2018 sig syst",0.67,0.25,1.0,1.0)
            p16rat   = ROOT.TPad("p16","2016 sig syst rat",0.0,0.0,0.33,0.25)
            p17rat   = ROOT.TPad("p17","2017 sig syst rat",0.33,0.0,0.67,0.25)
            p18rat   = ROOT.TPad("p18","2018 sig syst rat",0.67,0.0,1.0,0.25)
            
            #Draw it
            year_rats  = []
            year_pads  = [p16,p17,p18]
            rat_pads   = [p16rat,p17rat,p18rat]
            
            tcycomp.Draw()
            tcycomp.cd()
            devilabs = []
            ratups = []
            ratdwns = []
            for i,year in enumerate(year_strings):
                #final stuff before plotting
                CMS_lumi.lumi_13TeV = year
                devilabs.append(makeDifferenceTPave(yearup_plots[i],year_plots[i],yeardn_plots[i]))
                divup = yearup_plots[i].Clone()
                divup.Divide(yearup_plots[i],year_plots[i])
                divdn = yeardn_plots[i].Clone()
                divdn.Divide(yeardn_plots[i],year_plots[i])
                makeRatioReadable(divup,"deviated/nom")
                makeRatioReadable(divdn,"deviated/nom")
                divup.SetLineColor(ROOT.kRed)
                divdn.SetLineColor(ROOT.kBlue)
                divup.SetMarkerColor(ROOT.kRed)
                divdn.SetMarkerColor(ROOT.kBlue)
                divup.SetMarkerSize(1)
                divdn.SetMarkerSize(.8)
                ratups.append(divup)
                ratdwns.append(divdn)

                #Draw the stuff
                year_pads[i].Draw()
                year_pads[i].cd()
                year_pads[i].SetLogy()
                maxy = yearup_plots[i].GetMaximum()
                yearup_plots[i].SetMaximum(maxy*10)
                year_plots[i].SetLineColor(ROOT.kBlack)
                yearup_plots[i].SetLineColor(ROOT.kRed)
                yeardn_plots[i].SetLineColor(ROOT.kBlue)
                yearup_plots[i].Draw("hist")
                yeardn_plots[i].Draw("histsame")
                year_plots[i].Draw("histsame")
                devilabs[i].Draw()
                CMS_lumi.CMS_lumi(year_pads[i],4,13)
                year_pads[i].Update()
                tcycomp.cd()
                rat_pads[i].Draw()
                rat_pads[i].cd()
                ratups[i].Draw()
                ratdwns[i].Draw("same")
                tcycomp.cd()

            sigbyyear_comp = go.makeOutFile('Run2_161718_ZllHbbMET_byyear',chan+'_'+signalname+'_'+syst,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
            tcycomp.SaveAs(sigbyyear_comp)

            plotSignalSystematics(upsum,dnsum,testclasssum,signalname+"_"+syst)
        

    #validation plots and canvases
    #Old Nominal Signal
    #hsigori = signom[s]["tfile"].Get("h_zp_jigm")
    #hsig = hsigori.Clone()
    #hsig.Scale(signom[s]["scale"])
    #hsig = applyStatsUncToSignal(hsig,signom[s]["errdf"]["h_zp_jigm"]*signom[s]["scale"],signom[s]["scale"])
    #testmansum = test18.Clone()
    #testmansum.Add(test17)
    #testmansum.Add(test16)
    #testmansum.Rebin(2)
    #hsig.Rebin(2)
    
    #divadds = testclasssum.Clone()
    #divadds.Divide(testclasssum,testmansum)
    #makeRatioReadable(divadds,"classadd/manualadd")
    
    #divolds = testclasssum.Clone()
    #divolds.Divide(testclasssum,hsig)
    #makeRatioReadable(divolds,"Run2/2018-scaled")
    
    #test18.SetLineColor(ROOT.kViolet)
    #test17.SetLineColor(ROOT.kOrange)
    #test16.SetLineColor(ROOT.kTeal)
    #testmansum.SetLineColor(ROOT.kGray)
    #testclasssum.SetLineColor(ROOT.kBlack)
    #hsig.SetLineColor(ROOT.kRed)
    
    #testmansum.GetXaxis().SetTitle("RJR Zprime Estimator")
    #testmansum.GetYaxis().SetTitle("Events")
    #testclasssum.GetXaxis().SetTitle("RJR Zprime Estimator")
    #testclasssum.GetYaxis().SetTitle("Events")
    
    
    #TCanvases

        
#        #The tester one to validate the process.
#        ptest      = ROOT.TPad("pdy","dysr",0,0,0.33,1.0)
#        pverify    = ROOT.TPad("pvv","vvsr",0.33,0.3,0.66,1.0)
#        pverrat    = ROOT.TPad("pvv","vvsr",0.33,0.0,0.66,0.3)
#        precover   = ROOT.TPad("2018reco","2018reco",0.66,0.3,1.0,1.0)
#        precorat   = ROOT.TPad("2018reco","2018reco",0.66,0.0,1.0,0.3)
#        
#        tc.cd()
#        ptest.Draw()
#        ptest.cd()
#        #ptest.SetLogy()
#        CMS_lumi.extraText = "Year-by-year comparison"
#        testmansum.Draw('hist,e')
#        test18.Draw('histsame,e')
#        test17.Draw('histsame,e')
#        test16.Draw('histsame,e')
#        CMS_lumi.CMS_lumi(ptest,4,13)
#        ptest.Update()
#        tc.cd()
#        pverify.Draw()
#        pverify.cd()
#        CMS_lumi.extraText = "Manual versus. Class-based addition"
#        testmansum.Draw('e')
#        testclasssum.Draw('histsame,e')
#        CMS_lumi.CMS_lumi(pverify,4,13)
#        pverify.Update()
#        tc.cd()
#        pverrat.Draw()
#        pverrat.cd()
#        divadds.Draw("pe")
#        tc.cd()
#        precover.Draw()
#        precover.cd()
#        CMS_lumi.extraText = "Run 2 vs. Scaled 2018"
#        testclasssum.Draw('hist,e')
#        #test18.Scale(2.3)
#        test18.Draw("histsame,e")
#        hsig.Draw('histsame,e')
#        CMS_lumi.CMS_lumi(precover,4,13)
#        precover.Update()
#        tc.cd()
#        precorat.Draw()
#        precorat.cd()
#        divolds.Draw("pe")
#        
#        
#        sigcomp = go.makeOutFile('Run2_test_ZllHbbMET',chan+'_'+signalname+'_'+syst,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
#        tc.SaveAs(sigcomp)    
#


    
