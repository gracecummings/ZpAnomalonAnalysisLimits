import ROOT
import glob
import os
import sys
#import gecorg as go
import gecorg_test as go
import numpy as np
import pandas as pd
import configparser

def gatherUpDownHists(inifile,hname,reg,dynormup,dynormnom,dynormdwn):
    emptypt = inifile.Get(hname)
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

    print("norms not sued in this function for hist gathering.")
    hdyup  = systbkgsup.getAddedHist(empty19,"DYJetsToLL",reg,hname)#.Scale(dynormup)
    #hdyup.Scale(dynormup)
    hdydwn = systbkgsdwn.getAddedHist(empty20,"DYJetsToLL",reg,hname)#.Scale(dynormdwn)
    #hdydwn.Scale(dynormdwn)
    httup  = systbkgsup.getAddedHist(empty10,"TT",reg,hname)
    httdwn = systbkgsdwn.getAddedHist(empty11,"TT",reg,hname)
    hzzup  = systbkgsup.getAddedHist(empty12,"ZZTo2L2Q",reg,hname)
    hzzdwn = systbkgsdwn.getAddedHist(empty13,"ZZTo2L2Q",reg,hname)
    hwzup  = systbkgsup.getAddedHist(empty14,"WZTo2L2Q",reg,hname)
    hwzdwn = systbkgsdwn.getAddedHist(empty15,"WZTo2L2Q",reg,hname)
    hvvup  = hzzup.Clone()
    hvvdwn = hzzdwn.Clone()
    hvvup.Add(hwzup)
    hvvdwn.Add(hwzdwn)
    hdy    = bkgs.getAddedHist(empty21,"DYJetsToLL",reg,hname)#.Scale(dynormnom)
    #hdy.Scale(dynormnom)
    htt    = bkgs.getAddedHist(empty16,"TT",reg,hname)
    hzz    = bkgs.getAddedHist(empty17,"ZZTo2L2Q",reg,hname)
    hwz    = bkgs.getAddedHist(empty18,"WZTo2L2Q",reg,hname)
    hvv    = hzz.Clone()
    hvv.Add(hwz)

    histsup = [hdyup,httup,hvvup]
    histdwn = [hdydwn,httdwn,hvvdwn]
    histnom = [hdy,htt,hvv]

    return(histsup,histdwn,histnom)



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

def getGoodPlotRange(listofhists):
    maxlists = [x.GetMaximum() for x in listofhists]
    maxofh = max(maxlists)
    plotmax = maxofh*1.3
    return plotmax

def getRelativeDifference(hist,histnom):
    interest = hist.Integral()
    intnom = histnom.Integral()
    reldif = (interest-intnom)/intnom
    return reldif

def getDeviatedOverNominal(hist,histnom):
    interest = hist.Integral()
    intnom = histnom.Integral()
    devonom = interest/intnom
    return devonom

def getUncertaintyCombination(up,dwn):
    quad = (up*up+dwn*dwn)**(1/2)/2
    mid  = (abs(up)+abs(dwn))/2
    return round(quad,3),round(mid,3)

def getDifferenceSummaryNumberOnly(histup,hist,histdwn,name):
    up = getRelativeDifference(histup,hist)
    dwn = getRelativeDifference(histdwn,hist)
    return up,dwn


def getDifferenceSummary(histup,hist,histdwn,name):
    upnum,dwnnum = getDifferenceSummaryNumberOnly(histup,hist,histdwn,name)
    up = "relative  difference for "+name+" up variation: "+str(round(upnum,3))
    dwn = "relative difference for "+name+" dwn variation: "+str(round(dwnnum,3))
    return up,dwn

def getDeviatedOverNominalSummary(histup,hist,histdwn,name):
    upnum  = getDeviatedOverNominal(histup,hist)
    dwnnum = getDeviatedOverNominal(histdwn,hist)
    up = "Deviated over nominal  for "+name+" up variation: "+str(round(upnum,3))
    dwn = "Devoated over nominal for "+name+" dwn variation: "+str(round(dwnnum,3))
    return up,dwn

def makeDifferenceTPave(histup,hist,histdwn,name):
    #upnum,dwnnum = getDifferenceSummaryNumberOnly(histup,hist,histdwn,name)
    updiv  = getDeviatedOverNominal(histup,hist)
    dwndiv = getDeviatedOverNominal(histdwn,hist)
    #up = "deviated over nom, "+name+" up: "+str(round(upnum,3))
    #dwn = "deviated over nom, "+name+" dwn: "+str(round(dwnnum,3))
    updiv = "up deviated over nominal: "+str(round(updiv,3))
    dwndiv = "dwn deviated over nominal: "+str(round(dwndiv,3))
    lab = ROOT.TPaveText(.2,.73,.9,.85,"NBNDC")
    #lab.AddText(up)
    #lab.AddText(dwn)
    lab.AddText(updiv)
    lab.AddText(dwndiv)
    lab.SetFillColor(0)
    return lab

def makeRateUncertaintyTPave(quadunc,midunc):
    lab = ROOT.TPaveText(.5,.5,.85,.60,"NBNDC")
    #lab.AddText("rate unc = reldiff added in quadrature")
    lab.AddText("rel diff add in quad, div by 2 = "+str(quadunc))
    lab.AddText("haldway between rel diff = "+str(midunc))
    lab.SetFillColor(0)
    return lab


if __name__=='__main__':

    #will replace with command line options
    zptcut  = '100.0'
    hptcut  = '300.0'
    metcut  = '75.0'
    btagwp  = '0.8'
    chan    = 'mumu'
    sigxs   = 1.0
    rebindiv = 2

    limrangelow = 1400
    limrangehigh = 3000


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

    ####Rename and restucture
    htt = newNameAndStructure(htt,"TT",rebindiv,limrangelow,limrangehigh)
    hdy = newNameAndStructure(hdy,"DY",1,limrangelow,limrangehigh)
    hvv = newNameAndStructure(hvv,"VV",rebindiv,limrangelow,limrangehigh)

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
        httup = newNameAndStructure(httup,"TT_"+syst+"Up",rebindiv,limrangelow,limrangehigh)
        httdwn = newNameAndStructure(httdwn,"TT_"+syst+"Down",rebindiv,limrangelow,limrangehigh)
        hvvup = newNameAndStructure(hvvup,"VV_"+syst+"Up",rebindiv,limrangelow,limrangehigh)
        hvvdwn = newNameAndStructure(hvvdwn,"VV_"+syst+"Down",rebindiv,limrangelow,limrangehigh)
        hdyup = newNameAndStructure(hdyup,"DY_"+syst+"Up",1,limrangelow,limrangehigh)
        hdydwn = newNameAndStructure(hdydwn,"DY_"+syst+"Down",1,limrangelow,limrangehigh)
        
        #Debug for background
        #Extra plots
        dynormup = np.load(config.get(syst,'pathup')+'/Run2_161718_dynormalization_'+config.get(syst,'strup')+'_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
        dynormnom = np.load(config.get(syst,'pathnom')+'/Run2_161718_dynormalization_'+config.get(syst,'strnom')+'_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
        dynormdwn = np.load(config.get(syst,'pathdwn')+'/Run2_161718_dynormalization_'+config.get(syst,'strdwn')+'_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]

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

        print(dynormup)
        print(dynormdwn)
        print(dynormnom)
        
        hdyuppt  = systbkgsup.getAddedHist(empty19,"DYJetsToLL","sr","h_h_pt")#.Scale(dynormup)
        #hdyuppt.Scale(dynormup)
        hdydwnpt = systbkgsdwn.getAddedHist(empty20,"DYJetsToLL","sr","h_h_pt")#.Scale(dynormdwn)
        #hdydwnpt.Scale(dynormdwn)
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
        #hdypt.Scale(dynormnom)
        httpt    = bkgs.getAddedHist(empty16,"TT","sr","h_h_pt")
        hzzpt    = bkgs.getAddedHist(empty17,"ZZTo2L2Q","sr","h_h_pt")
        hwzpt    = bkgs.getAddedHist(empty18,"WZTo2L2Q","sr","h_h_pt")
        hvvpt    = hzzpt.Clone()
        hvvpt.Add(hwzpt)

        #for met plots
        hmetsup,hmetsdwn,hmets = gatherUpDownHists(tf1,"h_met","sr",dynormup,dynormnom,dynormdwn)
        
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
        htt.GetYaxis().SetRangeUser(0,getGoodPlotRange([httup,htt,httdwn]))
        hdy.GetYaxis().SetRangeUser(0,getGoodPlotRange([hdyup,hdy,hdydwn]))
        hvv.GetYaxis().SetRangeUser(0,getGoodPlotRange([hvvup,hvv,hvvdwn]))
        hdyuppt.GetYaxis().SetRangeUser(0,getGoodPlotRange([hdyuppt,hdypt,hdydwnpt]))
        httuppt.GetYaxis().SetRangeUser(0,getGoodPlotRange([httuppt,httpt,httdwnpt]))
        hvvuppt.GetYaxis().SetRangeUser(0,getGoodPlotRange([hvvuppt,hvvpt,hvvdwnpt]))
            
        dyupx = hdy.GetXaxis()
        dyupx.SetTitle("Z' Jigsaw Mass Estimator")
        dyuptx = hdyuppt.GetXaxis()
        dyuptx.SetTitle("Higgs Candidate pT")
        ttupx = htt.GetXaxis()
        ttupx.SetTitle("Z' Jigsaw Mass Estimator")
        ttuptx = httuppt.GetXaxis()
        ttuptx.SetTitle("Higgs Candidate pT")
        vvupx = hvv.GetXaxis()
        vvupx.SetTitle("Z' Jigsaw Mass Estimator")
        vvuptx = hvvuppt.GetXaxis()
        vvuptx.SetTitle("Higgs Candidate pT")
        
        hdy.SetStats(0)
        htt.SetStats(0)
        hvv.SetStats(0)
        hdyuppt.SetStats(0)
        httuppt.SetStats(0)
        hvvuppt.SetStats(0)
        
        #Relevant Info
        dyyieldup,dyyielddwn = getDifferenceSummary(hdyup,hdy,hdydwn,"DY")
        ttyieldup,ttyielddwn = getDifferenceSummary(httup,htt,httdwn,"TT")
        vvyieldup,vvyielddwn = getDifferenceSummary(hvvup,hvv,hvvdwn,"VV")

        dyup,dydwn = getDifferenceSummaryNumberOnly(hdyup,hdy,hdydwn,"DY")
        ttup,ttdwn = getDifferenceSummaryNumberOnly(httup,htt,httdwn,"TT")
        vvup,vvdwn = getDifferenceSummaryNumberOnly(hvvup,hvv,hvvdwn,"VV")

        dyupdiv,dydwndiv =getDeviatedOverNominalSummary(hdyup,hdy,hdydwn,"DY")
        ttupdiv,ttdwndiv =getDeviatedOverNominalSummary(httup,htt,httdwn,"TT")
        vvupdiv,vvdwndiv =getDeviatedOverNominalSummary(hvvup,hvv,hvvdwn,"VV")

        dyquad,dymid = getUncertaintyCombination(dyup,dydwn)
        ttquad,ttmid = getUncertaintyCombination(ttup,ttdwn)
        vvquad,vvmid = getUncertaintyCombination(vvup,vvdwn)

        print(dyyieldup)
        print(dyupdiv)
        print(dyyielddwn)
        print(dydwndiv)
        print("  when added in quadrature the uncertainty is {0}.\n   When taking the mid point, the uncertainty is {1}.".format(dyquad,dymid))
        print(ttyieldup)
        print(ttupdiv)
        print(ttyielddwn)
        print(ttdwndiv)
        print("  when added in quadrature the uncertainty is {0}.\n   When taking the mid point, the uncertainty is {1}.".format(ttquad,ttmid))
        print(vvyieldup)
        print(vvupdiv)
        print(vvyielddwn)
        print(vvdwndiv)
        print("  when added in quadrature the uncertainty is {0}.\n    When taking the mid point, the uncertainty is {1}.".format(vvquad,vvmid))
        
        #Strings to print
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
        dydifflab = makeDifferenceTPave(hdyup,hdy,hdydwn,"DY")
        ttdifflab = makeDifferenceTPave(httup,htt,httdwn,"TT")
        vvdifflab = makeDifferenceTPave(hvvup,hvv,hvvdwn,"VV")

        dyptdifflab = makeDifferenceTPave(hdyuppt,hdypt,hdydwnpt,"DY")
        ttptdifflab = makeDifferenceTPave(httuppt,httpt,httdwnpt,"TT")
        vvptdifflab = makeDifferenceTPave(hvvuppt,hvvpt,hvvdwnpt,"VV")

        dyunclab = makeRateUncertaintyTPave(dyquad,dymid)
        ttunclab = makeRateUncertaintyTPave(ttquad,ttmid)
        vvunclab = makeRateUncertaintyTPave(vvquad,vvmid)
        
        dyjiglab = ROOT.TPaveText(.6,.3,.9,.4,"NBNDC")
        dyjiglab.AddText(dystrup)
        dyjiglab.AddText(dystrnom)
        dyjiglab.AddText(dystrdwn)
        dyjiglab.SetFillColor(0)
        dyjiglab.GetLine(0).SetTextColor(ROOT.kRed)
        dyjiglab.GetLine(1).SetTextColor(ROOT.kBlack)
        dyjiglab.GetLine(2).SetTextColor(ROOT.kBlue)
        
        dyptlab = ROOT.TPaveText(.6,.3,.9,.4,"NBNDC")
        dyptlab.AddText(dystruppt)
        dyptlab.AddText(dystrnompt)
        dyptlab.AddText(dystrdwnpt)
        dyptlab.SetFillColor(0)
        dyptlab.GetLine(0).SetTextColor(ROOT.kRed)
        dyptlab.GetLine(1).SetTextColor(ROOT.kBlack)
        dyptlab.GetLine(2).SetTextColor(ROOT.kBlue)
        
        ttjiglab = ROOT.TPaveText(.6,.3,.9,.4,"NBNDC")
        ttjiglab.AddText(ttstrup)
        ttjiglab.AddText(ttstrnom)
        ttjiglab.AddText(ttstrdwn)
        ttjiglab.SetFillColor(0)
        ttjiglab.GetLine(0).SetTextColor(ROOT.kRed)
        ttjiglab.GetLine(1).SetTextColor(ROOT.kBlack)
        ttjiglab.GetLine(2).SetTextColor(ROOT.kBlue)
        
        ttptlab = ROOT.TPaveText(.6,.3,.9,.4,"NBNDC")
        ttptlab.AddText(ttstruppt)
        ttptlab.AddText(ttstrnompt)
        ttptlab.AddText(ttstrdwnpt)
        ttptlab.SetFillColor(0)
        ttptlab.GetLine(0).SetTextColor(ROOT.kRed)
        ttptlab.GetLine(1).SetTextColor(ROOT.kBlack)
        ttptlab.GetLine(2).SetTextColor(ROOT.kBlue)
        
        vvjiglab = ROOT.TPaveText(.6,.3,.9,.4,"NBNDC")
        vvjiglab.AddText(vvstrup)
        vvjiglab.AddText(vvstrnom)
        vvjiglab.AddText(vvstrdwn)
        vvjiglab.SetFillColor(0)
        vvjiglab.GetLine(0).SetTextColor(ROOT.kRed)
        vvjiglab.GetLine(1).SetTextColor(ROOT.kBlack)
        vvjiglab.GetLine(2).SetTextColor(ROOT.kBlue)

        vvptlab = ROOT.TPaveText(.6,.3,.9,.4,"NBNDC")
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
        hdy.Draw('hist')
        hdyup.Draw('histsame')
        hdydwn.Draw('histsame')
        dyjiglab.Draw()
        dydifflab.Draw()
        dyunclab.Draw()
        pdy.Update()
        tcbkg.cd()
        ptt.Draw()
        ptt.cd()
        htt.Draw('hist')
        httup.Draw('histsame')
        httdwn.Draw('histsame')
        ttjiglab.Draw()
        ttdifflab.Draw()
        ttunclab.Draw()
        ptt.Update()
        tcbkg.cd()
        pvv.Draw()
        pvv.cd()
        hvv.Draw('hist')
        hvvup.Draw('histsame')
        hvvdwn.Draw('histsame')
        vvjiglab.Draw()
        vvdifflab.Draw()
        vvunclab.Draw()
        pvv.Update()
        tcbkg.cd()
        ptdy.Draw()
        ptdy.cd()
        hdyuppt.Draw('hist')
        hdypt.Draw('histsame')
        hdydwnpt.Draw('histsame')
        dyptlab.Draw()
        dyptdifflab.Draw()
        ptdy.Update()
        tcbkg.cd()
        pttt.Draw()
        pttt.cd()
        httuppt.Draw('hist')
        httpt.Draw('histsame')
        httdwnpt.Draw('histsame')
        ttptlab.Draw()
        ttptdifflab.Draw()
        pttt.Update()
        tcbkg.cd()
        ptvv.Draw()
        ptvv.cd()
        hvvuppt.Draw('hist')
        hvvpt.Draw('histsame')
        hvvdwnpt.Draw('histsame')
        vvptlab.Draw()
        vvptdifflab.Draw()
        ptvv.Update()
        tcbkg.cd()
        tcbkg.Update()
        
        bkgcomp = go.makeOutFile('Run2_161718_ZllHbbMET',chan+'_bkg_'+syst,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        tcbkg.SaveAs(bkgcomp)
    
        #met debug plot
        metup = hmetsup[0].Clone()
        metup.Add(hmetsup[1])
        metup.Add(hmetsup[2])
        metdwn = hmetsdwn[0].Clone()
        metdwn.Add(hmetsdwn[1])
        metdwn.Add(hmetsdwn[2])
        meth = hmets[0].Clone()
        meth.Add(hmets[1])
        meth.Add(hmets[2])

        metup.SetLineColor(ROOT.kRed)
        metdwn.SetLineColor(ROOT.kBlue)
        meth.SetLineColor(ROOT.kBlack)
        metup.GetYaxis().SetRangeUser(0,getGoodPlotRange([metup,meth,metdwn]))
        meth.GetYaxis().SetRangeUser(0,getGoodPlotRange([metup,meth,metdwn]))
        metdwn.GetYaxis().SetRangeUser(0,getGoodPlotRange([metup,meth,metdwn]))
        metup.GetXaxis().SetTitle("MET")
        metdwn.GetXaxis().SetTitle("MET")
        meth.GetXaxis().SetTitle("MET")
        meth.SetStats(0)
        metup.SetStats(0)
        metdwn.SetStats(0)

        print("nominal met: ",meth.Integral())
        print("met up: : ",metup.Integral())
        print("met dwn: ",metdwn.Integral())
    
        tcmet  = ROOT.TCanvas("tcmet","metcheck",600,400)
        metup.Draw("hist")
        meth.Draw("histsame")
        metdwn.Draw("histsame")

        metcomp = go.makeOutFile('Run2_161718_ZllHbbMET',chan+'_met_'+syst,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        tcmet.SaveAs(metcomp)
