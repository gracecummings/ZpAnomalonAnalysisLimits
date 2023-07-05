import ROOT
import glob
import os
import sys
#import gecorg as go
import gecorg_test as go
import numpy as np
import pandas as pd
import configparser
import tdrstyle
import CMS_lumi

tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = "137.6 fb^{-1}"
CMS_lumi.writeExtraText = 1

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

def makeRatioReadable(divadds,ratiostring):
    divadds.GetXaxis().SetTitleSize(0.11*2)
    divadds.GetXaxis().SetTitleOffset(0.65)
    divadds.GetXaxis().SetLabelSize(0.075*2)
    divadds.GetYaxis().SetTitle(ratiostring)
    divadds.GetYaxis().SetTitleSize(0.11*2)
    divadds.GetYaxis().SetTitleOffset(.45)
    divadds.GetYaxis().SetLabelSize(0.08*2)
    divadds.GetYaxis().SetLabelOffset(0.02)
    divadds.GetYaxis().SetNdivisions(503)
    divadds.SetMinimum(0.)
    divadds.SetMarkerStyle(8)
    divadds.SetMaximum(2.)

def makeRatioPlots(ul,dl,nl):
    ratsup = []
    ratsdn = []
    for i,hist in enumerate(nl):
        nom = hist
        up  = ul[i]
        dn  = dl[i]
        ratup = up.Clone()
        ratdn = dn.Clone()
        ratup.Divide(up,nom)
        ratdn.Divide(dn,nom)
        makeRatioReadable(ratup,"up/nom")
        makeRatioReadable(ratdn,"dn/nom")
        ratup.SetLineColor(ROOT.kRed)
        ratdn.SetLineColor(ROOT.kBlue)
        ratup.SetMarkerColor(ROOT.kRed)
        ratdn.SetMarkerColor(ROOT.kBlue)
        ratup.SetMarkerSize(1)
        ratdn.SetMarkerSize(.8)
        ratsup.append(ratup)
        ratsdn.append(ratdn)

    return(ratsup,ratsdn)


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
    #maxlists = [x.GetMaximum() for x in listofhists]
    maxlists = [x.GetBinContent(x.FindBin(x.GetMaximum())+1) for x in listofhists]
    maxofh = max(maxlists)
    plotmax = maxofh*2
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
    lab = ROOT.TPaveText(.4,.6,.9,.75,"NBNDC")
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

if __name__=='__main__':

    #will replace with command line options
    zptcut  = '100.0'
    hptcut  = '300.0'
    metcut  = '75.0'
    btagwp  = '0.8'
    chan    = 'mumu'
    sigxs   = 1.0
    rebindiv = 2

    limrangelow = 1800
    limrangehigh = 10000


    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()

    ####load in the files with the nominal distributions
    bkgs = go.backgrounds('mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref',zptcut,hptcut,metcut,btagwp,'systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco')
    sig  = go.signal(config.get('nominal','pathsignom'),zptcut,hptcut,metcut,btagwp,sigxs,systr=config.get('nominal','strnom'))
    #dyEst = ROOT.TFile(config.get('nominal','pathnom')+'/Run2_161718_dy_extraploationalphat_'+config.get('nominal','strnom')+'__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    #dyEst = ROOT.TFile('analysis_output_ZpAnomalon/2023-01-23/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    dyEst = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-09/alpha_method_ttbar_likli/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')

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
    hdy    = dyEst.Get("extrphistnoerrs").Clone()
    htt = bkgs.getAddedHist(empty1,"TT","sr","h_zp_jigm")
    hzz  = bkgs.getAddedHist(empty2,"ZZTo2L2Q","sr","h_zp_jigm")
    print("hzz integral: ",hzz.Integral())
    hwz  = bkgs.getAddedHist(empty3,"WZTo2L2Q","sr","h_zp_jigm")
    print("hwz integral: ",hwz.Integral())
    hdytr= bkgs.getAddedHist(empty41,"DYJetsToLL","sr","h_zp_jigm")
    hvv  = hzz.Clone()
    hvv.Add(hwz)

    ####Troubleshooting
    httchecksb = bkgs.getAddedHist(empty42,"TT","sb","h_zp_jigm")
    hzzchecksb  = bkgs.getAddedHist(empty43,"ZZTo2L2Q","sb","h_zp_jigm")
    hwzchecksb  = bkgs.getAddedHist(empty44,"WZTo2L2Q","sb","h_zp_jigm")
    hdytrchecksb= bkgs.getAddedHist(empty45,"DYJetsToLL","sb","h_zp_jigm")

    ####Rename and restucture
    htt = newNameAndStructure(htt,"TT",1,limrangelow,limrangehigh)
    hdy = newNameAndStructure(hdy,"DY",1,limrangelow,limrangehigh)
    hvv = newNameAndStructure(hvv,"VV",rebindiv,limrangelow,limrangehigh)

    newbinedges = makeBinLowEdges(hvv,2800)
    #print(newbinedges)


    hvv = hvv.Rebin(len(newbinedges)-1,"VV",newbinedges)
    htt = htt.Rebin(len(newbinedges)-1,"TT",newbinedges)
    hdy = hdy.Rebin(len(newbinedges)-1,"DY",newbinedges)


    #for syst in systs:
    #    if "nominal" == syst:
    #        continue
    #    if "qcdscale" not in syst:
    #        continue
    #    print("------- Looking at systematic ",syst)
    #systbkgsup  = go.backgrounds(config.get(syst,'pathup'),zptcut,hptcut,metcut,btagwp,config.get(syst,'strup'))
    #systbkgsdwn = go.backgrounds(config.get(syst,'pathdwn'),zptcut,hptcut,metcut,btagwp,config.get(syst,'strdwn'))

    syst = "vvxsunc"
    
    #if rebindiv == 2:
    #dyEstup     = ROOT.TFile('analysis_output_ZpAnomalon/2023-01-23/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    dyEstup     = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-09/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_vvuncup_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    dyEstdwn    = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-09/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_vvuncdown_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')


    ####Prepping holders####
    empty4 = empty.Clone()
    empty5 = empty.Clone()
    empty6 = empty.Clone()
    empty7 = empty.Clone()
    empty8 = empty.Clone()
    empty9 = empty.Clone()
        
    ####Gathering the Systematic
    #Background
    httup  = bkgs.getAddedHist(empty4,"TT","sr","h_zp_jigm")
    httdwn = bkgs.getAddedHist(empty5,"TT","sr","h_zp_jigm")
    hzzup  = bkgs.getAddedHistXSErr(empty6,"ZZTo2L2Q","sr","h_zp_jigm",1)
    print('zz up: ',hzzup.Integral())
    hzzdwn = bkgs.getAddedHistXSErr(empty7,"ZZTo2L2Q","sr","h_zp_jigm",-1)
    print('zz dn: ',hzzdwn.Integral())
    hwzup  = bkgs.getAddedHistXSErr(empty8,"WZTo2L2Q","sr","h_zp_jigm",1)
    hwzdwn = bkgs.getAddedHistXSErr(empty9,"WZTo2L2Q","sr","h_zp_jigm",-1)
    print("wz up: ",hwzup.Integral())
    print("wz dn: ",hwzdwn.Integral())

    hvvup  = hzzup.Clone()
    hvvup.Add(hwzup)
    hvvdwn = hzzdwn.Clone()
    hvvdwn.Add(hwzdwn)
    hdyup  = dyEstup.Get("extrphistnoerrs").Clone()
    hdydwn = dyEstdwn.Get("extrphistnoerrs").Clone()

    print("vv up: ",hvvup.Integral())
    print("vv dn: ",hvvdwn.Integral())


    
    #Rename and Restructure
    httup = newNameAndStructure(httup,"TT_"+syst+"Up",rebindiv,limrangelow,limrangehigh)
    httdwn = newNameAndStructure(httdwn,"TT_"+syst+"Down",rebindiv,limrangelow,limrangehigh)
    hvvup = newNameAndStructure(hvvup,"VV_"+syst+"Up",rebindiv,limrangelow,limrangehigh)
    hvvdwn = newNameAndStructure(hvvdwn,"VV_"+syst+"Down",rebindiv,limrangelow,limrangehigh)
    hdyup = newNameAndStructure(hdyup,"DY_"+syst+"Up",1,limrangelow,limrangehigh)
    hdydwn = newNameAndStructure(hdydwn,"DY_"+syst+"Down",1,limrangelow,limrangehigh)
    httup = httup.Rebin(len(newbinedges)-1,"TT_"+syst+"Up",newbinedges)
    httdwn = httdwn.Rebin(len(newbinedges)-1,"TT_"+syst+"Down",newbinedges)
    hvvup = hvvup.Rebin(len(newbinedges)-1,"VV_"+syst+"Up",newbinedges)
    hvvdwn = hvvdwn.Rebin(len(newbinedges)-1,"VV_"+syst+"Down",newbinedges)
    hdyup = hdyup.Rebin(len(newbinedges)-1,"DY_"+syst+"Up",newbinedges)
    hdydwn = hdydwn.Rebin(len(newbinedges)-1,"DY_"+syst+"Down",newbinedges)

    print("ratio of zz up to zz nominal : ",hzzup.Integral()/hzz.Integral())
    print("ratio of zz dwn to zz nominal: ",hzzdwn.Integral()/hzz.Integral())
    print("ratio of wz up to zz nominal : ",hwzup.Integral()/hwz.Integral())
    print("ratio of zz dwn to zz nominal: ",hwzdwn.Integral()/hwz.Integral())
    
    #Debug for background
    #Extra plots
    #dynormup = np.load('analysis_output_ZpAnomalon/2023-05-08/Run2_161718_dynormalization_alphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_vvup_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
    #dynormnom = np.load(config.get('nominal','pathnom')+'/Run2_161718_dynormalization_alphat_'+config.get('nominal','strnom')+'_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
    #dynormdwn = np.load('analysis_output_ZpAnomalon/2023-05-08/Run2_161718_dynormalization_alphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_vvdwn_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]

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
    
    htt.GetYaxis().SetRangeUser(0,getGoodPlotRange([httup,htt,httdwn]))
    hdy.GetYaxis().SetRangeUser(0,getGoodPlotRange([hdyup,hdy,hdydwn]))
    hvv.GetYaxis().SetRangeUser(0,getGoodPlotRange([hvvup,hvv,hvvdwn]))
    
    dyupx = hdy.GetXaxis()
    dyupx.SetTitle("RJR Z' Jigsaw Mass Estimator (GeV)")
    ttupx = htt.GetXaxis()
    ttupx.SetTitle("RJR Z' Jigsaw Mass Estimator (GeV)")
    vvupx = hvv.GetXaxis()
    vvupx.SetTitle(" RJR Z' Jigsaw Mass Estimator (GeV)")
    hdy.GetYaxis().SetTitle("Events")
    hvv.GetYaxis().SetTitle("Events")
    hdy.SetStats(0)
    htt.SetStats(0)
    hvv.SetStats(0)
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
        
    #Labels
    dydifflab = makeDifferenceTPave(hdyup,hdy,hdydwn,"DY")
    ttdifflab = makeDifferenceTPave(httup,htt,httdwn,"TT")
    vvdifflab = makeDifferenceTPave(hvvup,hvv,hvvdwn,"VV")
    
    dyunclab = makeRateUncertaintyTPave(dyquad,dymid)
    ttunclab = makeRateUncertaintyTPave(ttquad,ttmid)
    vvunclab = makeRateUncertaintyTPave(vvquad,vvmid)
    
    dyjiglab = ROOT.TPaveText(.6,.4,.9,.55,"NBNDC")
    dyjiglab.AddText(dystrup)
    dyjiglab.AddText(dystrnom)
    dyjiglab.AddText(dystrdwn)
    dyjiglab.SetFillColor(0)
    dyjiglab.GetLine(0).SetTextColor(ROOT.kRed)
    dyjiglab.GetLine(1).SetTextColor(ROOT.kBlack)
    dyjiglab.GetLine(2).SetTextColor(ROOT.kBlue)
    
    ttjiglab = ROOT.TPaveText(.6,.3,.9,.4,"NBNDC")
    ttjiglab.AddText(ttstrup)
    ttjiglab.AddText(ttstrnom)
    ttjiglab.AddText(ttstrdwn)
    ttjiglab.SetFillColor(0)
    ttjiglab.GetLine(0).SetTextColor(ROOT.kRed)
    ttjiglab.GetLine(1).SetTextColor(ROOT.kBlack)
    ttjiglab.GetLine(2).SetTextColor(ROOT.kBlue)
    
    vvjiglab = ROOT.TPaveText(.6,.4,.9,.55,"NBNDC")
    vvjiglab.AddText(vvstrup)
    vvjiglab.AddText(vvstrnom)
    vvjiglab.AddText(vvstrdwn)
    vvjiglab.SetFillColor(0)
    vvjiglab.GetLine(0).SetTextColor(ROOT.kRed)
    vvjiglab.GetLine(1).SetTextColor(ROOT.kBlack)
    vvjiglab.GetLine(2).SetTextColor(ROOT.kBlue)
   
    ratsup,ratsdn = makeRatioPlots([hdyup,httup,hvvup],[hdydwn,httdwn,hvvdwn],[hdy,htt,hvv])
        
    #TCanvases
    tcbkg  = ROOT.TCanvas("tcbkg","backgrounds",1200,500)
    pdy    = ROOT.TPad("pdy","dysr",0,0.2,0.5,1.0)
    pdydiv = ROOT.TPad("pdydiv","dysrdiv",0,0.0,0.5,0.2)
    pvv    = ROOT.TPad("pvv","vvsr",0.5,0.2,1.0,1.0)
    pvvdiv = ROOT.TPad("pvvdiv","vvsrdiv",0.5,0.0,1.0,0.2)
    ptdy   = ROOT.TPad("ptdy","dysr",0,0.0,0.5,0.5)
    ptvv   = ROOT.TPad("ptvv","vvsr",0.5,0.0,1.0,0.5)

    CMS_lumi.extraText = "Z+jets estimation, internal"

    tcbkg.cd()
    pdy.Draw()
    pdy.cd()
    hdy.Draw('hist')
    hdyup.Draw('histsame')
    hdydwn.Draw('histsame')
    dyjiglab.Draw()
    dydifflab.Draw()
    #dyunclab.Draw()
    CMS_lumi.CMS_lumi(pdy,4,13)
    pdy.Update()
    tcbkg.cd()
    pdydiv.Draw()
    pdydiv.cd()
    ratsup[0].Draw()
    ratsdn[0].Draw("same")
    tcbkg.cd()
    #ptt.Draw()
    #ptt.cd()
    #htt.Draw('hist')
    #httup.Draw('histsame')
    #httdwn.Draw('histsame')
    #ttjiglab.Draw()
    #ttdifflab.Draw()
    #ttunclab.Draw()
    #ptt.Update()
    tcbkg.cd()
    CMS_lumi.extraText = "VV estimation, internal"
    pvv.Draw()
    pvv.cd()
    hvv.Draw('hist')
    hvvup.Draw('histsame')
    hvvdwn.Draw('histsame')
    vvjiglab.Draw()
    vvdifflab.Draw()
    #vvunclab.Draw()
    CMS_lumi.CMS_lumi(pvv,4,13)
    pvv.Update()
    tcbkg.cd()
    pvvdiv.Draw()
    pvvdiv.cd()
    ratsup[-1].Draw()
    ratsdn[-1].Draw("same")
    tcbkg.cd()
    tcbkg.Update()
        
    bkgcomp = go.makeOutFile('Run2_161718_ZllHbbMET',chan+'_bkg_'+syst,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tcbkg.SaveAs(bkgcomp)
