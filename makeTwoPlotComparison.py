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

def makeRatioPlotsSimple(ul,dl):
    rats = []
    for i,hist in enumerate(ul):
        nom = hist
        dn  = dl[i]
        rat = nom.Clone()
        rat.Divide(nom,dn)
        makeRatioReadable(rat,"May08/Jan23")
        rat.SetMarkerSize(1)
        rat.SetMarkerSize(.8)
        rats.append(rat)
    return(rats)


#def newNameAndStructure(hist,name,rebindiv,limrangelow,limrangehigh):
#    hist.Rebin(rebindiv)
#    nbins = hist.GetNbinsX()
#    binw  = hist.GetBinWidth(1)
#    newbins = [limrangelow+x*binw for x in range(int((limrangehigh-limrangelow)/binw))]
#   nh = ROOT.TH1F(name,name,len(newbins),limrangelow,limrangehigh)
#    for b,le in enumerate(newbins):
#        bnum = hist.FindBin(le)
#        bincontent = hist.GetBinContent(bnum)
#        binerror   = hist.GetBinError(bnum)
#        nh.SetBinContent(b+1,bincontent)
#        nh.SetBinError(b+1,binerror)
#    return nh

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

    dyEst1 = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-08/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    dyEst2 = ROOT.TFile('analysis_output_ZpAnomalon/2023-01-23/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')

    keys = dyEst1.GetListOfKeys()

    for key in keys:
       h1 = dyEst1.Get(key.GetName())
       h2 = dyEst2.Get(key.GetName())

       if "TH1D" not in str(type(h1)):
          #Strings to print
          #h1des  = "May 08, integral: "+str(round(h1.Integral(),3))
          #h2des  = "Jan 23, integral: "+str(round(h2.Integral(),3))
          #dyjiglab = ROOT.TPaveText(.6,.4,.9,.55,"NBNDC")
          #dyjiglab.AddText(h1des)
          #dyjiglab.AddText(h2des)
          #dyjiglab.SetFillColor(0)
          #dyjiglab.GetLine(0).SetTextColor(ROOT.kRed)
          #dyjiglab.GetLine(1).SetTextColor(ROOT.kBlue)
          h1.SetLineColor(ROOT.kRed)
          h2.SetLineColor(ROOT.kBlue)
          
          #TCanvases
          tcbkg  = ROOT.TCanvas("tcbkg","backgrounds",1200,500)
          pdy    = ROOT.TPad("pdy","dysr",0,0.2,0.5,1.0)
          pdydiv = ROOT.TPad("pdydiv","dysrdiv",0,0.0,0.5,0.2)
          pvv    = ROOT.TPad("pvv","vvsr",0.5,0.2,1.0,1.0)
          pvvdiv = ROOT.TPad("pvvdiv","vvsrdiv",0.5,0.0,1.0,0.2)
          ptdy   = ROOT.TPad("ptdy","dysr",0,0.0,0.5,0.5)
          ptvv   = ROOT.TPad("ptvv","vvsr",0.5,0.0,1.0,0.5)
          
          CMS_lumi.extraText = "internal"
          
          tcbkg.cd()
          pdy.Draw()
          pdy.cd()
          h1.Draw()
          h2.Draw('same')
          #dyjiglab.Draw()
          #dyunclab.Draw()
          CMS_lumi.CMS_lumi(pdy,4,13)
          pdy.Update()
          tcbkg.cd()
          pdydiv.Draw()
          pdydiv.cd()
          #rats[0].Draw()
          tcbkg.cd()
          tcbkg.cd()
          pvv.Draw()
          pvv.SetLogy()
          pvv.cd()
          h1.Draw()
          h2.Draw('same')
          #ttjiglab.Draw()
          #vvunclab.Draw()
          CMS_lumi.CMS_lumi(pvv,4,13)
          pvv.Update()
          tcbkg.cd()
          pvvdiv.Draw()
          pvvdiv.cd()
          #rats[-1].Draw()
          tcbkg.cd()
          tcbkg.Update()
        
          bkgcomp = go.makeOutFile('Run2_161718_ZllHbbMET',chan+'_bkg_dycompmayjanhell_'+key.GetName(),'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
          tcbkg.SaveAs(bkgcomp)


       else:
          h1rebin = h1.Clone()
          h2rebin = h2.Clone()

          h1rebin = go.newNameAndStructure(h1rebin,"May08Rebin_"+key.GetName(),1,limrangelow,limrangehigh)
          h2rebin = go.newNameAndStructure(h2rebin,"Jan23Rebin_"+key.GetName(),1,limrangelow,limrangehigh)
       
          newbinedges = np.array([1800.,2000.,2200.,2400.,2600.,2800.,10000.])
          
          h1rebin.Rebin(len(newbinedges)-1,"May08Rebin_"+key.GetName(),newbinedges)
          h2rebin.Rebin(len(newbinedges)-1,"Jan23Rebin_"+key.GetName(),newbinedges)
          
          #Plotting parameters
          h1.SetLineColor(ROOT.kRed)
          h1rebin.SetLineColor(ROOT.kRed)
          h2.SetLineColor(ROOT.kBlue)
          h2rebin.SetLineColor(ROOT.kBlue)
          
          h1.GetYaxis().SetRangeUser(0,getGoodPlotRange([h1,h2]))
          h2.GetYaxis().SetRangeUser(0,getGoodPlotRange([h1,h2]))
          h1rebin.GetYaxis().SetRangeUser(0,getGoodPlotRange([h1rebin,h2rebin]))
          h2rebin.GetYaxis().SetRangeUser(0,getGoodPlotRange([h1rebin,h2rebin]))
          
          dyupx = h1.GetXaxis()
          dyupx.SetTitle("RJR Z' Jigsaw Mass Estimator (GeV)")
          ttupx = h2.GetXaxis()
          ttupx.SetTitle("RJR Z' Jigsaw Mass Estimator (GeV)")
          vvupx = h1rebin.GetXaxis()
          vvupx.SetTitle(" RJR Z' Jigsaw Mass Estimator (GeV)")
          wwupx = h2rebin.GetXaxis()
          wwupx.SetTitle(" RJR Z' Jigsaw Mass Estimator (GeV)")
          h1.GetYaxis().SetTitle("Events")
          h1rebin.GetYaxis().SetTitle("Events")
          h2.GetYaxis().SetTitle("Events")
          h2rebin.GetYaxis().SetTitle("Events")
          h1.SetStats(0)
          h2.SetStats(0)
          h1rebin.SetStats(0)
          h2rebin.SetStats(0)
          
          #Strings to print
          h1des  = "May 08, integral: "+str(round(h1.Integral(),3))
          h2des  = "Jan 23, integral: "+str(round(h2.Integral(),3))
          h1desre  = "May 08 rebin, integral: "+str(round(h1rebin.Integral(),3))
          h2desre  = "Jan 23 rebin, integral: "+str(round(h2rebin.Integral(),3))
          
          dyjiglab = ROOT.TPaveText(.6,.4,.9,.55,"NBNDC")
          dyjiglab.AddText(h1des)
          dyjiglab.AddText(h2des)
          dyjiglab.SetFillColor(0)
          dyjiglab.GetLine(0).SetTextColor(ROOT.kRed)
          dyjiglab.GetLine(1).SetTextColor(ROOT.kBlue)
          
          ttjiglab = ROOT.TPaveText(.6,.3,.9,.4,"NBNDC")
          ttjiglab.AddText(h1desre)
          ttjiglab.AddText(h2desre)
          ttjiglab.SetFillColor(0)
          ttjiglab.GetLine(0).SetTextColor(ROOT.kRed)
          ttjiglab.GetLine(1).SetTextColor(ROOT.kBlue)
          
          rats = makeRatioPlotsSimple([h1,h1rebin],[h2,h2rebin])
          
          #TCanvases
          tcbkg  = ROOT.TCanvas("tcbkg","backgrounds",1200,500)
          pdy    = ROOT.TPad("pdy","dysr",0,0.2,0.5,1.0)
          pdydiv = ROOT.TPad("pdydiv","dysrdiv",0,0.0,0.5,0.2)
          pvv    = ROOT.TPad("pvv","vvsr",0.5,0.2,1.0,1.0)
          pvvdiv = ROOT.TPad("pvvdiv","vvsrdiv",0.5,0.0,1.0,0.2)
          ptdy   = ROOT.TPad("ptdy","dysr",0,0.0,0.5,0.5)
          ptvv   = ROOT.TPad("ptvv","vvsr",0.5,0.0,1.0,0.5)
          
          CMS_lumi.extraText = "internal"
          
          tcbkg.cd()
          pdy.Draw()
          pdy.cd()
          h1.Draw('hist')
          h2.Draw('histsame')
          dyjiglab.Draw()
          #dyunclab.Draw()
          CMS_lumi.CMS_lumi(pdy,4,13)
          pdy.Update()
          tcbkg.cd()
          pdydiv.Draw()
          pdydiv.cd()
          rats[0].Draw()
          tcbkg.cd()
          tcbkg.cd()
          pvv.Draw()
          pvv.cd()
          h1rebin.Draw('hist')
          h2rebin.Draw('histsame')
          ttjiglab.Draw()
          #vvunclab.Draw()
          CMS_lumi.CMS_lumi(pvv,4,13)
          pvv.Update()
          tcbkg.cd()
          pvvdiv.Draw()
          pvvdiv.cd()
          rats[-1].Draw()
          tcbkg.cd()
          tcbkg.Update()
          
          bkgcomp = go.makeOutFile('Run2_161718_ZllHbbMET',chan+'_bkg_dycompmayjanhell_'+key.GetName(),'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
          tcbkg.SaveAs(bkgcomp)
