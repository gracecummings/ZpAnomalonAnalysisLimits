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
CMS_lumi.lumi_13TeV = "137 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"

def getDeviatedOverNominal(hist,histnom):
    interest = hist.Integral()
    intnom = histnom.Integral()
    devonom = interest/intnom
    return devonom

def makeDeviatedTPave(histup,hist,histdwn,name):
    updiv  = getDeviatedOverNominal(histup,hist)
    dwndiv = getDeviatedOverNominal(histdwn,hist)
    updiv = "up deviated over nominal: "+str(round(updiv,3))
    dwndiv = "dwn deviated over nominal: "+str(round(dwndiv,3))
    lab = ROOT.TPaveText(.3,.3,.9,.5,"NBNDC")
    lab.AddText(updiv)
    lab.AddText(dwndiv)
    lab.SetFillColor(0)
    return lab

def getDeviatedOverNominalSummary(histup,hist,histdwn,name):
    upnum  = getDeviatedOverNominal(histup,hist)
    dwnnum = getDeviatedOverNominal(histdwn,hist)
    up = "Deviated over nominal  for "+name+" up variation: "+str(round(upnum,4))
    dwn = "Devoated over nominal for "+name+" dwn variation: "+str(round(dwnnum,4))
    return up,dwn

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

    pathbkg  = 'mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref'#This should be to mumu MC
    pathdata = 'emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref'#This should be to emu
    zptcut  = args.zptcut#'150.0'
    hptcut  = args.hptcut#'300.0'
    metcut  = args.metcut#'200.0'
    btagwp  = args.btagwp#'0.8'


    systr = 'systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco'

    emuscale = np.load(pathdata+'/Run2_2016_2017_2018_ttemunormalization_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
    emuunc   = np.load(pathdata+'/Run2_2016_2017_2018_ttemunormalization_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[1]
    jecup    = np.load('emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref_jec/Run2_2016_2017_2018_ttemunormalization_systjecup_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
    jecdn    = np.load('emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref_jec/Run2_2016_2017_2018_ttemunormalization_systjecdwn_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
    metdn    = np.load('emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref_unclmet/Run2_2016_2017_2018_ttemunormalization_systuncldwn_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
    metup    = np.load('emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref_unclmet/Run2_2016_2017_2018_ttemunormalization_systunclup_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
    btgup    = np.load('emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref_btag/Run2_2016_2017_2018_ttemunormalization_systnominal_hem_kf_btagup_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
    btgdn    = np.load('emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref_btag/Run2_2016_2017_2018_ttemunormalization_systnominal_hem_kf_btagdwn_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]

    print("Scaling the emu yield by ",emuscale)

    bkgs  = go.backgrounds(pathbkg,zptcut,hptcut,metcut,btagwp,systr)
    data  = go.run2(pathdata,zptcut,hptcut,metcut,btagwp,'systnominal_btagnom_muidnom')

    scales = {"nominal":emuscale,"Up":emuscale+emuunc,"Down":emuscale-emuunc}#"jecUp":jecup,"jecDown":jecdn,"unclmetUp":metup,"unclmetDown":metdn,"btagUp":btgup,"btagDown":btgdn}

    systscales = {"jec":[jecup,jecdn],"unclmet":[metup,metdn],"btag":[btgup,btgdn]}
    
    tf1 = ROOT.TFile(bkgs.a18ttsr[0])
    emptymzp = tf1.Get('h_zp_jigm').Clone()
    emptymzp.Reset("ICESM")
    emptysdm = tf1.Get('h_h_sd').Clone()
    emptysdm.Reset("ICESM")

    srplots = {}#key+scale then that holds
    trplots = {}#a list with mc,data,div

    #garwood poisson errors
    hsrdat1 = data.getAddedHist(emptymzp.Clone(),"sr","h_zp_jigm")#This is the emu channel,so we are good
    hsrdat1.Rebin(2)
    unchists = []
    for b in range(hsrdat1.GetNbinsX()+1):
        #if b == 0:
        #    continue
        hname = "TT_TT_StatsUncBin"+str(b)
        hup   = hsrdat1.Clone()
        hdn   = hsrdat1.Clone()
        c     = hsrdat1.GetBinContent(b)
        alpha = 1.- 0.682689492
        if c > 0:
            errup = ROOT.Math.gamma_quantile_c(alpha/2,int(c)+1,1)-c
            errdn = c - ROOT.Math.gamma_quantile_c(alpha/2,int(c),1)
            hup.SetBinContent(b,c+errup)
            hdn.SetBinContent(b,c+errdn)
            if (c - errdn) < 0:
                print("Uh oh! Unphysical Errors!")
        else:
            errup = ROOT.Math.gamma_quantile_c(alpha/2,int(c)+1,1)-c
            hup.SetBinContent(b,c+errup)
        hup.SetName("Up_"+hname)
        hdn.SetName("Down_"+hname)
        unchists.append(hup)
        unchists.append(hdn)
    uncout = go.makeOutFile('Run2_161718','emu_extrapStats_'+systr,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    uncoutf = ROOT.TFile(uncout,"recreate")
    for h in unchists:
        h.Scale(emuscale)
        h.Write()
    uncoutf.Close()

    
    for scale in scales.keys():
        emus = scales[scale]
        print('Checking the effects of ',scale)
        print(' using scale of ',emus)

        emptymzp1 = emptymzp.Clone()
        emptymzp2 = emptymzp.Clone()
        emptymzp3 = emptymzp.Clone()
        emptymzp4 = emptymzp.Clone()
        
        emptysdm1 = emptysdm.Clone()
        emptysdm2 = emptysdm.Clone()
        
        #SB would be the input to the alpha method. SR would be the SR estimation
        #hsbtt = bkgs.getAddedHist(emptymzp1,"TT","sb","h_zp_jigm")
        #hsbdat = data.getAddedHist(emptymzp2,"sb","h_zp_jigm")
        hsrtt = bkgs.getAddedHist(emptymzp3,"TT","sr","h_zp_jigm")
        hsrdat = data.getAddedHistPoissonErrors(emptymzp4,"sr","h_zp_jigm")#This is the emu channel,so we are good
        #hsbtt.Rebin(2)
        #hsbdat.Rebin(2)
        hsrtt.Rebin(2)
        hsrdat.Rebin(2)

        #TR is the dataset for the background template for the alpha method norm
        htrtt = bkgs.getAddedHist(emptysdm1,"TT","tr","h_h_sd")
        htrdat = data.getAddedHistPoissonErrors(emptysdm2,"tr","h_h_sd")

    
        #Scale'em
        #hsbdat.Scale(emuscale)
        hsrdat.Scale(emus)
        htrdat.Scale(emus)

        #for b in range(hsrdat.GetNbinsX()+1):
        #    print("Bin {0} content {1} error {2}".format(b,hsrdat.GetBinContent(b),hsrdat.GetBinErrorUp(b)))

        #make'em purty
        bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
        bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
        #setPlotParamsMC(hsbtt,"Z' Jigsaw Mass Estimator M_{Zp}",bkgcols[1])
        setPlotParamsMC(hsrtt,"Z' Jigsaw Mass Estimator M_{Zp} emuscale_"+scale,bkgcols[1])
        setPlotParamsMC(htrtt,"Higgs Candidate soft drop mass emuscale_"+scale,bkgcols[1])
        setPlotParamsData(hsrdat,"Z' Jigsaw Mass Estimator M_{Zp} emuscale_"+scale)
        #setPlotParamsData(hsbdat,"Z' Jigsaw Mass Estimator M_{Zp}")
        setPlotParamsData(htrdat,"Higgs Candidate soft drop mass emuscale_"+scale)
        
        #make ratios
        #hdivttsb = makeRatioHist(hsbdat,hsbtt,"emudata/mc")
        hdivttsr = makeRatioHist(hsrdat,hsrtt,"emudata/mc")
        hdivtttr = makeRatioHist(htrdat,htrtt,"emudata/mc")

        srplots[scale] = [hsrtt,hsrdat,hdivttsr]
        trplots[scale] = [htrtt,htrdat,hdivtttr]
        
    tc = ROOT.TCanvas("tc","tt channel overlays",1500,600)
    histpad1 = ROOT.TPad("histpad1","pad1",0,.2,.5,1)
    ratpad1 = ROOT.TPad("ratpad1","ratio1",0,0,.5,.2)
    histpad2 = ROOT.TPad("histpad2","pad2",0.5,.2,1,1)
    ratpad2  = ROOT.TPad("ratpad2","ratio2",.5,0,1,.2)

    ratlinemzp1 = ROOT.TLine(hdivttsr.GetBinLowEdge(1),1,hdivttsr.GetBinWidth(1)*hdivttsr.GetNbinsX(),1)
    ratlinemsd = ROOT.TLine(hdivtttr.GetBinLowEdge(1),1,hdivtttr.GetBinWidth(1)*hdivtttr.GetNbinsX(),1)
    ttleg  = ROOT.TLegend(0.55,0.60,0.9,0.8)
    #ttleg.AddEntry(hsbtt,"TT, $\mu\mu$","f" )
    #ttleg.AddEntry(hsbdat,"data, $e\mu$ scaled","ep" )
    ttleg.SetBorderSize(0)
    ttleg.AddEntry(trplots["nominal"][0],"TT $\mu\mu$ MC","f")
    ttleg.AddEntry(trplots["nominal"][1],"$e\mu$ data w/ scaling","pe")

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
    trplots["nominal"][0].Draw("HIST")
    CMS_lumi.CMS_lumi(histpad1,4,13)
    trplots["nominal"][1].Draw("SAME")
    ttleg.Draw()
    descriptxt0.Draw()
    histpad1.Update()
    tc.cd()
    ratpad1.Draw()
    ratpad1.cd()
    trplots["nominal"][2].Draw()
    ratlinemsd.Draw()
    tc.cd()
    histpad2.Draw()
    histpad2.cd()
    srplots["nominal"][0].Draw("HIST")
    CMS_lumi.CMS_lumi(histpad2,4,13)
    srplots["nominal"][1].Draw("SAME")
    ttleg.Draw()
    descriptxt2.Draw()
    histpad2.Update()
    tc.cd()
    ratpad2.Draw()
    ratpad2.cd()
    srplots["nominal"][2].Draw()
    ratlinemzp1.Draw()
    
    extrpplots = go.makeOutFile('Run2_161718','emu_scaled_set_on_mumu_mc_'+systr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(extrpplots)

    upyield,downyield = getDeviatedOverNominalSummary(srplots["Up"][1],srplots["nominal"][1],srplots["Down"][1],"emu")
    devp = makeDeviatedTPave(srplots["Up"][1],srplots["nominal"][1],srplots["Down"][1],"emu")

    print(upyield)
    print(downyield)

    tc1 = ROOT.TCanvas("tc1","tt up/dwn",750,600)

    uncsleg = ROOT.TLegend(0.55,0.6,0.9,0.8)
    uncsleg.SetBorderSize(0)

    tc1.cd()
    srplots["nominal"][1].Draw('hist')
    cols = {"nominal":ROOT.kBlack,"Up":ROOT.kRed,"Down":ROOT.kBlue}
    for scale in srplots.keys():
        h = srplots[scale][1]
        h.SetLineColor(cols[scale])
        uncsleg.AddEntry(h,"emu extrap scale "+scale,"l")
        h.Draw('histsame')
    uncsleg.Draw()
    devp.Draw()
    tc.cd()

    uncsplots = go.makeOutFile('Run2_161718','emu_scaled_uncertaintyiWithScale_'+systr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc1.SaveAs(uncsplots)

    #make extrapolation root file
    outf = go.makeOutFile('Run2_161718','emu_extrapolation_'+systr,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    outtf = ROOT.TFile(outf,"recreate")

    for key in srplots.keys():#this will be a dict with all the syst plots as well
        h = srplots[key][1]
        h.SetName("TT_scale"+key)
        h.Write()

    outtf.Close()
        
        
    for syst in systscales.keys():
        print("Checking the differences in ",syst)
        up = systscales[syst][0]
        dn = systscales[syst][1]
        emptymzp4 = emptymzp.Clone()
        print("      up: ",up)
        print("     nom: ",scales["nominal"])
        print("     dwn: ",dn)
        
        hsrdat = data.getAddedHist(emptymzp4,"sr","h_zp_jigm")#This is the emu channel,so we are good

        hsrdatnom = hsrdat.Clone()
        hsrdatup  = hsrdat.Clone()
        hsrdatdwn = hsrdat.Clone()

        hsrdatnom.Scale(scales["nominal"])
        hsrdatup.Scale(up)
        hsrdatdwn.Scale(dn)

        upy,dwny = getDeviatedOverNominalSummary(hsrdatup,hsrdatnom,hsrdatdwn,syst)
        print("    ",upy)
        print('    ',dwny)

