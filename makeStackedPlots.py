import argparse
import ROOT
import glob
import os
import gecorg_test as go
import pandas as pd
import numpy as np
from datetime import date
from ROOT import kOrange, kViolet, kCyan, kGreen, kPink, kAzure, kMagenta, kBlue, kBird
from math import sqrt
import tdrstyle
import CMS_lumi

#tdrstyle.setTDRStyle()
#CMS_lumi.lumi_13TeV = "101.27 fb^{-1}"
#CMS_lumi.writeExtraText = 1
#CMS_lumi.extraText = "Simulation Preliminary"

def regionFormatter(regionstr):
    regdict = {'sideband':'sb','signalr':'sr','totalr':'tr'}
    formatted = regdict[regionstr]
    return formatted

def lumiFormatter(yearlist):
    lumidict = {16:36.31,17:41.53,18:59.74}
    lumi = 0
    for year in yearlist:
        lumi += lumidict[year]

    lumi = round(lumi,2)
    lumistr = str(lumi)+" fb^{-1}"
    return lumistr

def yearFormatter(yearlist):
    yearstr =''
    for year in yearlist:
        yearstr = yearstr+str(year)
    return yearstr
    


if __name__=='__main__':
    #build module objects
    parser = argparse.ArgumentParser()

    #Define parser imputs
    parser.add_argument("-x","--xsec", type=float,help = "desired siganl cross section in fb")
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-dir","--directory", type=str,help = "date folder with output")
    parser.add_argument("-r","--region",help="region of phase space: totalr,sideband, or signalr")
    parser.add_argument("-y","--year", type=float,help = "year of samples eg. 2017 -> 17")
    parser.add_argument("-s","--syst",type=str,help="systematic string")
    parser.add_argument("-d","--data", type=bool,help = "plotting data?")
    args = parser.parse_args()

    #Get command line parameters
    sig_xsec  = args.xsec
    zptcut    = args.zptcut
    hptcut    = args.hptcut
    metcut    = args.metcut
    btagwp    = args.btagwp
    year      = args.year
    regname   = args.region
    lumi      = 0
    pathplots = args.directory
    systr     = args.syst
    plot_data = args.data
    sigdivsor = 5
    valid = False
    islog = False
    
    #Select Plotting years and region
    years = [16,17,18]
    if year:
        years = [int(year)]
    yearstr = yearFormatter(years)
    reg = regionFormatter(regname)
    if plot_data:
        print('Plotting the data!')
    else:
        print('You know you are not plotting data, right?')
        
    #Gather Input
    bkgs  = go.backgrounds(pathplots,zptcut,hptcut,metcut,btagwp,systr)
    data  = go.run2(pathplots,zptcut,hptcut,metcut,btagwp,'systnominal_btagnom_muidnom')
    sigs =  go.signal(pathplots,zptcut,hptcut,metcut,btagwp,sig_xsec,years,systr)
    datareg = 'sb'
    if valid:
        #dynorm = np.load('analysis_output_ZpAnomalon/2022-05-17/Run2_161718_dynormalization_systnominal_kfnom_btagnom_muidnom_elidnom_elreconom_signalblind_Zptcut100.0_Hptcut300.0_metcut0.0_btagwp0.8.npy')[0]
        bkgs = go.validation('mumu_2022-03-31_ProperREOIDSF_validationOfAlphaMethodExtrap',zptcut,hptcut,metcut,btagwp,'systnominal_kfnom_btagnom_muidnom_elidnom_elreconom')
        data = go.validation('mumu_2022-03-31_ProperREOIDSF_validationOfAlphaMethodExtrap',zptcut,hptcut,metcut,btagwp,systr)
        datareg = 'vr'

    dynorm = 1
    if len(years) == 2:#dynorms only matter for composite years
        dynorm = np.load(pathplots+'/Run2_2017_2018_dynormalization_'+systr+'_signalblind_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.npy')[0]
    elif len(years) == 3:
            dynorm = np.load(pathplots+'/Run2_161718_dynormalization_alphat_'+systr+'_signalblind_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.npy')[0]

    #Colors, Naming, general style
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = lumiFormatter(years)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Simulation Preliminary"
    stkpadydims = [0.0,1.]
    ratpadydims = [0.0,0.0]
    tcanvasdims = [600,560]

    #Gather plots
    testyear = years[0]#picks first year in list, so desired year if only one
    testfile = bkgs.bkgs["DYJetsToLL"][testyear][reg][0][0]#stacked plots should always have DY
    testtfile = ROOT.TFile(testfile)
    keys = testtfile.GetListOfKeys()
    siginfo = sigs.getPreppedSig(reg,sig_xsec,years)
    sigcolors = go.colsFromPalette(siginfo,ROOT.kCMYK)
    siginfo = sorted(siginfo,key = lambda sig: (sig["mzp"],sig["mnd"])) 

    #names and param. To Do: expand to include plot limits for linear scale
    #titles = {
    #    "h_z_pt":["Z p_{T} (GeV)",0,40,1],
    #    "h_z_eta":["\eta_{Z}",0,100,1],
    #    "h_z_phi":["\phi_{Z}",0,90,2],
    #    "h_z_phiw":["\phi_{Z}",0,90,2],
    #    "h_z_m":["m_{Z} (GeV)",0,60,1],
    #    "h_h_pt":["Higgs p_{T} (GeV)",0,60,1],
    #    "h_h_eta":["\eta_{Higss}",0,130,1],
    #    "h_h_phi":["\phi_{Higgs}",0,70,2],
    #    "h_h_phiw":["\phi_{Higgs}",0,70,2],
    #    "h_h_m":["m_{h} (GeV)",0,50,1],
    #    "h_h_sd":["Higgs Soft Drop Mass (GeV)",0,75,1],#45 normally max
    #    "h_met":["p_{T}^{miss} (GeV)",0,100,1],
    #    "h_met_phi":["\phi p_{T}^{miss}",0,80,2],
    #    "h_met_phiw":["\phi p_{T}^{miss}",0,80,2],
    #    "h_zp_jigm":["Jigsaw Mass Estimator Z'",0,60,2],
    #    "h_nd_jigm":["Jigsaw Mass Estimator ND",0,60,1],
    #    "h_ns_jigm":["Jigsaw Mass Estimator NS",0,100,1],
    #    "h_btag":["btag operating point",0,70,1],
    #    "h_dphi_zh":["\Delta\phi_{ZH}",0,80,2],
    #    "h_dphi_zmet":["\Delta\phi_{ZMET}",0,60,2],
    #    "h_dphi_hmet":["\Delta\phi_{HMET}",0,60,2],
    #    "h_dr_zh":["\Delta R(ZH)",0,170,1],
    #    "h_dr_lmuh":["\Delta R(lmu,H)",0,200,1],
    #    "h_dr_slmuh":["\Delta R(slmu,H)",0,200,1],
    #    "h_dr_slmulmu":["\Delta R(slmu,lmu)",0,60,1],
    #    "h_dr_gz_gh":["\Delta R(genZ,genH)",0,250,1],
    #    "h_dr_lmu_gh":["\Delta R(lmu,genH)",0,250,1],
    #    "h_dr_slmu_gh":["\Delta R(slmu,genH)",0,250,1],
    #    "h_LMu_pt":["leading \mu p_{T} (GeV)",0,40,1],
    #    "h_LMu_phi":["\phi_{leading \mu}",0,100,2],
    #    "h_LMu_eta":["\eta_{leading \mu} ",0,100,1],
    #    "h_sLMu_pt":["subleading \mu p_{T} (GeV)",0,40,1],
    #    "h_sLMu_phi":["\phi_{subleading \mu}",0,100,2],
    #    "h_sLMu_eta":["\eta_{subleading \mu}",0,100,1],
    #}

    #For antibtagged alpha method region
    titles = {
        "h_z_pt":["Z p_{T} (GeV)",0,40,2],
        "h_z_eta":["\eta_{Z}",0,100,1],
        "h_z_phi":["\phi_{Z}",0,90,2],
        "h_z_phiw":["\phi_{Z}",0,90,2],
        "h_z_m":["m_{Z} (GeV)",0,100,1],
        "h_h_pt":["Higgs p_{T} (GeV)",0,60,1],
        "h_h_eta":["\eta_{Higss}",0,130,1],
        "h_h_phi":["\phi_{Higgs}",0,70,2],
        "h_h_phiw":["\phi_{Higgs}",0,70,2],
        "h_h_m":["m_{h} (GeV)",0,50,1],
        "h_h_sd":["Higgs Soft Drop Mass (GeV)",0,70,1],#45 normally max
        "h_metxy":["p_{T}^{miss} xy corrected (GeV)",0,100,1],
        "h_metxy_phi":["\phi p_{T}^{miss} xy corrected",0,80,2],
        "h_metxy_phiw":["\phi p_{T}^{miss} xy corrected",0,80,2],
        "h_met":["p_{T}^{miss} (GeV)",0,100,1],
        "h_met_phi":["\phi p_{T}^{miss}",0,80,2],
        "h_met_phiw":["\phi p_{T}^{miss}",0,80,2],
        "h_zp_jigm":["Jigsaw Mass Estimator Z'",0,60,2],
        "h_nd_jigm":["Jigsaw Mass Estimator ND",0,60,1],
        "h_ns_jigm":["Jigsaw Mass Estimator NS",0,100,1],
        "h_btag":["btag operating point",0,70,1],
        "h_dphi_zh":["\Delta\phi_{ZH}",0,80,2],
        "h_dphi_zmet":["\Delta\phi_{ZMET}",0,60,2],
        "h_dphi_hmet":["\Delta\phi_{HMET}",0,60,2],
        "h_dr_zh":["\Delta R(ZH)",0,170,1],
        "h_dr_lmuh":["\Delta R(lmu,H)",0,200,1],
        "h_dr_slmuh":["\Delta R(slmu,H)",0,200,1],
        "h_dr_slmulmu":["\Delta R(slmu,lmu)",0,60,1],
        "h_dr_gz_gh":["\Delta R(genZ,genH)",0,250,1],
        "h_dr_lmu_gh":["\Delta R(lmu,genH)",0,250,1],
        "h_dr_slmu_gh":["\Delta R(slmu,genH)",0,250,1],
        "h_LMu_pt":["leading \mu p_{T} (GeV)",0,40,1],
        "h_LMu_phi":["\phi_{leading \mu}",0,100,2],
        "h_LMu_eta":["\eta_{leading \mu} ",0,100,1],
        "h_sLMu_pt":["subleading \mu p_{T} (GeV)",0,40,1],
        "h_sLMu_phi":["\phi_{subleading \mu}",0,100,2],
        "h_sLMu_eta":["\eta_{subleading \mu}",0,100,1],
        "h_totmom":["Total Scalar Transverse Momentum",0,80,4],
        "h_metfrac":["MET/(Total Scalar Transverse Momentum)",0,100,2],
        "h_zptfrac":["p_{T}(Z)/(Total Scalar Transverse Momentum)",0,100,2],
        "h_hptfrac":["p_{T}(H)/(Total Scalar Transverse Momentum)",0,100,2],
        
    }


    #make the plots
    for key in keys:
        hname = key.GetName()
        print("working to make stacked plot of ",hname)
        #if "h_zp_jigm" not in hname:
        #    continue
        if ("gh" in hname) or ("gz" in hname) and plot_data:
            continue

        if ("h_h_sd" not in hname) and ("totalr" in regname):
            continue
        
        #Make holder histograms
        h = testtfile.Get(hname)
        if (not isinstance(h,ROOT.TH1)) or ('h_weights' in hname):
            continue
        empty = h.Clone()
        empty.Reset("ICESM")#creates an empty hist with same structure
        empty1 = empty.Clone()
        empty2 = empty.Clone()
        empty3 = empty.Clone()
        empty4 = empty.Clone()
        empty5 = empty.Clone()
        empty6 = empty.Clone()

        #Gather histograms
        hdy  = bkgs.getAddedHist(empty2,"DYJetsToLL",reg,hname,years = years)
        htt  = bkgs.getAddedHist(empty3,"TT",reg,hname,years = years)
        hzz  = bkgs.getAddedHist(empty4,"ZZTo2L2Q",reg,hname,years = years)
        hwz  = bkgs.getAddedHist(empty5,"WZTo2L2Q",reg,hname,years = years)

        #colors

        hdy.Scale(dynorm)
        hdy.SetFillColor(bkgcols[0])
        hdy.SetLineColor(bkgcols[0])
        hdy.Rebin(titles[hname][3])
        htt.SetFillColor(bkgcols[1])
        htt.SetLineColor(bkgcols[1])
        htt.Rebin(titles[hname][3])
        hwz.SetFillColor(bkgcols[2])
        hwz.SetLineColor(bkgcols[2])
        hwz.Rebin(titles[hname][3])
        hzz.SetFillColor(bkgcols[3])
        hzz.SetLineColor(bkgcols[3])
        hzz.Rebin(titles[hname][3])

        if ("gz" in hname) or ("gh" in hname):
            hdy.SetFillColor(0)
            hdy.SetLineColor(0)
            htt.SetFillColor(ROOT.kWhite)
            htt.SetLineColor(ROOT.kWhite)
            hwz.SetFillColor(ROOT.kWhite)
            hwz.SetLineColor(ROOT.kWhite)
            hzz.SetFillColor(ROOT.kWhite)
            hzz.SetLineColor(ROOT.kWhite)
        
        #Bkg Stack for Plotting
        hsbkg = ROOT.THStack("hsbkg","")
        hsbkg.Add(hzz)
        hsbkg.Add(hwz)
        hsbkg.Add(htt)
        hsbkg.Add(hdy)
        hsbkg.SetMinimum(titles[hname][1])
        hsbkg.SetMaximum(titles[hname][2])
        if islog:
            hsbkg.SetMinimum(0.01)
            hsbkg.SetMaximum(10000000)
        
        #Make added hist for ratio plotting
        hbkg = hzz.Clone()
        hbkg.Add(hwz)
        #hbkg.Add(hzz)
        hbkg.Add(htt)
        hbkg.Add(hdy)

        #legend
        leg = ROOT.TLegend(0.5,0.45,0.88,0.80)
        leg.SetBorderSize(0)
        leg.AddEntry(hdy,"DYJetsToLL","f")
        leg.AddEntry(htt,"TT","f")
        leg.AddEntry(hwz,"WZTo2L2Q","f")
        leg.AddEntry(hzz,"ZZTo2L2Q","f")
        
        #Data Histograms
        if plot_data:
            #Make the data histograms
            #hdat = data.getAddedHist(empty6,reg,hname,years=years)
            if valid:
                hdat = data.getAddedHistData(empty6,datareg,hname,years=years)
            if not valid:
                hdat = data.getAddedHist(empty6,datareg,hname,years=years)
            hdat.Rebin(titles[hname][3])
            hdat.SetBinErrorOption(1)
            hdat.SetMarkerStyle(8)
            hdat.SetMarkerSize(0.7)
            hdat.SetMarkerColor(ROOT.kBlack)
            leg.AddEntry(hdat,"Data")
            stkpadydims = [0.3,1.]
            ratpadydims = [0.0,0.3]
            tcanvasdims = [600,800]

            #Division for ratio plots
            hdiv = hdat.Clone()
            hdiv.Divide(hdat,hbkg)

            #ratline = ROOT.TLine(hbkg.GetBinLowEdge(1),1,hbkg.GetBinWidth(1)*titles[hname][3]*hbkg.GetNbinsX(),1)
            ratline = ROOT.TLine(hbkg.GetBinLowEdge(1),1,hbkg.GetBinLowEdge(hbkg.GetNbinsX())+hbkg.GetBinWidth(hbkg.GetNbinsX()),1)

        #Plotting itself

        tc = ROOT.TCanvas("tc",hname,tcanvasdims[0],tcanvasdims[1])
        p1 = ROOT.TPad("p1","stack_"+hname,0,stkpadydims[0],1.0,stkpadydims[1])
        p1.SetLeftMargin(0.15)
        p1.SetRightMargin(0.05)
        p2 = ROOT.TPad("p2","signif_"+hname,0,ratpadydims[0],1.0,ratpadydims[1])
        p2.SetRightMargin(.05)
        p2.SetLeftMargin(0.15)
        p2.SetBottomMargin(0.25)
        p2.SetTopMargin(0.05)

        #Prepare first pad for stack
        p1.Draw()
        p1.cd()
        #p1.SetLogy()

        #Draw the stack
        hsbkg.Draw("HIST")#add PFC for palette drawing
        hsbkg.GetXaxis().SetTitle(titles[hname][0])
        hsbkg.GetXaxis().SetTitleSize(0.05)
        hsbkg.GetXaxis().SetTitleOffset(1.1)
        hsbkg.GetXaxis().SetLabelSize(0.04)
        hsbkg.GetYaxis().SetTitle("Events")
        hsbkg.GetYaxis().SetTitleSize(0.05)
        hsbkg.GetYaxis().SetLabelSize(0.04)
        CMS_lumi.CMS_lumi(p1,4,13)


        #Print some parameters:
        print("The total background is: ",hbkg.Integral())
        
        #Draw the Signal
        for s,sig in enumerate(siginfo):
            print(sig["name"])
            if s % 5 != 0:
                continue
            if "Zp5500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8" not in sig["name"]:
                continue
            name = sig["name"]
            signame = 'holder'
            if "Tune" in name:
                strippedname = name.split("_Tune")[0]
                signame = strippedname.replace("-","")
            else:
                signame = name.replace("-","")

            hsig = sig["tfile"].Get(hname)
            hsig.Rebin(titles[hname][3])
            hsig.Scale(sig['scale'])
            hsig.SetLineColor(sigcolors[s])
            hsig.SetStats(0)
            print("The total signal for {0} is: {1}".format(signame,hsig.Integral()))


        
            if s == 0:
                hsig.Draw("histsame")
                leg.AddEntry(hsig,signame+" "+str(sig_xsec/1000)+" pb","l")
            elif s % sigdivsor == 0:
                hsig.Draw("histsame")
                leg.AddEntry(hsig,signame+" "+str(sig_xsec/1000)+" pb","l")
            #hsig.Draw("histsame")
            #leg.AddEntry(hsig,signame+" "+str(sig_xsec/1000)+" pb","l")
            p1.cd()

        leg.Draw()
        p1.Update()
        tc.cd()

        if plot_data:
            print("The total data is: ",hdat.Integral())
            p1.cd()
            hdat.Draw("histsame,pe")
            tc.cd()
            p2.Draw()
            p2.cd()
            hdiv.Draw()
            hdiv.GetXaxis().SetTitle("bin center")
            hdiv.GetXaxis().SetTitleSize(0.11)
            hdiv.GetXaxis().SetTitleOffset(0.65)
            hdiv.GetXaxis().SetLabelSize(0.075)
            hdiv.GetYaxis().SetTitle("data/MC")
            hdiv.GetYaxis().SetTitleSize(0.11)
            hdiv.GetYaxis().SetTitleOffset(.45)
            hdiv.GetYaxis().SetLabelSize(0.08)
            hdiv.GetYaxis().SetLabelOffset(0.02)
            hdiv.GetYaxis().SetNdivisions(503)
            hdiv.SetMinimum(0.)
            hdiv.SetMaximum(2.)
            ratline.Draw()

        tc.cd()
        
        #Save the plot
        pngname = go.makeOutFile(hname,'ratio_'+yearstr+'_'+regname,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        tc.SaveAs(pngname)


