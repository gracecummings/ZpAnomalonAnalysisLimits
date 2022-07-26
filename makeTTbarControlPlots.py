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

def regionFormatter(regionstr):
    regdict = {'sideband':'sb','signalr':'sr','totalr':'tr'}
    formatted = regdict[regionstr]
    return formatted

def lumiFormatter(yearlist):
    lumidict = {16:35.9,17:41.53,18:59.74}
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
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-dir","--directory", type=str,help = "folder with input")
    parser.add_argument("-r","--region",help="region of phase space: totalr,sideband, or signalr")
    parser.add_argument("-y","--year", type=float,help = "year of samples eg. 2017 -> 17")
    parser.add_argument("-s","--syst",type=str,help="systematic string")
    parser.add_argument("-d","--data", type=bool,help = "plotting data?")
    args = parser.parse_args()

    #Get command line parameters
    zptcut    = args.zptcut
    hptcut    = args.hptcut
    metcut    = args.metcut
    btagwp    = args.btagwp
    year      = args.year
    regname   = args.region
    pathplots = args.directory
    systr     = args.syst
    plot_data = args.data

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
    #print(pathplots)
    ttmcemu  = go.backgrounds(pathplots,zptcut,hptcut,metcut,btagwp,systr)
    dataemu  = go.run2(pathplots,zptcut,hptcut,metcut,btagwp,'systnominal_btagnom_muidnom')

    #print(ttmcemu.bkgs)
    #print(dataemu.data[int(year)])
    
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
    testfile = ttmcemu.bkgs["TT"][testyear][reg][0][0]#you are making ttbar control plots
    testtfile = ROOT.TFile(testfile)
    keys = testtfile.GetListOfKeys()

    #names and param.
    #titles = {
    #    "h_z_pt":["Z p_{T} (GeV)",0,800,1],
    #    "h_z_eta":["\eta_{Z}",0,1200,1],
    #    "h_z_phi":["\phi_{Z}",0,700,2],
    #    "h_z_phiw":["\phi_{Z}",0,700,2],
    #    "h_z_m":["m_{Z} (GeV)",0,500,1],
    #    "h_h_pt":["Higgs p_{T} (GeV)",0,1100,1],
    #    "h_h_eta":["\eta_{Higss}",0,1100,1],
    #    "h_h_phi":["\phi_{Higgs}",0,600,2],
    #    "h_h_phiw":["\phi_{Higgs}",0,600,2],
    #    "h_h_m":["m_{h} (GeV)",0,450,1],
    #    "h_h_sd":["Higgs Soft Drop Mass (GeV)",0,350,1],#45 normally max
    #    "h_met":["p_{T}^{miss} (GeV)",0,2200,1],
    #    "h_met_phi":["\phi p_{T}^{miss}",0,750,2],
    #    "h_met_phiw":["\phi p_{T}^{miss}",0,750,2],
    #    "h_zp_jigm":["Jigsaw Mass Estimator Z'",0,1700,2],
    #    "h_nd_jigm":["Jigsaw Mass Estimator ND",0,1600,1],
    #    "h_ns_jigm":["Jigsaw Mass Estimator NS",0,3000,1],
    #    "h_btag":["btag operating point",0,2000,1],
    #    "h_dphi_zh":["\Delta\phi_{ZH}",0,900,2],
    #    "h_dphi_zmet":["\Delta\phi_{ZMET}",0,350,2],
    #    "h_dphi_hmet":["\Delta\phi_{HMET}",0,550,2],
    #    "h_dr_zh":["\Delta R(ZH)",0,1500,1],
    #    "h_dr_lmuh":["\Delta R(lmu,H)",0,800,1],
    #    "h_dr_slmuh":["\Delta R(slmu,H)",0,800,1],
    #    "h_dr_slmulmu":["\Delta R(slmu,lmu)",0,5200,1],
    #    "h_dr_gz_gh":["\Delta R(genZ,genH)",0,250,1],
    #    "h_dr_lmu_gh":["\Delta R(lmu,genH)",0,250,1],
    #    "h_dr_slmu_gh":["\Delta R(slmu,genH)",0,250,1],
    #    "h_LMu_pt":["leading \mu p_{T} (GeV)",0,5500,1],
    #    "h_LMu_phi":["\phi_{leading \mu}",0,5500,2],
    #    "h_LMu_eta":["\eta_{leading \mu }",0,5500,1],
    #    "h_sLMu_pt":["subleading \mu p_{T} (GeV)",0,5500,1],
    #    "h_sLMu_phi":["\phi_{subleading \mu}",0,5500,2],
    #    "h_sLMu_eta":["\eta_{subleading \mu}",0,5500,1],
    #    "h_dr_leadleph":["\Delta R(leading lepton,H)",0,800,1],
    #    "h_dr_sleadleph":["\Delta R(subleading lepton,H)",0,800,1],
    #    "h_dr_leps":["\Delta R(leptons)",0,800,1],
    #    "h_leadlep_pt":["leading lepton p_{T} (GeV)",0,5500,1],
    #    "h_sleadlep_pt":["subleading lepton p_{T} (GeV)",0,5500,1],
    #    "h_leadlep_phi":["\phi_{leading lepton}",0,5500,2],
    #    "h_sleadlep_phi":["\phi_{subleading lepton}",0,5500,2],
    #    "h_leadlep_eta":["\eta_{leading lepton}",0,5500,1],
    #    "h_sleadlep_eta":["\eta_{subleading lepton}",0,5500,1],
    #    "h_electron_pt":["electron p_{T} (GeV)",0,800,1],
    #    "h_electron_phi":["\phi_{electron}",0,700,2],
    #    "h_electron_eta":["\eta_{electron}",0,1200,1],
    #    "h_muon_pt":["muon p_{T} (GeV)",0,800,1],
    #    "h_muon_phi":["\phi_{muon}",0,700,2],
    #    "h_muon_eta":["\eta_{muon}",0,1200,1],
    #}

    titles = {
        "h_z_pt":["Z p_{T} (GeV)",0,30,2],
        "h_z_eta":["\eta_{Z}",0,40,1],
        "h_z_phi":["\phi_{Z}",0,40,2],
        "h_z_phiw":["\phi_{Z}",0,90,2],
        "h_z_m":["m_{Z} (GeV)",0,30,2],
        "h_h_pt":["Higgs p_{T} (GeV)",0,30,1],
        "h_h_eta":["\eta_{Higss}",0,40,1],
        "h_h_phi":["\phi_{Higgs}",0,40,2],
        "h_h_phiw":["\phi_{Higgs}",0,40,2],
        "h_h_m":["m_{h} (GeV)",0,30,1],
        "h_h_sd":["Higgs Soft Drop Mass (GeV)",0,30,1],#45 normally max
        "h_met":["p_{T}^{miss} (GeV)",0,50,1],
        "h_met_phi":["\phi p_{T}^{miss}",0,30,2],
        "h_met_phiw":["\phi p_{T}^{miss}",0,80,2],
        "h_zp_jigm":["Jigsaw Mass Estimator Z'",0,30,2],
        "h_nd_jigm":["Jigsaw Mass Estimator ND",0,60,1],
        "h_ns_jigm":["Jigsaw Mass Estimator NS",0,100,1],
        "h_btag":["btag operating point",0,40,1],
        "h_dphi_zh":["\Delta\phi_{ZH}",0,50,2],
        "h_dphi_zmet":["\Delta\phi_{ZMET}",0,40,2],
        "h_dphi_hmet":["\Delta\phi_{HMET}",0,40,2],
        "h_dr_zh":["\Delta R(ZH)",0,80,1],
        "h_dr_lmuh":["\Delta R(lmu,H)",0,200,1],
        "h_dr_slmuh":["\Delta R(slmu,H)",0,200,1],
        "h_dr_slmulmu":["\Delta R(slmu,lmu)",0,60,1],
        "h_dr_gz_gh":["\Delta R(genZ,genH)",0,250,1],
        "h_dr_lmu_gh":["\Delta R(lmu,genH)",0,250,1],
        "h_dr_slmu_gh":["\Delta R(slmu,genH)",0,250,1],
        "h_LMu_pt":["leading \mu p_{T} (GeV)",0,40,1],
        "h_LMu_phi":["\phi_{leading \mu",0,100,2],
        "h_LMu_eta":["\eta_{leading \mu ",0,100,1],
        "h_sLMu_pt":["subleading \mu p_{T} (GeV)",0,40,1],
        "h_sLMu_phi":["\phi_{subleading \mu}",0,100,2],
        "h_sLMu_eta":["\eta_{subleading \mu}",0,100,1],
        "h_electron_pt":["electron p_{T} (GeV)",0,20,1],
        "h_electron_phi":["\phi_{electron}",0,40,2],
        "h_electron_eta":["\eta_{electron}",0,60,1],
        "h_muon_pt":["muon p_{T} (GeV)",0,25,1],
        "h_muon_phi":["\phi_{muon}",0,40,2],
        "h_muon_eta":["\eta_{muon}",0,60,1],
        "h_dr_leps":["\Delta R(leptons)",0,60,1],
        "h_leadlep_pt":["leading lepton p_{T} (GeV)",0,20,1],
        "h_sleadlep_pt":["subleading lepton p_{T} (GeV)",0,30,1],
        "h_leadlep_phi":["\phi_{leading lepton}",0,60,2],
        "h_sleadlep_phi":["\phi_{subleading lepton}",0,60,2],
        "h_leadlep_eta":["\eta_{leading lepton}",0,60,1],
        "h_sleadlep_eta":["\eta_{subleading lepton}",0,60,1],
        "h_dr_leadleph":["\Delta R(leading lepton,H)",0,40,1],
        "h_dr_sleadleph":["\Delta R(subleading lepton,H)",0,40,1],
        "h_lmuon_pt":["leading muon p_{T} (GeV)",0,15,1],
        "h_lmuon_phi":["leading muom \phi",0,40,1],
        "h_lmuon_eta":["leading muon \phi",0,40,1],
        "h_slelectron_pt":["subleading electron p_{T} (GeV)",0,15,1],
        "h_slelectron_phi":["subleading electron \phi",0,40,1],
        "h_slelectron_eta":["subleading electron \eta",0,40,1],
        "h_lelectron_pt":["Leading electron p_{T} (GeV)",0,30,1],
        "h_lelectron_phi":["Leading electron \phi ",0,40,1],
        "h_lelectron_eta":["leading electron \eta ",0,40,1],
        "h_slmuon_pt":["subleading muon p_{T} ",0,30,1],
        "h_slmuon_phi":["subleading muon \phi ",0,40,1],
        "h_slmuon_eta":["subleading muon \eta ",0,40,1],
        "h_dr_leadmuonh":["\Delta R(leading muon,H)",0,40,1],
        "h_dr_subleadingeleh":["\Delta R(subleading electron,H)",0,40,1],
        "h_dr_leadeleh":["\Delta R(leading electron,H)",0,40,1],
        "h_dr_subleadingmuh":["\Delta R(subleading muon,H)",0,40,1],
        "h_metxy":["p_{T}^{miss} xy corrected (GeV)",0,100,1],
        "h_metxy_phi":["\phi p_{T}^{miss} xy corrected",0,80,2],
        "h_metxy_phiw":["\phi p_{T}^{miss} xy corrected",0,80,2],

    }



        #make the plots
    for key in keys:
        hname = key.GetName()
        print("Making control plot of ",hname)

        if ("gh" in hname) or ("gz" in hname) and plot_data:
            continue
        
        #Make holder histograms
        h = testtfile.Get(hname)
        if (not isinstance(h,ROOT.TH1)) or ('h_weights' in hname):
            continue
        empty = h.Clone()
        empty.Reset("ICESM")#creates an empty hist with same structure
        empty3 = empty.Clone()
        empty4 = empty.Clone()

        #Gather histograms
        htt  = ttmcemu.getAddedHist(empty3,"TT",reg,hname,years = years)
        htt.SetFillColor(bkgcols[1])
        htt.SetLineColor(bkgcols[1])
        htt.Rebin(titles[hname][3])
        print("Plot maximum: ",htt.GetMaximum())
        htt.GetYaxis().SetRangeUser(titles[hname][1],titles[hname][2]*1.5)
        print("Plot integral, mc: ",htt.Integral())


        #legend
        leg = ROOT.TLegend(0.5,0.7,0.88,0.80)
        leg.SetBorderSize(0)
        leg.AddEntry(htt,"TT, emu","f")

        #intrgral label
        yieldlabel = ROOT.TPaveText(0.5,0.6,0.88,0.7,"NBNDC")
        yieldlabel.SetFillColor(0)
        yieldlabel.AddText("Plot integral, mc: "+str(round(htt.Integral(),2)))

        if plot_data:
            #Make the data histograms
            hdat = dataemu.getAddedHist(empty4,reg,hname,years=years)
            #hdat = data.getAddedHist(empty6,'sb',hname,years=years)
            print("Plot integral, data: ",hdat.Integral())
            yieldlabel.AddText("Plot integral, data: "+str(round(hdat.Integral(),2)))
            hdat.Rebin(titles[hname][3])
            hdat.SetBinErrorOption(1)
            hdat.SetMarkerStyle(8)
            hdat.SetMarkerSize(0.7)
            hdat.SetMarkerColor(ROOT.kBlack)
            leg.AddEntry(hdat,"Data, emu")
            stkpadydims = [0.3,1.]
            ratpadydims = [0.0,0.3]
            tcanvasdims = [600,800]

            hdiv = hdat.Clone()
            hdiv.Divide(hdat,htt)

            ratline = ROOT.TLine(htt.GetBinLowEdge(1),1,htt.GetBinLowEdge(htt.GetNbinsX())+htt.GetBinWidth(htt.GetNbinsX()),1)

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
        
        htt.Draw("HIST")#add PFC for palette drawing
        htt.GetXaxis().SetTitle(titles[hname][0])
        htt.GetXaxis().SetTitleSize(0.05)
        htt.GetXaxis().SetTitleOffset(1.1)
        htt.GetXaxis().SetLabelSize(0.04)
        htt.GetYaxis().SetTitle("Events")
        htt.GetYaxis().SetTitleSize(0.05)
        htt.GetYaxis().SetLabelSize(0.04)
        CMS_lumi.CMS_lumi(p1,4,13)

        leg.Draw()
        yieldlabel.Draw()
        p1.Update()
        tc.cd()

        if plot_data:
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
        pngname = go.makeOutFile(hname,'emuratio_'+yearstr+'_'+regname,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        tc.SaveAs(pngname)

