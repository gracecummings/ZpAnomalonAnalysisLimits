import ROOT
import glob
import os
import sys
import tdrstyle
import CMS_lumi
import gecorg_test as go
import numpy as np
import pandas as pd

if __name__=='__main__':
    zptcut  = '150-175'
    hptcut  = '300.0'
    metcut  = '50-50'
    btagwp  = '0.8'
    chan    = 'mumu'
    sigxs   = 1.0
    zpbinwidth = 500
    ndbinwidth = 200

    #pathbase = 'limholder/limsJEC100GeVBins/limitOutput100GeVBins/Run2_2017_2018_ZllHbbMET_limits_nottathing_200_Zptcut150.0_Hptcut300.0_metcut200.0_btagwp0.8.pkl'
    #pathcomp = 'limholder/limsJEC200GeVBins/limitanalysis200GeV/Run2_2017_2018_ZllHbbMET_limits_nottathing_200_Zptcut150.0_Hptcut300.0_metcut200.0_btagwp0.8.pkl'
    pathbase = 'analysis_output_ZpAnomalon/2021-11-17/Run2_2017_2018_ZllHbbMET_limits_nottathing_200_Zptcut150.0_Hptcut300.0_metcut50.0_btagwp0.8.pkl'
    pathcomp = 'analysis_output_ZpAnomalon/2021-11-17/Run2_2017_2018_ZllHbbMET_limits_nottathing_200_Zptcut175.0_Hptcut300.0_metcut50.0_btagwp0.8.pkl'
    #title = "Zpt > {0}, MET cut 75 lims - 50 lims (lumi+statsonly)".format(zptcut)
    title = '(Zpt > 175, MET > 50) - (Zpt > 150, MET > 50) % diff'

    df0 = pd.read_pickle(pathbase)
    df1 = pd.read_pickle(pathcomp)

    zpmax = max(df0['mzp'])
    zpmin = min(df0['mzp'])
    ndmax = max(df0['mnd'])
    ndmin = min(df0['mnd'])
    numsig = len(df0)
    zpbins = int((zpmax-zpmin + zpbinwidth)/zpbinwidth)
    ndbins = int((ndmax-ndmin + ndbinwidth)/ndbinwidth)

    tc = ROOT.TCanvas("tc","lims",700,600)
    tc.SetLeftMargin(0.15)
    tc.SetRightMargin(0.15)
    hlim = ROOT.TH2F("hlim",title,zpbins,float((zpmin-zpbinwidth/2)),float((zpmax+zpbinwidth/2)),ndbins,float((ndmin-ndbinwidth/2)),float((ndmax+ndbinwidth/2)))
    hlim.SetStats(0)
    hlim.GetXaxis().SetTitle("m_{Z'} (GeV)")
    hlim.GetXaxis().SetTitleSize(0.05)
    hlim.GetXaxis().SetTitleOffset(.8)
    hlim.GetXaxis().SetNdivisions(405)
    hlim.GetXaxis().SetLabelSize(0.04)
    hlim.GetYaxis().SetTitle("m_{ND} (GeV)")
    hlim.GetYaxis().SetTitleSize(0.05)
    hlim.GetYaxis().SetTitleOffset(1.35)
    hlim.GetYaxis().SetLabelSize(0.04)
    hlim.GetZaxis().SetTitle("percent difference between limits")
    hlim.GetZaxis().SetTitleSize(0.04)
    hlim.GetZaxis().SetTitleOffset(.9)
    hlim.GetZaxis().SetLabelSize(0.025)
    #hlim.SetContour(len(coldiv),coldiv)
    ROOT.gStyle.SetPalette(ROOT.kPastel)

    #print(df0)
    #print(df1)

    df0 = df0.sort_values(by=["mzp","mnd"])
    df1 = df1.sort_values(by=["mzp","mnd"])
    
    #print(df0)
    #print(df1)

    #print(list(df0['mzp']))
    #print(list(df1['mzp']))

    df0 = df0.drop_duplicates()
    df1 = df1.drop_duplicates()



    print(df0)
    print(df1)

    df0 = df0.reset_index()
    df1 = df1.reset_index()

    print(df0)
    print(df1)

    print(list(df0['mzp']))
    print(list(df1['mzp']))
    print(list(df0['mnd']))
    print(list(df1['mnd']))
    
    if (list(df0['mzp']) == list(df1['mzp']) and list(df0['mnd']) == list(df1['mnd'])):
        print('Order of the entires is OK, beginning comp')
        dfdiff = df1['limit'] - df0['limit']
        dfdiffperc = dfdiff/df0['limit']*100
        #print(len(df0))
        print(dfdiff)

        for i in range(len(df0)):
            mzp = df0.iloc[i]['mzp']
            mnd = df0.iloc[i]['mnd']
            percdiff = dfdiffperc.iloc[i]
            zpbin = (mzp-zpmin)/zpbinwidth+1
            ndbin = (mnd-ndmin)/ndbinwidth+1
            hlim.SetBinContent(int(zpbin),int(ndbin),percdiff)

        tc.cd()
        hlim.Draw("colztext")
        tc.Update()

        plotname = go.makeOutFile('Run2_2017_2018_ZllHbbMET','limit_comp_nothing_stdlumi','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        tc.SaveAs(plotname)

    else:
        print("order of entries not okay, uncomment debug statements")
