import ROOT
import glob
import os
import sys
import tdrstyle
import CMS_lumi
import gecorg_test as go
import numpy as np
import pandas as pd
import argparse

#import configparser

tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = go.lumiFormatter([16,17,18])
CMS_lumi.writeExtraText = 0
#


def makeBaseDataframe(filelist):
    #makes a basic data frame with dummy values for the limit and limit error. Columns for siganl mass values.
    #real values are written later
    #could probably use a multi index slicer (.xs) somewhere
    limsinfo = [[int(x.split("/")[-1].split(".")[0].split("Combine")[-1].split("_")[0].split("Zp")[-1].split("ND")[0]),int(x.split("/")[-1].split(".")[0].split("Combine")[-1].split("_")[0].split("ND")[-1].split("NS")[0]),int(x.split("/")[-1].split(".")[0].split("Combine")[-1].split("_")[0].split("NS")[-1]),-1.0,-1.0,-1.0,-1.0,-1.0] for x in filelist]
    limsdf = pd.DataFrame(limsinfo,columns = ['mzp','mnd','mns','limit','limitup','limitdn','limit2up','limit2dn'])
    return limsdf
    


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

    
    #zptcut  = '150.0'
    #hptcut  = '300.0'
    #metcut  = '200.0'
    #btagwp  = '0.8'
    chan    = 'mumu'
    sigxs   = 1.0

    #limitpath= "limholder/limWithNothing/h"
    #limitpath= "limholder/limsNoJECs/h"
    #limitpath= "limholder/limsNoJECBlownUpLumi/h"
    #limitpath = "limholder/limsWithJECs/h"
    #limitpath = "limholder/limsJEC100GeVBins/h"
    #limitpath = "limholder/pfMETFullGridStatsAndLumi/h"#"analysis_output_ZpAnomalon/2021-11-17/h"
    #limitpath = "limholder/pfMETFullGridSyst_NoAutoMCStats/h"
    limitpath = "analysis_output_ZpAnomalon/2022-10-01/higgs"
    lims = glob.glob(limitpath+"*"+"Zptcut"+zptcut+"_Hptcut"+hptcut+"_metcut"+metcut+"_btagwp"+btagwp+".txt*")
    interpolatedlims = glob.glob("analysis_output_ZpAnomalon/2022-10-05/higgs*")
    extrainterp      = []#glob.glob("analysis_output_ZpAnomalon/2022-10-05/interpcombine/higgs*")
    descrip ='nottathing' 
    zpbinwidth = 500
    ndbinwidth = 200

    #Creates a dataframe first with just the masspoints
    #wanted to incldue tfiles, but live and learn
    limsdf = makeBaseDataframe(lims+interpolatedlims+extrainterp)

    #Find relvant params
    zpmax = max(limsdf['mzp'])
    zpmin = min(limsdf['mzp'])
    ndmax = max(limsdf['mnd'])
    ndmin = min(limsdf['mnd'])
    numsig = len(limsdf)
    zpbins = int((zpmax-zpmin + zpbinwidth)/zpbinwidth)
    ndbins = int((ndmax-ndmin + ndbinwidth)/ndbinwidth)
    
    print("Number of signal samples: ",numsig)
    print("Zprime range [{0},{1}]".format(zpmin,zpmax))
    print("ND range     [{0},{1}]".format(ndmin,ndmax))
    print("NS considered : ",limsdf['mns'][0])
    print("Number of Z' bins: ",zpbins)
    print("Number of ND bins: ",ndbins)

    #Make the Canvas
    tc = ROOT.TCanvas("tc","lims",700,600)
    tc.SetLeftMargin(0.15)
    tc.SetRightMargin(0.15)

    tc1 = ROOT.TCanvas("tc1","brzllim",600,600)

    #Prep the z-axis colors
    coldivl = [x*0.005 for x in range(24,100)]
    coldivl.insert(0,0.0)
    coldivl+= [x*0.05 for x in range(10,20)]
    coldiv = np.array(coldivl)

    #Make the 2D histogram
    hlim = ROOT.TH2F("hlim","",zpbins,float((zpmin-zpbinwidth/2)),float((zpmax+zpbinwidth/2)),ndbins,float((ndmin-ndbinwidth/2)),float((ndmax+ndbinwidth/2)))
    hlim.SetStats(0)
    hlim.GetXaxis().SetTitle("m_{Z'} (GeV)")
    hlim.GetXaxis().SetTitleSize(0.05)
    hlim.GetXaxis().SetNdivisions(405)
    hlim.GetXaxis().SetLabelSize(0.04)
    hlim.GetYaxis().SetTitle("m_{ND} (GeV)")
    hlim.GetYaxis().SetTitleSize(0.05)
    hlim.GetYaxis().SetTitleOffset(1.35)
    hlim.GetYaxis().SetLabelSize(0.04)
    hlim.GetZaxis().SetTitle("Median cross section upper limit fb (95% CL)")
    hlim.GetZaxis().SetTitleSize(0.04)
    hlim.GetZaxis().SetTitleOffset(1.1)
    hlim.GetZaxis().SetLabelSize(0.025)
    #hlim.SetMinimum(0.1)
    hlim.SetContour(len(coldiv),coldiv)#for a custom color bar
    
    for combout in lims+interpolatedlims+extrainterp:
        signame = combout.split("/")[-1].split(".")[0].split("Combine")[-1].split("_")[0]
        #print(signame)
        mzpstr = signame.split("Zp")[-1].split("ND")[0]
        mndstr = signame.split("ND")[-1].split("NS")[0]
        mnsstr = signame.split("NS")[-1]

        f = ROOT.TFile.Open(combout)
        tree = f.Get("limit")
        tree.GetEntry(2)#The 50% quantiles for median limit
        limit = tree.limit
        tree.GetEntry(3)#The 84% quantile
        limitup = tree.limit
        tree.GetEntry(1)#The 16% quantile
        limitdn = tree.limit
        tree.GetEntry(0)#The 2.5% quantile
        limitdn2= tree.limit
        tree.GetEntry(4)#The 97% quantile
        limitup2= tree.limit

        zpbin = int((int(mzpstr)-zpmin)/zpbinwidth+1)
        ndbin = int((int(mndstr)-ndmin)/ndbinwidth+1)
        hlim.SetBinContent(zpbin,ndbin,limit)

        #want to add the limit to the df
        #the df and lims are in the same order, could use integers that way
        #but, if for some reason they are not, need to be mindful
        zpdf = limsdf[limsdf['mzp']==int(mzpstr)]
        ndzpdf = zpdf[zpdf['mnd']==int(mndstr)]
        indx = ndzpdf.index#returns the row that satisfies that slice
        limsdf.loc[indx,'limit'] = limit#edits the main dataframe
        limsdf.loc[indx,'limitup'] = limitup
        limsdf.loc[indx,'limitdn'] = limitdn
        limsdf.loc[indx,'limit2up'] = limitup2
        limsdf.loc[indx,'limit2dn'] = limitdn2

        #print(signame)
        #print("    The median limit is: ",limit)
        #print("     The upper limit is: ",limitup)
        #print("     The lower limit is: ",limitdn)
        
        #print("    The Zpbin: ",hlim.GetXaxis().GetBinCenter(zpbin))
        #print("    The NDbin: ",hlim.GetYaxis().GetBinCenter(ndbin))
        #print("    The bin content is: ",hlim.GetBinContent(zpbin,ndbin))

    tc.cd()
    #tc.SetLogz()
    hlim.Draw('colztext')#creates the palette objects
    tc.Update()#allows the palette object to be grabbed.
    colbar = ROOT.TPaletteAxis(hlim.GetListOfFunctions().FindObject("palette"))
    colbar.SetX1NDC(0.87)
    colbar.SetX2NDC(0.9)
    colbar.SetY1NDC(0.15)
    colbar.SetY2NDC(0.95)
    #unfortunately, does not seem to move the predrawn one
    hlim.Draw('coltext')#draw wihout
    colbar.Draw()
    tc.Update()
    CMS_lumi.CMS_lumi(tc,4,0)
    tc.Modified()
    tc.Update()
    
    tc.Update()
    plotname = go.makeOutFile('Run2_161718_ZllHbbMET','limits_grid_mns'+str(limsdf['mns'][0]),'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(plotname)

    #Make brazilian flag plot
    #for each value of Z' mass
    #print(list(set(limsdf['mzp'])))
    for mzpt in list(set(limsdf['mzp'])):
        print("Making Brazlian Flag plot for mZp = "+str(mzpt))
        title = "limit per m_{Z'} = "+str(mzpt)
        df = limsdf[limsdf['mzp'] == mzpt]
        ndbins = len(df['mnd'])



        sf = df.sort_values(by=['mnd'])
        ndbinvals = sf['mnd']
        limits = sf['limit']
        limsup = sf['limitup']
        limsdn = sf['limitdn']
        lims2up = sf['limit2up']
        lims2dn = sf['limit2dn']

        if mzpt == 1500:
            print(limits)

        limits = np.array(limits)
        limsup = np.array(limsup)-np.array(limits)
        limsdn = np.array(limits)-np.array(limsdn)
        lims2up = np.array(lims2up)-np.array(limits)
        lims2dn = np.array(limits)-np.array(lims2dn)
        xerrs   = np.zeros(int(ndbins))

        ndbinvals = np.array(ndbinvals,dtype=float)
        
        #tg = ROOT.TGraphErrors(int(ndbins),ndbinvals,limits)
        mg = ROOT.TMultiGraph()
        mg.SetTitle(";m_{ND} (GeV);#sigma B A (fb)")
        tg = ROOT.TGraphAsymmErrors(int(ndbins),ndbinvals,limits,xerrs,xerrs,limsdn,limsup)
        g = ROOT.TGraphErrors(int(ndbins),ndbinvals,limits)
        tg2 = ROOT.TGraphAsymmErrors(int(ndbins),ndbinvals,limits,xerrs,xerrs,lims2dn,lims2up)
        leg = ROOT.TLegend(0.60,0.55,0.93,0.8)
        empt = ROOT.TObject()
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        g.SetLineColor(ROOT.kBlack)
        g.SetLineWidth(3)
        g.SetLineStyle(9)

        tg.SetFillColor(ROOT.kGreen+1)
        tg.SetMarkerStyle(8)
        tg2.SetFillColor(ROOT.kOrange)
        tg2.SetMarkerStyle(8)

        leg.SetHeader("95% CL Upper Limits")
        leg.AddEntry(empt,"m_{Z'} = "+str(mzpt)+" GeV","")
        leg.AddEntry(g,"Median Expected","l")
        leg.AddEntry(tg,"68% expected","f")
        leg.AddEntry(tg2,"95% expected","f")
        mg.Add(tg2)
        mg.Add(tg)
        mg.Add(g)

        mg.SetMinimum(0.01)
        mg.SetMaximum(max(limits)*1000)
        #mg.GetXaxis().SetTitle("m_{ND} (GeV)")
        #mg.GetYaxis().SetTitle("\sigma (fb)")
                               
        
        tc1.cd()
        tc1.SetLogy()
        mg.Draw("AC3")
        leg.Draw()
        CMS_lumi.writeExtraText = 1
        CMS_lumi.extraText = "Simulation Preliminary"
        CMS_lumi.CMS_lumi(tc1,4,13)
        tc1.Modified()
        tc1.Update()

        plotname = go.makeOutFile('Run2_161718_ZllHbbMET','limits_Zp'+str(mzpt)+"_NS"+str(limsdf['mns'][0]),'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        tc1.SaveAs(plotname)
        tg.Delete()
        tc1.Clear()

    pklname = go.makeOutFile('Run2_161718_ZllHbbMET','limits_'+descrip+'_'+str(limsdf['mns'][0]),'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    limsdf.to_pickle(pklname)

