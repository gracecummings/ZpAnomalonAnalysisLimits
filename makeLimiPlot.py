import ROOT
import glob
import os
import sys
import tdrstyle
import CMS_lumi
import gecorg_test as go
import numpy as np
import pandas as pd
import itertools
#import configparser

tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = "101.27 fb^{-1}"
CMS_lumi.writeExtraText = 0
#CMS_lumi.extraText = "Simulation Preliminary"


if __name__=='__main__':

    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'
    chan    = 'mumu'
    sigxs   = 1.0

    limitpath= "limholder/limWithNothing/h"
    #lims = glob.glob("limholder/limsWithJECs/h*")
    lims = glob.glob(limitpath+"*")
    zpbinwidth = 500
    ndbinwidth = 200

    #Creates list of lists with masses, names, and tfiles seperated
    #allows for pulling info based on a single particle mass
    limsinfo = [[x,int(x.split("/")[-1].split(".")[0].split("Combine")[-1].split("Zp")[-1].split("ND")[0]),int(x.split("/")[-1].split(".")[0].split("Combine")[-1].split("ND")[-1].split("NS")[0]),int(x.split("/")[-1].split(".")[0].split("Combine")[-1].split("NS")[-1])] for x in lims]
    limsdf = pd.DataFrame(limsinfo,columns = ['path','mzp','mnd','mns'])
    limsdf['path'].astype("string")

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
    coldivl = [x*0.02 for x in range(10,18)]
    coldivl.insert(0,0.0)
    coldivl = coldivl+[0.4,0.45,0.65,0.8,0.95,1.1]
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
    hlim.GetZaxis().SetTitle("Median cross section upper limit (95% CL)")
    hlim.GetZaxis().SetTitleSize(0.04)
    hlim.GetZaxis().SetTitleOffset(.9)
    hlim.GetZaxis().SetLabelSize(0.025)
    

    hlim.SetContour(len(coldiv),coldiv)
    
    for combout in lims:
        signame = combout.split("/")[-1].split(".")[0].split("Combine")[-1]
        mzpstr = signame.split("Zp")[-1].split("ND")[0]
        mndstr = signame.split("ND")[-1].split("NS")[0]
        mnsstr = signame.split("NS")[-1]

        f = ROOT.TFile.Open(combout)
        tree = f.Get("limit")
        tree.GetEntry(2)#The 50% quantiles for median limit
        limit = tree.limit
        zpbin = int((int(mzpstr)-zpmin)/zpbinwidth+1)
        ndbin = int((int(mndstr)-ndmin)/ndbinwidth+1)
        hlim.SetBinContent(zpbin,ndbin,limit)

        #print(signame)
        #print("    The median limit is: ",limit)
        #print("    The quantile is: ",tree.quantileExpected)
        #print("    The Zpbin: ",hlim.GetXaxis().GetBinCenter(zpbin))
        #print("    The NDbin: ",hlim.GetYaxis().GetBinCenter(ndbin))
        #print("    The bin content is: ",hlim.GetBinContent(zpbin,ndbin))
    
    tc.cd()
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
    plotname = go.makeOutFile('Run2_2017_2018_ZllHbbMET','limits_'+str(limsdf['mns'][0]),'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(plotname)

    #Make brazilian flag plot
    #for each value of Z' mass
    #print(list(set(limsdf['mzp'])))
    for mzpt in list(set(limsdf['mzp'])):
        print("Making Brazlian Flag plot for mZp = "+str(mzpt))
        title = "limit per m_{Z'} = "+str(mzpt)
        df = limsdf[limsdf['mzp'] == mzpt]
        #ndmax = max(limsdf['mnd'])
        #ndmin = min(limsdf['mnd'])
        ndbins = len(df['mnd'])
        #ndbinvals = np.zeros(ndbins)
        ndbinvals = []
        limits = []
        #1sigs  = []
        #2sigs 

        #tg = ROOT.TGraphErrors(ndbins,
        #hbf = ROOT.TH1F("hbf",title,ndbins,float((ndmin-ndbinwidth/2)),float((ndmax+ndbinwidth/2)))
        #hbf.SetStats(0)
        #hbf.GetXaxis().SetTitle("m_{ND} (GeV)")
        #hbf.GetYaxis().SetTitle("median xs upper limit (95% CL) fb")
        #print(hbf)
        for i,mndt in enumerate(sorted(df['mnd'])):
            print("Checking mnd = "+str(mndt))
            gdf = df[df['mnd'] == mndt]
            fileweneed = glob.glob(limitpath+"iggsCombineZp"+str(mzpt)+"ND"+str(mndt)+"*")
            f = ROOT.TFile.Open(fileweneed[0])
            tree = f.Get("limit")
            tree.GetEntry(2)#The 50% quantiles for median limit
            limit = tree.limit
            limits.append(limit)
            ndbinvals.append(float(mndt))
            #print("  Limit = ",limit)

        limits = np.array(limits)
        ndbinvals = np.array(ndbinvals)
        print(ndbins)
        print(limits)
        print(ndbinvals)
        print(min(limits))
        print(max(limits))
        tg = ROOT.TGraphErrors(int(ndbins),ndbinvals,limits)
        #tg.SetMinimum(min(limits)-min(limits)*.2)
        #tg.SetMaximum(max(limits)+max(limits)*.2)
        tg.SetMinimum(0)
        tg.SetMaximum(4)
        tc1.cd()
        tg.Draw()
        plotname = go.makeOutFile('Run2_2017_2018_ZllHbbMET','limits_Zp'+str(mzpt)+"_"+str(limsdf['mns'][0]),'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        #plotname = go.makeOutFile('Run2_2017_2018_ZllHbbMET','limits_Zp'+str(mzpt)+"_NS"+str(limsdf['mns'][0]),'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        tc1.SaveAs(plotname)


