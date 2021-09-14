import ROOT
import glob
import os
import sys
#import gecorg as go
import gecorg_test as go
import numpy as np
import pandas as pd
import itertools
#import configparser

if __name__=='__main__':

    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'
    chan    = 'mumu'
    sigxs   = 1.0

    lims = glob.glob("limholder/limsWithJECs/*")
    zpbinwidth = 500
    ndbinwidth = 200

    #Creates list of lists with masses, names, and tfiles seperated
    #allows for pulling info based on a single particle mass
    limsinfo = [[x,int(x.split("/")[-1].split(".")[0].split("Combine")[-1].split("Zp")[-1].split("ND")[0]),int(x.split("/")[-1].split(".")[0].split("Combine")[-1].split("ND")[-1].split("NS")[0]),int(x.split("/")[-1].split(".")[0].split("Combine")[-1].split("NS")[-1])] for x in lims]
    limsdf = pd.DataFrame(limsinfo,columns = ['path','mzp','mnd','mns'])

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

    hlim = ROOT.TH2F("hlim","",zpbins,float((zpmin-zpbinwidth/2)),float((zpmax+zpbinwidth/2)),ndbins,float((ndmin-ndbinwidth/2)),float((ndmax+ndbinwidth/2)))
    hlim.SetStats(0)
    
    tc = ROOT.TCanvas("tc","lims",600,600)

    hlim.Fill(1,1,1)
    hlim.Fill(1,2,2)
    print(hlim.GetBinContent(1,1))
    
    for combout in lims:
        signame = combout.split("/")[-1].split(".")[0].split("Combine")[-1]
        mzpstr = signame.split("Zp")[-1].split("ND")[0]
        mndstr = signame.split("ND")[-1].split("NS")[0]
        mnsstr = signame.split("NS")[-1]

        f = ROOT.TFile.Open(combout)
        tree = f.Get("limit")
        limit = tree.GetEntry(2)#The 50% quantiles for median limit
        zpbin = int((int(mzpstr)-zpmin)/zpbinwidth+1)
        ndbin = int((int(mndstr)-ndmin)/ndbinwidth+1)
        hlim.Fill(zpbin,ndbin,4)
        print(hlim.GetBinContent(zpbin,ndbin))

        print(signame)
        print("    The median limit is: ",tree.limit)
        print("    The quantile is: ",tree.quantileExpected)
        print("    The Zpbin: ",hlim.GetXaxis().GetBinCenter(zpbin))
        print("    The NDbin: ",hlim.GetYaxis().GetBinCenter(ndbin))

    
    tc.cd()
    hlim.Draw('colz')
    plotname = go.makeOutFile('Run2_2017_2018_ZllHbbMET','limits_'+str(limsdf['mns'][0]),'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(plotname)
