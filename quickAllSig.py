import ROOT
import glob
import os
import sys
#import gecorg as go
import gecorg_test as go
import numpy as np
import pandas as pd
import configparser


if __name__=='__main__':
    
    #will replace with command line options
    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'
    chan    = 'mumu'
    sigxs   = 1.0


    sig  = go.signal('SignalInputs/',zptcut,hptcut,metcut,btagwp,sigxs,101.27)
    siginfo = sig.prepsigsr
    sigcolors = go.colsFromPalette(siginfo,ROOT.kCMYK)
    tcAllSig = ROOT.TCanvas("tcAllSig","AllSig",600,600)
    sigAllname = go.makeOutFile("allsig_2018_NS200",'nominaljec','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))

    siginfo = sorted(siginfo,key = lambda sig: (sig["mzp"],sig["mnd"]))                    
    leg = ROOT.TLegend(0.50,0.50,0.88,0.88)
    for s,sig in enumerate(siginfo):
        name = sig["name"]
        print(name)
        signame = 'holder'
        if "Tune" in name:
            strippedname = name.split("_Tune")[0]
            signame = strippedname.replace("-","")
        else:
            signame = name.replace("-","")
        hsig = sig["tfile"].Get("h_zp_jigm")
        hsig.SetLineColor(sigcolors[s])
        hsig.SetMaximum(700)
        hsig.SetStats(0)

        tcAllSig.cd()
        if s == 0:
            hsig.Draw("hist")
            leg.AddEntry(hsig,signame,"l")
        elif s % 3 == 0:
            hsig.Draw("histsame")
            leg.AddEntry(hsig,signame,"l")
        tcAllSig.Update()

    leg.SetBorderSize(0)
    leg.Draw()
    tcAllSig.Update()
    tcAllSig.SaveAs(sigAllname)
