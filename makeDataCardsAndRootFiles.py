import ROOT
import glob
import os
import sys
import gecorg as go
import numpy as np
import pandas as pd
import configparser

if __name__=='__main__':

    #will replace with command line options
    zptcut  = '150.0'
    hptcut  = '300.0'
    metcut  = '200.0'
    btagwp  = '0.8'

    bkgs = go.backgrounds('BkgInputs/',zptcut,hptcut,metcut,btagwp)
    sig  = go.signal('BkgInputs/',zptcut,hptcut,metcut,btagwp,10,101.27)#path tbd

    #Get ttbar, VV straight from selections
    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty2 = empty.Clone()
    empty3 = empty.Clone()

    htt = bkgs.getAddedHist(empty,"TT","sr","h_zp_jigm")
    hzz  = bkgs.getAddedHist(empty2,"ZZTo2L2Q","sr","h_zp_jigm")
    hwz  = bkgs.getAddedHist(empty3,"WZTo2L2Q","sr","h_zp_jigm")
    hvv  = hzz.Clone()
    hvv.Add(hwz)

    htt.SetName("ttbar")
    hvv.SetName("VV")
    
    prepName = go.makeOutFile('TestCombineInput','gatheredinfo','.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    prepFile = ROOT.TFile(prepName,"recreate")
    htt.Write()
    hvv.Write()

    
    
