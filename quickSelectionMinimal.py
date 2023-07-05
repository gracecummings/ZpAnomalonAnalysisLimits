import uproot as up4
import uproot3 as up3
import pandas as pd
import numpy as np
import boost_histogram as bh
import argparse
import glob
import gecorg as go


parser = argparse.ArgumentParser()

def boostUnc(values,weights,nbins,binstart,binstop):
    boosth = bh.Histogram(bh.axis.Regular(bins=nbins,start=binstart,stop=binstop),storage=bh.storage.Weight())
    boosth.fill(values,weight=weights)
    boostvar = boosth.view().variance
    boosterr = np.sqrt(boostvar)
    return boosterr

def wrapPhi(phi):
    if phi < 0:
        wphi = -1*phi
    else:
        wphi = phi
    return wphi

def wrapDeltaPhi(dphi):
    if dphi > 3.14159:
        dp = 2*3.14159 - dphi
    else:
        dp = dphi
    return dp
 

def deltaPhi(v1phi,v2phi):
    dp = abs(v1phi-v2phi)
    wdp = dp.map(wrapDeltaPhi)
    #This returns a df with the same number of events as input
    #If a cut is introduced, need to carry weight column
    return wdp

def deltaR(v1phi,v2phi,v1eta,v2eta):
    dR = ((v2phi-v1phi)**2+(v2eta-v1eta)**2)**(1/2)
    return dR

if __name__=='__main__':
    parser.add_argument("-f","--sample",help = "sample file")
    parser.add_argument("-o","--output",help = "output file name")
    parser.add_argument("-b","--btagger",help = "btagger selection, need name part of key")
    parser.add_argument("-wp","--btagWP",type=float,default=0.6,help = "doulbeB tagger working point (default M1)")
    parser.add_argument("-zpt","--zPtCut",type=float,default = 100.0,help = "pT cut on Z")
    parser.add_argument("-hpt","--hPtCut",type=float,default = 250.0,help = "pT cut on h")
    parser.add_argument("-met","--metPtCut",type=float,default = 50.0,help = "pT cut on met")
    parser.add_argument("-sdm","--sdmCut",type=float,default = 10.0,help = "lowest soft drop mass cut")
    parser.add_argument("-date","--date",type=str,help = "where are your topiary plots?")
    parser.add_argument("-sr","--signalregion",type=bool,help = "do you want a signal region plot?")
    parser.add_argument("-c","--comboregion",type=bool,help = "do you want combined SR and SB?")
    parser.add_argument("-v","--validation",type=bool,help = "validation region bounds?")
    parser.add_argument("-syst","--systematics",type=str)
    args = parser.parse_args()

    samp   = args.sample
    sdmcut = args.sdmCut
    zptcut = args.zPtCut
    hptcut = args.hPtCut
    metcut = args.metPtCut
    btaggr = args.btagger
    btagwp = args.btagWP
    sr     = args.signalregion
    comb   = args.comboregion
    valid  = args.validation
    systl  = args.systematics

    #inputfiles = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/'+samp+'*_topiary*.root')
    #inputfiles = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/'+samp+'*_topiary*systnominal*.root')
    inputfiles = glob.glob(args.date+'/'+samp+'*_topiary*systnominal*.root')

    for fjec in inputfiles:
        #for jet systematics
        samp = fjec.split("/")[-1]
        front = samp.split("_Zpt")[0]
        jectype = front.split("_")[-1]
        samp = front.split("_topiary")[0]
        #inputfiles = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/'+samp+'*_topiary*'+jectype+'*.root')
        inputfiles = glob.glob(args.date+'/'+samp+'*_topiary*'+jectype+'*.root')
    
        print("Doing selections on:")
        print("    ",inputfiles[:1])
        print("    ",samp)
        print("    Concerning jet systematics:")
        print("    ",jectype)
        stype,year = go.sampleType(samp)

        if not valid:
            if sr and stype != 0:
                print("    using signal region selections")
            elif comb and stype != 0:
                print("    using full region selections")
            else:
                print("    using sideband selections")

        if valid:
            print("    Doing validation of alpha method cuts")
            if sr:
                print("    'signalr' labeled events are in soft drop mass bands (55,70]")
            else:
                print("    'sideband' labeled events are in soft drop mass bands [30,55],[150,5000]")

        metstr = ''
        branches = [b'MET',
                    b'METPhi',
                    b'METclean',
                    b'METPhiclean',
                    b'GenMET',
                    b'GenMETPhi',
                    b'event_weight_kf',
                    b'event_weight_btag'
        ]

    
        #events = up3.pandas.iterate(inputfiles[:1],'PreSelection;1',branches=branches)
        tree = up3.open(inputfiles[0])['PreSelection;1']
        events = tree.pandas.df(branches=branches)
        eventweights = events['event_weight_kf']*events['event_weight_btag']

        #lets make some histograms.
        rootfilename  = go.makeOutFile(samp,'upout_metcomp','.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))#need to update for btagger
        rootOutFile   = up3.recreate(rootfilename,compression = None)

        rootOutFile["h_pfmet_phi"]    = np.histogram(events['METPhi'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
        rootOutFile["h_pfmet"]        = np.histogram(events['MET'],bins=39,range=(50,2000),weights=eventweights)
        rootOutFile["h_metclean_phi"] = np.histogram(events['METPhiclean'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
        rootOutFile["h_metclean"]     = np.histogram(events['METclean'],bins=39,range=(50,2000),weights=eventweights)
        rootOutFile["h_genmet_phi"]   = np.histogram(events['GenMETPhi'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
        rootOutFile["h_genmet"]       = np.histogram(events['GenMET'],bins=39,range=(50,2000),weights=eventweights)
        
        #Book Keeping
        f = up3.open(inputfiles[0])
        rootOutFile["hnevents"]      = str(f['hnorigevnts'].values[0])
        rootOutFile["hnevents_passing"] = str(eventweights.sum())#str(len(btdf))
        print(eventweights.sum())

