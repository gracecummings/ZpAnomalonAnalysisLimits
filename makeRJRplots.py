import uproot as up4
import uproot3 as up3
import pandas as pd
import numpy as np
import argparse
import glob
import gecorg_test as go

if __name__=='__main__':
    parser.add_argument("-f","--sample",help = "sample file")
    parser.add_argument("-o","--output",help = "output file name")
    parser.add_argument("-b","--btagger",help = "btagger selection, need name part of key")
    parser.add_argument("-wp","--btagWP",type=float,default=0.6,help = "doulbeB tagger working point (default M1)")
    parser.add_argument("-zpt","--zPtCut",type=float,default = 100.0,help = "pT cut on Z")
    parser.add_argument("-hpt","--hPtCut",type=float,default = 250.0,help = "pT cut on h")
    parser.add_argument("-met","--metPtCut",type=float,default = 50.0,help = "pT cut on met")
    parser.add_argument("-sdm","--sdmCut",type=float,default = 10.0,help = "lowest soft drop mass cut")
    parser.add_argument("-d","--directory",type=str,help = "where are your topiary plots?")
    parser.add_argument("-sr","--signalregion",type=bool,help = "do you want a signal region plot?")
    parser.add_argument("-tot","--comboregion",type=bool,help = "do you want combined SR and SB?")
    parser.add_argument("-v","--validation",type=bool,help = "validation region bounds?")
    parser.add_argument("-syst","--systematics",type=str)
    parser.add_argument("-a","--alphar",type=bool,help = "alpha ratio regions?")
    parser.add_argument("-c","--chan",type=str)
    args = parser.parse_args()

    sampname   = args.sample
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
    alphatest = args.alphar
    channel = args.chan
    topdir = args.directory

    inputfiles = glob.glob(topdir+'/'+sampname+'*_topiary_'+channel+'*.root')
            branches = [b'ZCandidate_*',
                    b'RunNum',
                    b'LumiBlockNum',
                    b'EvtNum',
                    b'hCandidate_*',
                    b'metsuable',
                    b'metphiusable',
                    b'ZPrime_mass_est',
                    b'ND_mass_est',
                    b'NS_mass_est',
                    b'event_weight*',
        ]
        inputfiles = glob.glob(topdir+'/'+samp+'*_topiary_'+channel+'_*'+jectype+'*.root')
        print("Doing selections on:")
        print("    ",inputfiles[:1])
        print("    ",samp)
        #events = up3.pandas.iterate(inputfiles[:1],'PreSelection;1',branches=branches)
        tree = up3.open(inputfiles[0])['PreSelection;1']
        events = tree.pandas.df(branches=branches)
        goodname = samp
        zptdf  = events[events['ZCandidate_pt'] > zptcut]
        metdf   = zptdf[zptdf['metsuable'] > metcut]
        hptdf  = metdf[metdf['hCandidate_pt'] > hptcut]
        btdf   = hptdf[hptdf['hCandidate_'+btaggr] > float(btagwp)]

        #Actual Analysis
        #if not (valid or alphatest):
        if not valid:
            srup   = btdf[btdf['hCandidate_sd'] > 70.]
            bldf   = srup[srup['hCandidate_sd'] < 150.]#full blinded region
            srdf   = bldf[bldf['hCandidate_sd'] > 110.]#Higgs Peak
            lowsbh = btdf[btdf['hCandidate_sd'] <= 70.]
            lowsb  = lowsbh[lowsbh['hCandidate_sd'] > 30.]
            highsb = btdf[btdf['hCandidate_sd'] >= 150.]
            sbdf   = pd.concat([lowsb,highsb])
            totdf = pd.concat([bldf,lowsb,highsb])

