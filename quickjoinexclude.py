import uproot as up4
import uproot3 as up3
import pandas as pd
import numpy as np
import boost_histogram as bh
import argparse
import glob
import gecorg_test as go

parser = argparse.ArgumentParser()

if __name__=='__main__':
    parser.add_argument("-f","--sample",help = "sample file")
    parser.add_argument("-d","--directory",type=str,help = "where are your topiary plots?")
    parser.add_argument("-sr","--signalregion",type=bool,help = "do you want a signal region plot?")
    parser.add_argument("-tot","--comboregion",type=bool,help = "do you want combined SR and SB?")
    parser.add_argument("-v","--validation",type=bool,help = "validation region bounds?")
    parser.add_argument("-syst","--systematics",type=str)
    parser.add_argument("-c","--chan",type=str)
    args = parser.parse_args()

    samp   = args.sample
    topdir = args.directory
    channel = args.chan

    if (("Run" in samp) and ("emu" in channel)):
        if (("Single" not in samp) or ("EGamma" not in samp)):
            print("Going to combine datasets")
            muf = glob.glob(topdir+'/'+samp+'*SingleMuon*_topiary_'+channel+'*.root')
            euf = glob.glob(topdir+'/'+samp+'*SingleElectron*_topiary_'+channel+'*.root')
            print(muf)
            print(euf)

            branches = [b'RunNum',
                        b'LumiBlockNum',
                        b'EvtNum',
                        ]

            mtree = up3.open(muf[0])['PreSelection;1']
            mevents = mtree.pandas.df(branches=branches)
            etree = up3.open(euf[0])['PreSelection;1']
            eevents = etree.pandas.df(branches=branches)
            print(mevents)
            print(eevents)

            frames = [mevents,eevents]
            mixdf = pd.concat(frames)

            print(mixdf)

            mixdf.drop_duplicates(subset = ['RunNum','LumiBlockNum','EvtNum'])

            print(mixdf)
