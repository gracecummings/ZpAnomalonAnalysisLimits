import pandas as pd
#import uproot3 as up3
import gecorg_test as go
import numpy as np
import ROOT 

def makeHistogramSymErrs(dfgen):
    #get hist params
    #1.) Check if the binning is regular
    diffs = dfgen['ptmax']-dfgen['ptmin']
    meanval = diffs[:-1].mean()
    if (meanval - diffs[0]) != 0:
        print("Binning not regular, need non-uproot machinery")

    #2.) Maeke the histogram
    nbins = len(diffs)
    lowedge = dfgen['ptmin'][0]
    highedge = dfgen['ptmin'][nbins-1]+diffs[0]#replaces last high edge (inf) one bin width higher
    histname = str(dfgen['year'][0])+'sf'
    print("Making histgram ",histname)
    histout = ROOT.TH1F(histname,histname,nbins,lowedge,highedge)
    for i in range(nbins):
        print("The min pt from df: ",dfgen['ptmin'][i])
        print("    sf: ",dfgen['sf'][i])
        print("  sferr: ",dfgen['sfhigh'][i])
        print("    The min pt from the histogram: ",histout.GetBinLowEdge(i+1))
        #Make sure the errors are symmetric
        errdiff = dfgen['sfhigh'][i]-dfgen['sflow'][i]
        if errdiff != 0:
            print("Bin errors are not symmetric, aborting")
            break
        histout.SetBinContent(i+1,dfgen['sf'][i])
        histout.SetBinError(i+1,dfgen['sfhigh'][i])
        print("Hist content: ",histout.GetBinContent(i+1))
        print("Hist errors: ",histout.GetBinError(i+1))

    histout.Write()

def makeHistogramsErrAndSfSep(dfgen):
    #get hist params
    #1.) Check if the binning is regular
    diffs = dfgen['ptmax']-dfgen['ptmin']
    meanval = diffs[:-1].mean()
    binedges = np.array(dfgen['ptmin'],dtype=float)
    binedges = np.append(binedges,5000.0)

    #2.) Make the histograms
    nbins = len(diffs)
    lowedge = dfgen['ptmin'][0]
    highedge = dfgen['ptmin'][nbins-1]+diffs[0]#replaces last high edge (inf) one bin width higher
    hnamesf = str(dfgen['year'][0])+'sf'
    hnameunchigh = str(dfgen['year'][0])+'uncUp'
    hnameunclow = str(dfgen['year'][0])+'uncDown'
    print("Making histgrams for ",str(dfgen['year'][0]))

    houtsf = ROOT.TH1F(hnamesf,hnamesf,nbins,binedges)
    houtuncup = ROOT.TH1F(hnameunchigh,hnameunchigh,nbins,binedges)
    houtuncdwn = ROOT.TH1F(hnameunclow,hnameunclow,nbins,binedges)

    for i in range(nbins):
        print("    The min pt from the histogram: ",houtsf.GetBinLowEdge(i+1))
        houtsf.SetBinContent(i+1,dfgen['sf'][i])
        houtuncup.SetBinContent(i+1,dfgen['sfhigh'][i])
        houtuncdwn.SetBinContent(i+1,dfgen['sflow'][i])
        print("Hist sf content: ",houtsf.GetBinContent(i+1))
        print("sf unc up content: ",houtuncup.GetBinContent(i+1))
        print("sf unc dwn content: ",houtuncdwn.GetBinContent(i+1))

    houtsf.Write()
    houtuncup.Write()
    houtuncdwn.Write()

if __name__=='__main__':

    #for ttbar sf
    sfttdf = pd.read_csv('deepak8_hbbSF_final_TTSF.csv')
    ttsfrootname = go.makeOutFile('DeepAK8MassDecorrelZHbbvQCD','ttscalefactors','.root','','','','')
    ttsfroot = ROOT.TFile(ttsfrootname,"recreate")
    df16 = sfttdf[sfttdf['year'] == 2016]
    df17 = sfttdf[sfttdf['year'] == 2017].copy()
    df17 = df17.reset_index(drop=True)#start the indexing at 0 for ease
    df18 = sfttdf[sfttdf['year'] == 2018].copy()
    df18 = df18.reset_index(drop=True)
    
    makeHistogramsErrAndSfSep(df16)
    makeHistogramsErrAndSfSep(df17)
    makeHistogramsErrAndSfSep(df18)

    ttsfroot.Write()
    ttsfroot.Close()

    
    #For non-ttbar errors

    #read in the sf
    sfdf = pd.read_csv('deepak8_hbbSF_final_noTTSF.csv')

    #make the root file to hold the sf in th1s
    sfrootname = go.makeOutFile('DeepAK8MassDecorrelZHbbvQCD','scalefactors','.root','','','','')
    sfroot = ROOT.TFile(sfrootname,"recreate")

    df16 = sfdf[sfdf['year'] == 2016]
    df16lp = df16[df16['wp'] == 0.8]
    df17 = sfdf[sfdf['year'] == 2017]
    df17lp = df17[df17['wp'] == 0.8].copy()
    df17lp = df17lp.reset_index(drop=True)#start the indexing at 0 for ease
    df18 = sfdf[sfdf['year'] == 2018]
    df18lp = df18[df18['wp'] == 0.8].copy()
    df18lp = df18lp.reset_index(drop=True)

    #makeHistogramSymErrs(df16lp)
    #makeHistogramSymErrs(df17lp)
    #makeHistogramSymErrs(df18lp)

    makeHistogramsErrAndSfSep(df16lp)
    makeHistogramsErrAndSfSep(df17lp)
    makeHistogramsErrAndSfSep(df18lp)

    sfroot.Write()
    sfroot.Close()
    
