#This script transfer work John Hakala did for a 2D interpolator. Method is his, as are the functions. G. Cummings made it into a user friendly script
import array
from scipy.interpolate import RBFInterpolator as rbf
import ROOT as r
import numpy as np
import os
from datetime import date
import glob

def makeOutFile(sampstring,descrip,ftype,zptcut,hptcut,metcut,btagwp):
    if not os.path.exists("analysis_output_ZpAnomalon/"+str(date.today())+"/"):
        os.makedirs("analysis_output_ZpAnomalon/"+str(date.today())+"/")
    outFile = "analysis_output_ZpAnomalon/"+str(date.today())+"/"+sampstring+"_"+descrip+"_Zptcut"+zptcut+"_Hptcut"+hptcut+"_metcut"+metcut+"_btagwp"+btagwp+ftype
    return outFile

def getRealSignals(f):
    hists = {}
    tf=r.TFile(f)
    histKeys = [key.GetName() for key in tf.GetListOfKeys() if not "_lumixs" in key.GetName()]
    for histKey in histKeys:
        massesKey = tuple(float(x) for x in histKey.replace("Zp", "").split("NS")[0].split("ND"))
        hist = tf.Get(histKey).Clone()
        hist.SetDirectory(0)
        r.SetOwnership(hist, r.kFALSE)
        hists[massesKey] = hist
    return hists

def getQuantileInterps(mZp,mND,quantrbf):
    interpQuantiles = {}
    for quantileKey in quantrbf.keys():
        interpQuantiles[quantileKey]=quantrbf[quantileKey]([[mZp, mND]])
    return interpQuantiles

if __name__=="__main__":

    files = glob.glob("analysis_output_ZpAnomalon/2022-10-02/Run2_161718_ZllHbbMET_syst*")

    for f in files:
        ####file with the histograms to seed the interpolations
        #f = "analysis_output_ZpAnomalon/2022-10-05/Run2_161718_ZllHbbMET_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100_Hptcut300_metcut75_btagwp8E-10.root"
        systr = f.split("Run2_161718_ZllHbbMET_")[-1].split("_Zptcut100_Hptcut300_metcut75_btagwp8E-10.root")[0]
        print("Doing signal interpolations for ",systr)
        #systr = "systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco"
    
        ####Get the real signal histograms from a file
        #returns a dictionary of hists with the masses in a float
        hists = getRealSignals(f)
        masspoints = set(hists.keys())
        totalgrid = set([tuple([x,400+y*200]) for x in range(1500,6000,500) for y in range(0,int(x/400)-1)])
        missing = totalgrid - masspoints
        
        ####Build quantile info
        # the precision here can be lowered
        nQuantiles = 10000
        empty = [0. for i in range(0,nQuantiles)]
        # % corresponding to the quantiles
        xes = [(i)*1/float(nQuantiles) for i in range(0, nQuantiles)]
        arrX = array.array('d', xes)
        quantiles  = {}
        normalizations = {}
        
        mZps = []
        mNDs = []
        norms = []
        
        #gather the info for each generated quantile from the original hists
        #print("Check hist number, in debug mode")
        #debug = 0
        for key in hists:
            #if debug > 3:
            #    break
            #debug+=1
            quantiles[key]=array.array('d', empty)
            #get values at each quantile, for each signal point
            #quantile value stored in `quantiles` dictionary
            hists[key].GetQuantiles(nQuantiles, quantiles[key], arrX)
            #normalizations[key] = hists[key].GetSumOfWeights()#og way
            mZps.append(key[0])
            mNDs.append(key[1])
            norms.append(hists[key].GetSumOfWeights())
        
        ###Build arrays of the quantile info
        quantileArrays = {}
        for iQuantile in range(0, nQuantiles):
            # % quantile is the key
            # value is a list of the quantile positiion for each original signal
            #`v` is the list of quantile outputs for each signal
            quantileArrays[xes[iQuantile]]=[v[iQuantile] for v in quantiles.values()]

        ###Shift stuff to be in an intuitive way to plot
        points = np.transpose([mZps, mNDs])

        ###Build the interpolation objets for shape and norm
        normRbf  = rbf(points, norms)
        rbfs = {}
        for quantileKey in quantileArrays.keys():
            #for each quantile, prepare to do a 2D interpolation
            rbfs[quantileKey]=rbf(points, quantileArrays[quantileKey], None, 0.0)

        ###Find the interpolated values we need
        pointsToTest=list(missing)#[(2500., 800.), (1800., 800.), (1950., 875.), (2100., 825.), (2100., 950.), (3250., 750.), (3750., 900.), (3500., 1350.)]

        ####Define params to make the interpolated histograms
        nBinsInterp = 100
        upEdgeInterp = 10000
        htest = list(hists.values())[0]#takes the first histogram in hist array
        rebinFactor = int(htest.GetXaxis().GetBinWidth(htest.GetNbinsX()//2) / (upEdgeInterp/float(nBinsInterp)))#complicated way to get 2. Finds the bin width of a random bin in the orignal hist, and finds the relationship to the new bin width

        interpHists = {}
        cans= []

        #namer = 0
        outfname = makeOutFile("interpolation",systr,".root","100.0","300.0","75.0","0.8")
        outf = r.TFile(outfname,"RECREATE")
        print("Writing interpolations to ",outfname)

        for iPoint, pointToTest in enumerate(pointsToTest):
            #cans.append(r.TCanvas())
            #cans[-1].cd()
            newQuantiles = getQuantileInterps(*pointToTest,rbfs)
            #nFromHist = list(hists.values())[0].GetSumOfWeights()
            nToGenerate = 10000
            scaleFactor = normRbf([[pointToTest[0], pointToTest[1]]])/float(nToGenerate)

            interpHists[pointToTest] = r.TH1F(f"interpHist{iPoint}", "interpolated hist", nBinsInterp, 0, upEdgeInterp)
            interpQuantiles = [q[0] for q in newQuantiles.values()]
            for i in range(0, len(interpQuantiles)-1):
                for j in range(0,nToGenerate//nQuantiles):
                    interpHists[pointToTest].Fill(r.gRandom.Uniform(interpQuantiles[i], interpQuantiles[i+1]))
            interpHists[pointToTest].Scale(scaleFactor)
            interpHists[pointToTest].Rebin(int(rebinFactor))
            #interpHists[pointToTest].Draw("HIST SAME")
            #cans[-1].SaveAs("canvas"+str(namer)+".png")
            interpHists[pointToTest].SetName("Zp"+str(round(pointToTest[0]))+"ND"+str(round(pointToTest[1]))+"NS200")
            interpHists[pointToTest].Write()
            #namer+=1

