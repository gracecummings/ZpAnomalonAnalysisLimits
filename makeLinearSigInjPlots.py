import ROOT
import matplotlib.pyplot as plt
import numpy as np
import glob
import argparse
import os
from datetime import date

ROOT.gROOT.SetBatch(ROOT.kTRUE)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-mp","--masspoint", type=str,help = "mass point string, ex: 'Zp5500ND1800NS200'")
    parser.add_argument("-date","--date",type=str,help = "date directory where stuff is stored")
    parser.add_argument("-n","--tag",type=str,help = "tag to append to filename")
    args = parser.parse_args()


    mp = args.masspoint
    tag = args.tag
    infs = glob.glob("analysis_output_ZpAnomalon/"+args.date+"/"+mp+"*_expectedsignal*signalinjectionplots_Zptcut100_Hptcut300_metcut75_btagwp8E-10.root")
    
    pdict = {}
    injs = []

    for f in infs:
        inj = f.split("_expectedsignal")[-1].split("_")[0]
        pdict[float(inj)] = f
        injs.append(float(inj))

    injs = sorted(injs)
    rmeas = []
    errs  = []

    for p in injs:
        print('Injected signal: ',p)
        tf = ROOT.TFile(pdict[p])
        #h = tf.Get("hrge")
        h = tf.Get("hr")
        f1 = h.GetFunction("f1")
        mu = f1.GetParameter(1)
        sigma = abs(f1.GetParameter(2))
        print("The mean of the guassian: ",mu)
        print('The sigma of the gaussian: ',sigma)
        rmeas.append(mu)
        errs.append(sigma)
        
    x = np.linspace(-.2,round(max(injs)*2,1),5)
    
    fig = plt.figure()

    plt.plot(x,x,label="y=x")
    plt.errorbar(injs,rmeas,yerr=errs,label="r v. r injected, 1 sigma bars",fmt=".k")
    plt.legend()
    plt.title(mp)
    plt.xlabel("signal strength (r), injected")
    plt.ylabel("signal strenght (r), measured")

    if not os.path.exists("analysis_output_ZpAnomalon/"+str(date.today())+"/"):
        os.makedirs("analysis_output_ZpAnomalon/"+str(date.today())+"/")

    outFilepng = "analysis_output_ZpAnomalon/"+str(date.today())+"/linear_bias_test_"+mp+"_ge_"+tag+".png"
    outFilepdf = "analysis_output_ZpAnomalon/"+str(date.today())+"/linear_bias_test_"+mp+"_ge_"+tag+".pdf"
    
    plt.savefig(outFilepng)
    plt.savefig(outFilepdf)
    #plt.show()


