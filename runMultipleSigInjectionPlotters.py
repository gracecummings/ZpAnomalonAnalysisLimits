import gecorg_test as go
import ROOT
import argparse
import glob
import subprocess
from datetime import date

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--direc", type=str,help = "directory with files with the output of the injection test")
    parser.add_argument("-plot","--doplot", type=bool,help = "making the sig injection plots")
    parser.add_argument("-png","--dopngs", type=bool,help = "making the sig injection plot pngs")
    args = parser.parse_args()

    if args.direc and args.doplot:
        print("doing the pull plots with fits for sig injection")
        inj_res = glob.glob(args.direc+"/fitDiagnostics*expectedsignal*.root")

        for f in inj_res:
            inj_val = f.split("expectedsignal")[-1].split("_")[0]
            subprocess.run(["python","makeSignalInjectionTestPlots.py","-f",f,"-expsig",inj_val])
    elif args.direc and args.dopngs:
        print("making already made pull plots and making pngs")
        plot_fs = glob.glob(args.direc+"/Zp*expectedsignal*signalinjectionplots_*.root")

        for f in plot_fs:
            name1= f.split("_Zptcut100_Hptcut300_metcut75_btagwp8E-10.root")[0]
            print(args.direc)
            name = name1.split(args.direc)[-1]
            #lim = n1.split("_"+repstr)[0].split("expectedsignal")[-1]
            tf = ROOT.TFile(f)
            keys = tf.GetListOfKeys()
            nkeys = [key.GetName() for key in keys]
            keyset = set(nkeys)
            #badkeys = set(["hdrerr","hdrerrge","hrge","hrge","hdrge","hdrerrgehl"])
            badkeys = set(["hdrerr","hdrerrge","hrge"])
            names = {"hdr":"diff","hdrerrhl":"pull","hdrerrgehl":"gepull","hdrge":"gediff","hr":"rdist"}
            pltkeys = keyset-badkeys
            for key in pltkeys:
                hname = key
                h = tf.Get(hname)
                tc = ROOT.TCanvas("tc","tc",1060,800)
                tc.Draw()
                tc.cd()
                h.Draw()
                tc.cd()
                outname = go.makeOutFile(name,hname,".png","100","300","75","8E-10")
                #print(outname)
                tc.SaveAs(outname)

        
    else:
        print("you forgot to put a directory")    
