import gecorg_test as go
import ROOT
import argparse
import glob


if __name__=='__main__':
    #parser = argparse.ArgumentParser()
    #parser.add_argument("-f","--file", type=str,help = "fit diagnositic file with the output of the injection test")
    #parser.add_argument("-expsig","--expsig",type=float,help = "the signal strength of the expected signal")
    #parser.add_argument("-maxr","--maxr", type=float,help = "maxr")
    #parser.add_argument("-minr","--minr", type=float,help = "minr")
    #parser.add_argument("-b","--binwidth", type=float,help = "bin width")
    #args = parser.parse_args()
    print("Compiling fits")
    ROOT.gSystem.CompileMacro("basicfits.C","kc")#add f to force recompile
    ROOT.gSystem.Load("basicfits_C")


    print("Starting Analysis")
    f2ds = glob.glob("analysis_output_ZpAnomalon/2023-05-26/good2Dplots/*.root")

    #fits = []
    hists = []
    #errhists = []
    profs = []

    for f in f2ds:
        samp = f.split("-madgraph")[0].split("analysis_output_ZpAnomalon/2023-05-26/good2Dplots/")[-1]
        samp = samp.replace("-","_")
        samp = samp.replace(".","_")
        tf = ROOT.TFile(f)

        keys = tf.GetListOfKeys()
        keynames = [key.GetName() for key in keys]

        for key in keynames:
            if "nall" in key:
                continue
            print("Fitting {0} for {1}".format(key,samp))
            h = tf.Get(key)
            h.SetDirectory(0)
            h.SetName(samp+"_"+key+"_hist")
            prof = ROOT.profileHist2D(h,samp+"_"+key+"_tprof")
            prof.SetDirectory(0)
            #fit = ROOT.linearFit(h,samp+"_"+key+"_fit","QRE0+")
            #errhist = ROOT.linFitErrBandHist(h,samp+"_"+key+"_fit","RE0+")
            #errhist.SetDirectory(0)
            #fits.append(fit)
            hists.append(h)
            profs.append(prof)
            #errhists.append(errhist)
            #print(errhists)

        tf.Close()

    outname = go.makeOutFile("pileup_tranfers","recoquantToTrue",".root","","","","")
    fout = ROOT.TFile(outname,'recreate')

    for i,hist in enumerate(hists):
        #print(hist)
        #print(fits[i])
        #print(errhists[i])
        hist.Write()
        #fits[i].Write()
        #errhists[i].Write()
        profs[i].Write()
        
    fout.Close()
    
