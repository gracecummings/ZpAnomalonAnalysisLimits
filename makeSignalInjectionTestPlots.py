import gecorg_test as go
import ROOT
import argparse

ROOT.gROOT.SetBatch(ROOT.kTRUE)

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file", type=str,help = "fit diagnositic file with the output of the injection test")
    parser.add_argument("-expsig","--expsig",type=float,help = "the signal strength of the expected signal")
    parser.add_argument("-maxr","--maxr", type=float,help = "maxr")
    parser.add_argument("-minr","--minr", type=float,help = "minr")
    parser.add_argument("-b","--binwidth", type=float,help = "bin width")
    args = parser.parse_args()

    #command line inputs
    f = args.file
    expsig = args.expsig
    maxr = args.maxr
    minr = args.minr
    binw = args.binwidth
    if minr and maxr:
        binnumber = int(abs(maxr)+abs(minr)/binw)

    #open the files and get the basic objects
    tf = ROOT.TFile(f)
    tree = tf.Get("tree_fit_sb")
    outf = go.makeOutFile(f.split("fitDiagnostics")[-1].split(".root")[0],"signalinjectionplots",".root","100","300","75","8E-10")
    savef = ROOT.TFile(outf,"recreate")

    #define the histograms to fill
    #hr = ROOT.TH1F("hr","r {fit_status==0}",binnumber,minr,maxr)
    hr = ROOT.TH1F("hr","r {fit_status==0}",60,-10,10)
    #hdr = ROOT.TH1F("hdr","r-r_exp {fit_status==0}",20,-5,5)
    #hdr = ROOT.TH1F("hdr","r-r_exp {fit_status==0}",60,-10,10)
    hdr = ROOT.TH1F("hdr","r-r_exp {fit_status==0}",60,-10,10)
    #hdrerr = ROOT.TH1F("hdrerr","(r -r_exp)/rErr, standard error, {fit_status==0})",20,-5,5)
    hdrerrhl = ROOT.TH1F("hdrerrhl","(r -r_exp)/rErr, hilow error, {fit_status==0})",60,-10,10)
    hrge = ROOT.TH1F("hrge","r {fit_status>=0}",60,-10,-10)
    #hdrge = ROOT.TH1F("hdrge","r-r_exp {fit_status>=0}",20,-5,5)
    hdrge = ROOT.TH1F("hdrge","r-r_exp {fit_status>=0}",60,-10,10)
    #hdrerrge = ROOT.TH1F("hdrerrge","(r -r_exp)/rErr, standard error, {fit_status>=0}",20,-5,5)
    hdrerrgehl = ROOT.TH1F("hdrerrgehl","(r -r_exp)/rErr, hilow error, {fit_status>=0}",60,-10,10)
    #hhilowerr = ROOT.TH1F("hhilowerr","hilow error, {fit_status>=0}",60,-5,5)

    #fill the histograms
    #comparators = ["==",">="]
    #for comp in compartors:
    #    cutstring = "fit_status "+comp+" 0"

    tree.Draw("r>>hr","fit_status == 0","")
    tree.Draw("r>>hrge","fit_status >= 0","")
    deltarstr = "(r-"+str(expsig)+")>>hdr"
    deltarstrge = "(r-"+str(expsig)+")>>hdrge"
    tree.Draw(deltarstr,"fit_status == 0","")
    tree.Draw(deltarstrge,"fit_status >= 0","")
    #drerrstr = "(r-"+str(expsig)+")/rErr>>hdrerr"
    #drerrstrge = "(r-"+str(expsig)+")/rErr>>hdrerrge"
    #tree.Draw(drerrstr,"fit_status == 0","")
    #tree.Draw(drerrstrge,"fit_status >= 0","")

    rerrstr = "(rHiErr*((r-"+str(expsig)+")<0)+rLoErr*((r-"+str(expsig)+")>0))"#first hilow err string
    drerrstrhl = "(r-"+str(expsig)+")/"+rerrstr+">>hdrerrhl"
    drerrstrgehl = "(r-"+str(expsig)+")/"+rerrstr+">>hdrerrgehl"
    tree.Draw(drerrstrhl,"fit_status == 0","")
    tree.Draw(drerrstrgehl,"fit_status >= 0","")
    #tree.Draw(rerrstr+" >>hhilowerr","fit_status >= 0","")

    #hists = [hr,hdr,hdrerr,hrge,hdrge,hdrerrge,hdrerrhl,hdrerrgehl]
    hists = [hr,hdr,hrge,hdrge,hdrerrhl,hdrerrgehl]
    #ROOT.gStyle.SetOptFit(1111)
    f1 = ROOT.TF1("f1","[0]*exp(-0.5*((x-[1])/[2])^2)")
    
    for hist in hists:
        #hist.Fit("gaus")
        ROOT.gStyle.SetOptFit(1111)#fit window on all but first hist
        f1.SetParameter(0,hist.GetMaximum())
        f1.SetParameter(1,hist.GetBinLowEdge(hist.GetMaximumBin())+0.5*hist.GetBinWidth(hist.GetMaximumBin()))
        f1.SetParameter(2,hist.GetBinLowEdge(hist.FindLastBinAbove(hist.GetMaximum()/2))+0.5*hist.GetBinWidth(hist.FindBin(hist.GetMaximum()/2)))
        hist.Fit(f1)
        hist.GetYaxis().SetTitle("Number of toys")
        print("Making ",hist.GetName())
        if "drerr" in hist.GetName():
            hist.GetXaxis().SetTitle("r-r_{exp}/\sigma_{r}")
        elif "dr" in hist.GetName():
            hist.GetXaxis().SetTitle("r-r_{exp}")
    savef.Write()
    
