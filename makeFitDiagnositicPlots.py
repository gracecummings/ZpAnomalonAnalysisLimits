import ROOT
import argparse

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-fitdiag","--fitdiag", type=str,help = "fitdiagnositc output rootfile")
    parser.add_argument("-higgscomb","--higgscomb", type=str,help = "higgscombine output rootfile")
    parser.add_argument("-dir","--outdir", type=str,help = "directory to save the plots")
    args = parser.parse_args()

    #Finding a name
    #name from limit file
    protoname = args.higgscomb.split("higgsCombine")[-1]
    name = protoname.split(".FitDiagnostics.mH120.root")[0]
    #name from fit diag file
    protoname2 = args.fitdiag.split("fitDiagnostics")[-1]
    name2 = protoname2.split(".root")[0]
    if name == name2:
        print("Parsing values for ",name)
    else:
        print("!!!!!!!!!!You input files do not align, check before you trust the results!!!!!!!!!!")
    
    f = ROOT.TFile.Open(args.higgscomb)
    tree = f.Get("limit")
    tree.GetEntry(2)#The 50% quantiles for median limit
    limit = tree.limit

    print("The Limit found for this fit diagnositc run is: ",limit)
    print("Now making the plots of the post fit distributions.")

    ffits = ROOT.TFile.Open(args.fitdiag)
    treefits = ffits.Get("shapes_fit_b/Zp4000ND800NS200_mumu")
    keys = [x.GetName() for x in treefits.GetListOfKeys()]

    for key in keys:
        tc = ROOT.TCanvas("tc",key,533,450)
        h = treefits.Get(key)
        tc.Draw()
        tc.cd()
        h.Draw()
        pngname = args.outdir+"/"+name+"_"+key+"_postfit.png"
        tc.SaveAs(pngname)

    print("Now making the plots of the prefit distributions.")
    treefitspre = ffits.Get("shapes_prefit/Zp4000ND800NS200_mumu")
    keyspre = [x.GetName() for x in treefits.GetListOfKeys()]

    for key in keyspre:
        tc = ROOT.TCanvas("tc",key,533,450)
        h = treefitspre.Get(key)
        tc.Draw()
        tc.cd()
        h.Draw()
        pngname = args.outdir+"/"+name+"_"+key+"_prefit.png"
        tc.SaveAs(pngname)
