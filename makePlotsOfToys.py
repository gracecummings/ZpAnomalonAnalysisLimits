import ROOT
import gecorg_test as go
fname = "analysis_output_ZpAnomalon/2022-05-24/higgsCombineZp4000ND800NS200_20220523_expsig0.GenerateOnly.mH120.123456.root"

f = ROOT.TFile.Open(fname)

toydir = f.Get("toys")
keys = toydir.GetListOfKeys()
keys = [k.GetName() for k in keys]
toys  = [k for k in keys if "snapshot" not in k]
canvases = []

print("Beginning to plot, there are {0} toys".format(len(toys)))
for toy in toys:
    toycont = f.Get("toys/"+toy)
    x = ROOT.RooRealVar("CMS_th1x","CMS_th1x",0,8)
    fr = x.frame()
    toycont.plotOn(fr)
    can = ROOT.TCanvas(toy,toy)
    can.Draw()
    can.cd()
    fr.Draw()
    #canout = go.makeOutFile(fname.split("higgsCombine")[-1].split(".root")[0],toy,".png","","","","")
    #can.SaveAs(canout)
    canvases.append(can)

#f.Close()
of = go.makeOutFile(fname.split("higgsCombine")[-1].split(".root")[0],"plottedtoys",".root","","","","")
print("Finished plotting, savining plotted toys in ",of)
fout = ROOT.TFile(of,"recreate")
for c in canvases:
    c.Write()
fout.Close()
f.Close()
