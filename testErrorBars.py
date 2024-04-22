import ROOT

h = ROOT.TH1F("test","test",3,0,4)

h.SetBinContent(1,3)
h.SetBinContent(2,2)
h.SetBinContent(3,1)

g = h.Clone()

g.SetBinError(1,1)
g.SetBinError(2,0.5)
g.SetBinError(3,0.25)
g.SetLineColor(ROOT.kBlue-9)
g.SetFillColor(ROOT.kBlack)
g.SetFillStyle(3245)

h.SetFillColor(ROOT.kBlue-9)
h.SetLineColor(ROOT.kBlue-9)

tc = ROOT.TCanvas("testc","testc",600,600)
tc.Draw()
tc.cd()
h.Draw()
g.Draw('same,E2')

tc.SaveAs("test_hist.png")
