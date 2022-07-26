import ROOT
import gecorg_test as go

tf = ROOT.TFile("analysis_output_ZpAnomalon/2022-07-23/Run2_161718_ZllHbbMET_extended_poisson_mumu_Zp5500ND1800NS200_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root")

#samp = "Zp5500ND1800NS200"
#samp = "TT"
samp = "VV"
keys = tf.GetListOfKeys()
keys = [key.GetName() for key in keys]
#print(keys)
ttuncs = [key for key in keys if samp+"_"+samp+"_Stats" in key]
cols = go.colsFromPalette(ttuncs,ROOT.kRainBow)

tc = ROOT.TCanvas("tc","tc",1200,400)
tc.Divide(2,1)
tc.cd(1)
leg = ROOT.TLegend(0.1,0.1,0.9,0.9)
htt = tf.Get(samp)
htt.SetLineColor(ROOT.kBlack)
htt.SetLineWidth(3)
leg.AddEntry(htt,samp,"l")

maxes = []
for i,unc in enumerate(ttuncs):
    h = tf.Get(unc)
    h.SetStats(0)
    h.GetYaxis().SetRangeUser(0,.7)
    h.SetLineColor(cols[i])
    leg.AddEntry(h,unc,"l")
    if i == 0:
        h.Draw("hist")
    else:
        h.Draw("histsame")

    maxes.append(h.GetMaximum())

maxi = max(maxes)
print(maxi)

htt.GetYaxis().SetRangeUser(0,.7)
htt.Draw("histsame")
tc.cd(2)
leg.Draw()


tc.SaveAs("testofstatsuncn_newstats_"+samp+".png")

