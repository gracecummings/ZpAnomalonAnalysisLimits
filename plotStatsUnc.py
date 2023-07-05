import ROOT
import gecorg_test as go

tf = ROOT.TFile("analysis_output_ZpAnomalon/2023-05-09/Run2_161718_ZllHbbMET_vvsixpercshift_mumu_Zp2000ND400NS200_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root")

samp = "Zp2000ND400NS200"
#samp = "Zp4000ND800NS200"
#samp = "Zp5500ND1800NS200"
#samp = "VV"
#samp = "TT"
keys = tf.GetListOfKeys()
keys = [key.GetName() for key in keys]
#print(keys)
ttuncs = [key for key in keys if samp+"_"+samp+"_Stats" in key]
cols = go.colsFromPalette(ttuncs,ROOT.kRainBow)

tc = ROOT.TCanvas("tc","tc",1200,400)
tc.Divide(2,1)
tc.cd(1)

leg = ROOT.TLegend(0.1,0.1,0.8,0.9)
htt = tf.Get(samp)
htt2 = htt.Clone()
htt.SetLineColor(ROOT.kBlack)
htt.SetLineWidth(3)
leg.AddEntry(htt,samp,"l")

maxes = []
for i,unc in enumerate(ttuncs):
    h = tf.Get(unc)
    h.SetStats(0)
    #h.GetYaxis().SetRangeUser(-2,20)
    h.SetLineColor(cols[i])
    leg.AddEntry(h,unc,"l")
    if i == 0:
        h.Draw("hist")
    else:
        h.Draw("histsame")

    maxes.append(h.GetMaximum())

#maxi = max(maxes)
#print(maxi)

#htt.GetYaxis().SetRangeUser(-2,20)
htt.Draw("histsame")

tc.cd(2)
leg.Draw()


tc.SaveAs("testofstatsucn_1800lowbin_"+samp+".png")

