import glob
import ROOT

files = glob.glob("mumu_ProperREOIDSF_ProperTrigger/Run2018A-17Sep2018-v1.EGamma_topiary*")

tf = ROOT.TFile(files[0])
hzpt = tf.Get("hzbuild_pt")
hzpt.Reset("ICESM")
hzeta = tf.Get("hzbuild_eta")
hzeta.Reset("ICESM")
hztrigpt = tf.Get("hzpasstrig_pt")
hztrigpt.Reset("ICESM")
hztrigeta = tf.Get("hzpasstrig_eta")
hztrigeta.Reset("ICESM")

hlmupt = tf.Get("hlmubuild_pt")
hlmupt.Reset("ICESM")
hlmueta = tf.Get("hlmubuild_eta")
hlmueta.Reset("ICESM")
hlmutrigpt = tf.Get("hlmupasstrig_pt")
hlmutrigpt.Reset("ICESM")
hlmutrigeta = tf.Get("hlmupasstrig_eta")
hlmutrigeta.Reset("ICESM")

hists = [hzpt,hzeta,hztrigpt,hztrigeta,hlmupt,hlmueta,hlmutrigpt,hlmutrigeta]
histnames = ["hzbuild_pt","hzbuild_eta","hzpasstrig_pt","hzpasstrig_eta","hlmubuild_pt","hlmubuild_eta","hlmupasstrig_pt","hlmupasstrig_eta"]


for f in files:
    fr = ROOT.TFile(f)
    for h,hist in enumerate(hists):
        histname = histnames[h]
        hist.Add(fr.Get(histname))
        
for h in hists:
    h.Rebin(1)
        
hzpt.SetTitle("Z(\mu\mu) reconstructed, 2018 EGamma")
hzpt.GetXaxis().SetTitle("p_{T}(Z)")
hzeta.SetTitle("Z(\mu\mu) reconstructed, 2018 EGamma")
hzeta.GetXaxis().SetTitle("\eta(Z)")
hztrigpt.SetTitle("Z(\mu\mu) reconstructed AND pass trigger, 2018 EGamma")
hztrigpt.GetXaxis().SetTitle("p_{T}(Z)")
hztrigeta.SetTitle("Z(\mu\mu) reconstructed AND pass trigger, 2018 EGamma")
hztrigeta.GetXaxis().SetTitle("\eta(Z)")


hlmupt.SetTitle("Z(\mu\mu) reconstructed, 2018 EGamma")
hlmupt.GetXaxis().SetTitle("p_{T}(LMU)")
hlmueta.SetTitle("Z(\mu\mu) reconstructed, 2018 EGamma")
hlmueta.GetXaxis().SetTitle("\eta(LMU)")
hlmutrigpt.SetTitle("Z(\mu\mu) reconstructed AND pass trigger, 2018 EGamma")
hlmutrigpt.GetXaxis().SetTitle("p_{T}(LMU)")
hlmutrigeta.SetTitle("Z(\mu\mu) reconstructed AND pass trigger, 2018 EGamma")
hlmutrigeta.GetXaxis().SetTitle("\eta(LMU)")

zpteff = ROOT.TEfficiency(hztrigpt,hzpt)
zpteff.SetTitle("Event Efficiency versus Z p_{T}")
zetaeff = ROOT.TEfficiency(hztrigeta,hzeta)
zetaeff.SetTitle("Event Efficiency versus Z Eta")

lmupteff = ROOT.TEfficiency(hlmutrigpt,hlmupt)
lmupteff.SetTitle("Event Efficiency versus leading electron p_{T}")
lmuetaeff = ROOT.TEfficiency(hlmutrigeta,hlmueta)
lmuetaeff.SetTitle("Event Efficiency versus leading electron Eta")

tc = ROOT.TCanvas("tc","tc",1200,800)
tc.Divide(3,2)
tc.cd(1)
hzpt.Draw()
tc.cd(2)
hztrigpt.Draw()
tc.cd(3)
zpteff.Draw()
tc.cd(4)
hzeta.Draw()
tc.cd(5)
hztrigeta.Draw()
tc.cd(6)
zetaeff.Draw()
tc.SaveAs("zmumu_versus_efficiency_out_2018.png")

tc2 = ROOT.TCanvas("tc","tc",1200,800)
tc2.Divide(3,2)
tc2.cd(1)
hlmupt.Draw()
tc2.cd(2)
hlmutrigpt.Draw()
tc2.cd(3)
lmupteff.Draw()
tc2.cd(4)
hlmueta.Draw()
tc2.cd(5)
hlmutrigeta.Draw()
tc2.cd(6)
lmuetaeff.Draw()
tc2.SaveAs("lmu_versus_efficiency_out_2018.png")

tc3 = ROOT.TCanvas("tc3","tc3",1200,400)
tc3.Divide(2,1)
tc3.cd(1)
zpteff.Draw()
tc3.cd(2)
zetaeff.Draw()
tc3.SaveAs("efficiencies_v_zmumu_2018.png")
    
tc3 = ROOT.TCanvas("tc3","tc3",1200,400)
tc3.Divide(2,1)
tc3.cd(1)
lmupteff.Draw()
tc3.cd(2)
lmuetaeff.Draw()
tc3.SaveAs("efficiencies_v_lmumu_2018.png")
    
