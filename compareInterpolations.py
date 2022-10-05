import ROOT

tf0 = ROOT.TFile("analysis_output_ZpAnomalon/2022-10-05/interpolation_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root","READ")
tf1 = ROOT.TFile("~/Downloads/interpolated-hists.root","READ")

keysgec = tf0.GetListOfKeys()
keysog  = tf1.GetListOfKeys()
keysgec_named = [key.GetName() for key in keysgec]
keysog_named  = [key.GetName() for key in keysog]

gechists = {}

for key in keysgec_named:
    gechists[key] = tf0.Get(key)

tfout = ROOT.TFile("analysis_output_ZpAnomalon/2022-10-05/interpolation_divcheck.root","RECREATE")
for key in keysog_named:
    shortn = key.split("-")[-1]
    hist = tf1.Get(key)
    histdiv = hist.Clone()
    histdiv.Divide(histdiv,gechists[shortn])
    histdiv.SetTitle("Div. Original Interps w/ New:"+str(shortn))
    histdiv.Write()

    




