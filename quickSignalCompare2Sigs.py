import gecorg_test as go
import pandas as pd
import ROOT

signal = "Fall17.DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_upout_sideband_alphat"
signalerr = "Fall17.DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8_selected_errors_sideband_alphat"
htestname  = "h_zp_jigm"

#fnom = ROOT.TFile("mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref_alphatest/"+signal+"_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root")
#fnomerrs = "mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref_alphatest/"+signalerr+"_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.pkl"

#ftest = ROOT.TFile("analysis_output_ZpAnomalon/2023-05-30/pumapednominal/"+signal+"_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root")
#ftesterrs = "analysis_output_ZpAnomalon/2023-05-30/pumapednominal/"+signalerr+"_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.pkl"

ftest = ROOT.TFile("analysis_output_ZpAnomalon/pumapdwnwithnorm/"+signal+"_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root")
ftesterrs = "analysis_output_ZpAnomalon/pumapdwnwithnorm/"+signalerr+"_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.pkl"
fnom = ROOT.TFile("analysis_output_ZpAnomalon/pumapwithnorm/"+signal+"_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root")
fnomerrs = "analysis_output_ZpAnomalon/pumapwithnorm/"+signalerr+"_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.pkl"


#fnom = ROOT.TFile("analysis_output_ZpAnomalon/2023-05-24/pilupweights_nom/"+signal+"_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root")
#fnomerrs = "analysis_output_ZpAnomalon/2023-05-24/pilupweights_nom/"+signalerr+"_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.pkl"

hnom = fnom.Get(htestname)
htest = ftest.Get(htestname)

nomerrs = pd.read_pickle(fnomerrs)
testerrs = pd.read_pickle(ftesterrs)

for ibin in range(hnom.GetNbinsX()+1):
    if ibin == 0:
        continue
    else:
        hnom.SetBinError(ibin,nomerrs[htestname][ibin-1])
        htest.SetBinError(ibin,testerrs[htestname][ibin-1])

ROOT.gStyle.SetOptStat(0)
hnom.GetYaxis().SetTitle("Events")
hnom.GetXaxis().SetTitle("RJR Z Prime Estimator (GeV)")
htest.SetLineColor(ROOT.kRed)

#print(
print("Nominal integral:         ",hnom.Integral(hnom.FindBin(1800),hnom.FindBin(10000)))
print("Pileup weighted integral: ",htest.Integral(htest.FindBin(1800),htest.FindBin(10000)))
print("perc diff:                ",(hnom.Integral(hnom.FindBin(1800),hnom.FindBin(10000))-htest.Integral(htest.FindBin(1800),htest.FindBin(10000)))/hnom.Integral(hnom.FindBin(1800),hnom.FindBin(10000))*100)

hdiv = htest.Clone()
hdiv.Divide(htest,hnom)


hdiv.GetXaxis().SetTitleSize(0.07)
hdiv.GetXaxis().SetTitleOffset(0.65)
hdiv.GetXaxis().SetLabelSize(0.075)
hdiv.GetYaxis().SetTitle("weightsdn/weights")
hdiv.GetYaxis().SetTitleSize(0.11)
hdiv.GetYaxis().SetTitleOffset(.45)
hdiv.GetYaxis().SetLabelSize(0.08)
hdiv.GetYaxis().SetLabelOffset(0.02)
hdiv.GetYaxis().SetNdivisions(503)
hdiv.SetMinimum(0.)
hdiv.SetMarkerStyle(8)
hdiv.SetMaximum(2.)


leg = ROOT.TLegend(0.60,0.65,0.88,0.88)
leg.AddEntry(hnom," normed weights","l")
leg.AddEntry(htest,"normed hack-dn weights","l")


tc = ROOT.TCanvas("tc","tc",700,900)
p1 = ROOT.TPad("pvis","pvis",0.0,0.3,1.0,1.0)
p2 = ROOT.TPad("prat","prat",0.0,0.0,1.0,0.3)

lab = ROOT.TPaveText(0.4,0.2,0.88,0.60,"NBNDC")
lab.AddText(signal.split("UF0_")[-1].split("_Tune")[0])
lab.AddText("Normed weights integral: "+str(round(hnom.Integral(),2)))
lab.AddText("    norm PU-dn integral: "+str(round(htest.Integral(),2)))
lab.SetFillColor(0)


tc.cd()
p1.Draw()
p1.cd()
hnom.Draw("HIST")
htest.Draw("HISTSAME")
leg.Draw()
lab.Draw()
p1.Update()
tc.cd()
p2.Draw()
p2.cd()
hdiv.Draw()



pngname = go.makeOutFile(signal.replace(".","_"),"pileupweightcomp",".png","100","300","75","8E-10")
tc.SaveAs(pngname)



