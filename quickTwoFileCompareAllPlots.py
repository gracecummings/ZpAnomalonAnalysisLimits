import gecorg_test as go
import pandas as pd
import ROOT

fnom = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-09/alpha_method_ttbar_likli/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')

#fnom = ROOT.TFile('analysis_output_ZpAnomalon/pumapupwithnorm/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
ftest = ROOT.TFile('analysis_output_ZpAnomalon/pumapdwnwithnorm/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
#fnom = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-31/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_puw__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
#ftest = ROOT.TFile('analysis_output_ZpAnomalon/2023-06-01/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')#this is the normalized weight one


keys = fnom.GetListOfKeys()

for key in keys:
    htestname = key.GetName()

    if "par" in htestname:
        continue
    hnom = fnom.Get(htestname)
    htest = ftest.Get(htestname)
    
    #If the files are upouts
    #nomerrs = pd.read_pickle(fnomerrs)
    #testerrs = pd.read_pickle(ftesterrs)
    
    #for ibin in range(hnom.GetNbinsX()+1):
    #    if ibin == 0:
    #        continue
    #    else:
    #        hnom.SetBinError(ibin,nomerrs[htestname][ibin-1])
    #        htest.SetBinError(ibin,testerrs[htestname][ibin-1])
    
    ROOT.gStyle.SetOptStat(0)

    print(htestname)
    #print("   Nominal integral:         ",hnom.Integral(hnom.FindBin(1800),hnom.FindBin(10000)))
    #print("   Pileup weighted integral: ",htest.Integral(htest.FindBin(1800),htest.FindBin(10000)))
    #print("   perc diff:                ",(hnom.Integral(hnom.FindBin(1800),hnom.FindBin(10000))-htest.Integral(htest.FindBin(1800),htest.FindBin(10000)))/hnom.Integral(hnom.FindBin(1800),hnom.FindBin(10000))*100)
    #print(" doing the rebinning")

    hnom = go.newNameAndStructure(hnom,htestname,1,1800,10000)
    htest = go.newNameAndStructure(htest,htestname,1,1800,10000)
    newbinedges = go.makeBinLowEdges(hnom,2800)

    hnom = hnom.Rebin(len(newbinedges)-1,htestname,newbinedges)
    htest = htest.Rebin(len(newbinedges)-1,htestname,newbinedges)
    
    
    print("   Nominal integral:         ",hnom.Integral())
    print("   Pileup dwn-weighted integral: ",htest.Integral())
    print("   perc diff:                ",((hnom.Integral()-htest.Integral())/hnom.Integral())*100)
    
    hdiv = htest.Clone()
    hdiv.Divide(htest,hnom)
    
    
    hdiv.GetXaxis().SetTitleSize(0.07)
    hdiv.GetXaxis().SetTitleOffset(0.65)
    hdiv.GetXaxis().SetLabelSize(0.075)
    hdiv.GetYaxis().SetTitle("weightsdwn/noweights")
    hdiv.GetYaxis().SetTitleSize(0.11)
    hdiv.GetYaxis().SetTitleOffset(.45)
    hdiv.GetYaxis().SetLabelSize(0.08)
    hdiv.GetYaxis().SetLabelOffset(0.02)
    hdiv.GetYaxis().SetNdivisions(503)
    hdiv.SetMinimum(0.)
    hdiv.SetMarkerStyle(8)
    hdiv.SetMaximum(2.)
    
    hnom.GetYaxis().SetTitle("Events")
    hnom.GetXaxis().SetTitle("RJR Z Prime Estimator (GeV)")
    htest.SetLineColor(ROOT.kRed)
    
    leg = ROOT.TLegend(0.60,0.65,0.88,0.88)
    leg.AddEntry(hnom,"no weights","l")
    leg.AddEntry(htest,"mapdwn/norm weights","l")
    lab = ROOT.TPaveText(0.60,0.4,0.88,0.55,"NBNDC")
    lab.AddText(htestname)
    lab.AddText("No PU integral: "+str(round(hnom.Integral(),2)))
    lab.AddText("PUwtdwn  integral: "+str(round(htest.Integral(),2)))
    lab.AddText("  PUwtdwn/nominal: "+str(round(htest.Integral()/hnom.Integral(),2)))
    lab.SetFillColor(0)
    
    
    tc = ROOT.TCanvas("tc","tc",700,900)
    p1 = ROOT.TPad("pvis","pvis",0.0,0.3,1.0,1.0)
    p2 = ROOT.TPad("prat","prat",0.0,0.0,1.0,0.3)
    
    
    tc.cd()
    p1.Draw()
    p1.cd()
    hnom.Draw("HIST")
    htest.Draw("HISTSAME")
    leg.Draw()
    lab.Draw()
    p1.Update()
    tc.cd()
    tc.Update()
    p2.Draw()
    p2.cd()
    hdiv.Draw("pe")
    tc.Update()
    
    
    pngname = go.makeOutFile("alphar_"+htestname,"pileupweightcomp",".png","100","300","75","8E-10")
    tc.SaveAs(pngname)
    


