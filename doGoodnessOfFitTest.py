import ROOT
import gecorg_test as go
import tdrstyle
import CMS_lumi
import numpy as np

if __name__=='__main__':
    #Gather info
    #datf = ROOT.TFile('analysis_output_ZpAnomalon/2022-07-04/higgsCombineZp4000ND800NS200_GoFDataSB.GoodnessOfFit.mH120.root')#4.571
    #datf = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-14/higgsCombine_ZllHbbMET_dataObsFromSB_mumu_Zp4000ND800NS200_notoys.GoodnessOfFit.mH120.root')#1.74615
    #datf = ROOT.TFile('analysis_output_ZpAnomalon/2023-06-28/higgsCombinegof_unblind_mumu_Zp2000ND400NS200.GoodnessOfFit.mH120.root')
    ##datf = ROOT.TFile('analysis_output_ZpAnomalon/2023-06-30/higgsCombineunblind_mumu_Zp4000ND800NS200_KS_GoF_norranges_notoysfrequentist_notoysstraightmeasure.GoodnessOfFit.mH120.root')
    datf = ROOT.TFile('analysis_output_ZpAnomalon/2023-06-30/higgsCombineunblind_mumu_Zp4000ND800NS200_saturdated_GoF_norranges_straightmasure.GoodnessOfFit.mH120.root')
    tdat = datf.Get("limit")
    tdat.GetEntry(0)
    datstat = tdat.limit
    #mcf = ROOT.TFile('analysis_output_ZpAnomalon/2023-06-28/higgsCombinegof_unblind_mumu_Zp2000ND400NS200.GoodnessOfFit.mH120.123456.root')
    #mcf = ROOT.TFile('analysis_output_ZpAnomalon/2022-07-04/higgsCombineTest.GoodnessOfFit.mH120.123456.root')
    #mcf = ROOT.TFile('analysis_output_ZpAnomalon/2022-07-04/higgsCombineZp4000ND800NS200GoFMCtoys_rangeSet.GoodnessOfFit.mH120.123456.root')
    #mcf = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-14/higgsCombine_ZllHbbMET_dataObsFromSB_mumu_Zp4000ND800NS200_toysfreq_12345seed_lowerrange-100.GoodnessOfFit.mH120.12345.root')
    #mcf = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-14/higgsCombine_ZllHbbMET_dataObsFromSB_mumu_Zp4000ND800NS200_toysfreq_12345seed_lowerrange-100_bypass.GoodnessOfFit.mH120.12345.root')
    #mcf = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-14/higgsCombine_ZllHbbMET_dataObsFromSB_mumu_Zp4000ND800NS200_toysfreq_12345seed.GoodnessOfFit.mH120.12345.root')
    ##mcf = ROOT.TFile('analysis_output_ZpAnomalon/2023-06-30/higgsCombineunblind_mumu_Zp4000ND800NS200_KS_GoF_norranges_notoysfrequentist.GoodnessOfFit.mH120.123456.root')
    mcf = ROOT.TFile('analysis_output_ZpAnomalon/2023-06-30/higgsCombineunblind_mumu_Zp4000ND800NS200_saturated_GoF_norranges_toysfrequentist.GoodnessOfFit.mH120.123456.root')
    
    tmc = mcf.Get("limit")
    toyn = tmc.GetEntries()
    

    #Define Hists/things to plot
    hmc = ROOT.TH1F("hmc","Saturated test statistic "+str(toyn)+" toys",30,0,30)
    tstatl = ROOT.TLine(datstat,0,datstat,105.)
    tstatl.SetLineColor(ROOT.kRed)
    tstatl.SetLineWidth(5)
    lab = ROOT.TPaveText(0.15,0.2,0.35,0.3,"NBNDC")
    lab.AddText("Test Statistic From ")
    lab.AddText("Data fit = "+str(round(datstat,2)))
    lab.SetTextColor(ROOT.kRed)
    lab.SetFillColor(0)

    #Plots
    tmc.Draw("limit>>hmc","","")

    #Calculate P-value
    pval = hmc.Integral(hmc.FindBin(datstat),hmc.GetNbinsX())
    print(pval)
    pvalnorm = pval/hmc.Integral()
    print(pvalnorm)

    #integralabove =
    tpave = ROOT.TPaveText(0.5,0.5,0.7,0.7,"NBNDC")
    tpave.AddText("pval = "+str(pvalnorm))
    tpave.SetFillColor(0)
    #Draw
    tc = ROOT.TCanvas("tc","GoF Test",800,600)
    #tc.SetLogy()
    tc.cd()
    hmc.Draw()
    hmc.GetXaxis().SetTitle("Saturated Test Statistic")
    hmc.GetYaxis().SetTitle("toys per bin width 1")
    tstatl.Draw()
    lab.Draw()
    tpave.Draw()
    tc.Update()
    
    outf = go.makeOutFile("Run2_161718_gof_test","unblinded_ntoys"+str(toyn)+"_toysfreqToys_norangeincluded_saturated_0to300",".png","100","300","75","8E-10")
    tc.SaveAs(outf)
    
