import ROOT
import glob
import gecorg_test as go
import pandas as pd
import tdrstyle
import CMS_lumi

tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = "137.6 fb^{-1}"
CMS_lumi.writeExtraText = 1



if __name__=='__main__':

    zptcut  = '100.0'
    hptcut  = '300.0'
    metcut  = '75.0'
    btagwp  = '0.8'
    chan    = 'mumu'
    sigxs   = 1.0

    #sigofficial = glob.glob("mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref/Autumn18.ZpAnomalonHZ_UFO-Zp*NS200*upout_signalr*.root")
    sigofficial = go.signal_run2("mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref","100.0","300.0","75.0","0.8",1,[16,17,18],"systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco")
    #sigofficial = go.signal_run2("mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref","100.0","300.0","75.0","0.8")
    #siginterp   = ["analysis_output_ZpAnomalon/2022-10-05/interpolation_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root"]
    siginterp   = ["analysis_output_ZpAnomalon/2023-05-16/interpolation_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root"]

    #offdict = {"mzp":[],"mnd":[],"mns":[],"isoff":[],"tfile":[],"hist":[]}
    offdict = {"mzp":[],"mnd":[],"mns":[],"isoff":[],"hist":[]}

    tftest = ROOT.TFile(sigofficial.sig18sr[0])
    empty  = tftest.Get("h_zp_jigm")
    empty.Reset("ICESM")

    for fstr in sigofficial.sigsbyname:
        #signame = go.nameSignal(fstr)
        #signame = fstr.replace("-","")
        signame = fstr
        mzp,mnd,mns = go.massPoints(signame)
        #tf = ROOT.TFile(fstr)
        offdict["mzp"].append(mzp)
        offdict["mnd"].append(mnd)
        offdict["mns"].append(mns)
        offdict["isoff"].append(True)
        offdict["hist"].append(sigofficial.getAddedHist(signame,1,empty.Clone(),"sr","h_zp_jigm"))
        #offdict["tfile"].append(tf)
        #offdict["hist"].append("h_zp_jigm")

    siginterpf = ROOT.TFile(siginterp[0])
    interpkeys = siginterpf.GetListOfKeys()

    for key in interpkeys:
        name = key.GetName()
        mzp = name.split("Zp")[-1].split("ND")[0]
        mnd = name.split("ND")[-1].split("NS")[0]
        mns = name.split("NS")[-1]
        offdict["mzp"].append(int(mzp))
        offdict["mnd"].append(int(mnd))
        offdict["mns"].append(int(mns))
        offdict["isoff"].append(False)
        #offdict["tfile"].append(siginterpf)
        #offdict["hist"].append(name)
        offdict["hist"].append(siginterpf.Get(name))                       

    df = pd.DataFrame.from_dict(offdict)

    #mzps = set(df["mzp"])
    #mnds = set(df["mnd"])
    mzps = sorted(list(set(df["mzp"])))
    mnds = sorted(list(set(df["mnd"])))

    #mzps = list(df["mzp"])
    #mnds = list(df["mnd"])

    print(mzps)
    print(mnds)

    print(df)
    df = df.sort_values(by=["mzp","mnd"])
    print(df)

    for mass in mnds:
        CMS_lumi.extraText = "SR, m_{ND} = "+str(mass)
        dfint = df[df["mnd"] == mass]
        sigcolors = go.colsFromPalette(dfint,ROOT.kCMYK)
        tc = ROOT.TCanvas("tc"+str(mass),"tc"+str(mass),800,600)
        leg = ROOT.TLegend(0.65,0.45,0.88,0.75)
        leg.SetBorderSize(0)
        tc.cd()
        i = 0
        for index,row in dfint.iterrows():
            #h = row["tfile"].Get(row["hist"])
            h = row["hist"]
            h.SetLineColor(sigcolors[i])
            h.SetMaximum(3)
            h.GetXaxis().SetRangeUser(0,10000)
            h.SetStats(0)
            h.GetXaxis().SetTitle("RJR Z' Mass Estimator (GeV)")
            #h.GetYaxis().SetTitle("Events (unscaled) / 200 GeV")
            h.GetYaxis().SetTitle("Events / 200 GeV")
            h.GetYaxis().SetTitleOffset(1.4)
            leg.AddEntry(h,"Zp"+str(row["mzp"])+" ND"+str(row["mnd"])+" NS"+str(row["mns"]),"l")
            #if row["isoff"]:
            #h.Rebin(2)
            #else:
            if not row["isoff"]:
                h.SetLineWidth(2)
            #    h.SetLineStyle(9)
            if i == 0:
                h.Draw("hist")
            else:
                h.Draw("histsame")
            i+=1
        leg.Draw()
        CMS_lumi.CMS_lumi(tc,4,13)
                

        pngname = go.makeOutFile("signal_interp_plots",'MND'+str(mass),'.png',zptcut,hptcut,metcut,btagwp)
        tc.SaveAs(pngname)

    for mass in mzps:
        CMS_lumi.extraText = "SR, m_{Z'} = "+str(mass)
        dfint = df[df["mzp"] == mass]
        sigcolors = go.colsFromPalette(dfint,ROOT.kCMYK)
        tc = ROOT.TCanvas("tc"+str(mass),"tc"+str(mass),800,600)
        leg = ROOT.TLegend(0.65,0.45,0.88,0.75)
        leg.SetBorderSize(0)
        tc.cd()
        i = 0
        for index,row in dfint.iterrows():
            if int(row["mnd"]) > 1800:
                continue
            #h = row["tfile"].Get(row["hist"])
            h = row["hist"]

            print("mzp ",mass)
            print("mnd ",int(row["mnd"]))
            print("   integral: ",h.Integral(0,26))
            h.SetLineColor(sigcolors[i])
            #h.SetMaximum(1500)
            h.SetMaximum(3)
            h.GetXaxis().SetRangeUser(0,10000)
            h.SetStats(0)
            h.GetXaxis().SetTitle("RJR Z' Mass Estimator (GeV)")
            h.GetYaxis().SetTitle("Events (unscaled) / 200 GeV")
            h.GetYaxis().SetTitleOffset(1.4)
            leg.AddEntry(h,"Zp"+str(row["mzp"])+" ND"+str(row["mnd"])+" NS"+str(row["mns"]),"l")
            if not row["isoff"]:
                h.SetLineWidth(2)
            if i == 0:
                h.Draw("hist")
            else:
                h.Draw("histsame")
            i+=1
        leg.Draw()
        CMS_lumi.CMS_lumi(tc,4,13)
                

        pngname = go.makeOutFile("signal_interp_plots",'MZP'+str(mass),'.png',zptcut,hptcut,metcut,btagwp)
        tc.SaveAs(pngname)

        
    
        
        

    
