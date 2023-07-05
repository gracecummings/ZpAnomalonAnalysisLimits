import gecorg_test as go
import glob
import ROOT

#fdwns = glob.glob('mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref_jec/Run2_161718_dy_extraploationalphat_systjecdwn_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
#fnom = ROOT.TFile('mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
#fups = glob.glob('mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref_jec/Run2_161718_dy_extraploationalphat_systjecup_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')

fnom = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-09/alpha_method_ttbar_likli/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
fpu = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-31/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_puw__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')

#fdwns = glob.glob('mumu_2022-03-31_ProperREOIDSF_unclmetsyst/Run2_161718_dy_extraploationsystuncldwn_kfnom_btagnom_muidnom_elidnom_elreconom_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
#fnom = ROOT.TFile('mumu_2022-03-31_ProperREOIDSF/Run2_161718_dy_extraploationsystnominal_kfnom_btagnom_muidnom_elidnom_elreconom__Zptcut100.0_Hptcut300.0_metcut0.0_btagwp0.8.root')
#fups = glob.glob('mumu_2022-03-31_ProperREOIDSF_unclmetsyst/Run2_161718_dy_extraploationsystunclup_kfnom_btagnom_muidnom_elidnom_elreconom_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')


keys = fnom.GetListOfKeys()

#limlims = {"dyl":[30,250],"ttl":[30,400],"vvl":[30,250],"totalfit":[30,250],"flsb":[30,70],"fhsb":[150,250],"sbdatfit":[30,250],"totextrap":[30,250]}#lims for norm fits

limlims = {"dysbl":[1500,10000],"dysrl":[1500,10000],"ttsbl":[1500,10000],"vvsbl":[1500,10000],"datsbl":[1500,10000],"alpha":[1500,10000]}

#fdwns = sorted(fdwns)
#fups = sorted(fups)

#print(fdwns)
#print(fups)

#fsysts = list(zip(fdwns,fups))
canvases = []

ROOT.gROOT.SetBatch(ROOT.kTRUE)

for key in keys:
    hname = key.GetName()
    print(hname)
    if "extrap" in hname:
        continue

    upname = "pileupweights"
    up = fpu.Get(hname)
    nom = fnom.Get(hname)
    up.SetLineColor(ROOT.kRed)
    up.SetMarkerColor(ROOT.kRed)
    nom.SetLineColor(ROOT.kBlack)
    nom.SetMarkerColor(ROOT.kBlack)
    
    if "DYSR" in hname:
        nom.GetYaxis().SetRangeUser(0,4)
    if "dysrl" in hname:
        nom.GetYaxis().SetRangeUser(0,5)

    if "l" not in hname and ("fit" not in hname) and ("tot" not in hname) and ("fhsb" not in hname):
        intlabel = ROOT.TPaveText(.6,.3,.9,.45,"NBNDC")
        #intlabel.AddText("Dwn integral: "+str(round(dwn.Integral(),2)))    
        intlabel.AddText("Nom integral: "+str(round(nom.Integral(),2)))
        intlabel.AddText("PU W integral : "+str(round(up.Integral(),2)))
        intlabel.SetFillColor(0)
        
        intlabel2 = ROOT.TPaveText(.6,.15,.9,.25,"NBNDC")
        #intlabel2.AddText("Dwn over nominal: "+str(round(dwn.Integral()/nom.Integral(),2)))    
        intlabel2.AddText("PU W over nominal : "+str(round(up.Integral()/nom.Integral(),2)))
        intlabel2.SetFillColor(0)
                
    if ("l" in hname) or ("fit" in hname) or ("tot" in hname) or ("fhsb" in hname):
        #print("found a fit")
        intlabel = ROOT.TPaveText(.6,.35,.9,.50,"NBNDC")
        #intlabel.AddText("Dwn integral: "+str(round(dwn.Integral(limlims[hname][0],limlims[hname][1],0.0001),2)))    
        intlabel.AddText("Nom integral: "+str(round(nom.Integral(limlims[hname][0],limlims[hname][1],0.0001),2)))
        intlabel.AddText("PU W integral : "+str(round(up.Integral(limlims[hname][0],limlims[hname][1],0.0001),2)))
        intlabel.SetFillColor(0)
        
        intlabel2 = ROOT.TPaveText(.6,.2,.9,.3,"NBNDC")
        #intlabel2.AddText("Dwn over nominal: "+str(round(dwn.Integral(limlims[hname][0],limlims[hname][1],0.0001)/nom.Integral(limlims[hname][0],limlims[hname][1],0.0001),2)))    
        intlabel2.AddText("PU W over nominal : "+str(round(up.Integral(limlims[hname][0],limlims[hname][1],0.0001)/nom.Integral(limlims[hname][0],limlims[hname][1],0.0001),2)))
        intlabel2.SetFillColor(0)
        
        #intlabel.GetLine(0).SetTextColor(ROOT.kBlue)
        intlabel.GetLine(0).SetTextColor(ROOT.kBlack)
        intlabel.GetLine(1).SetTextColor(ROOT.kRed)
        

                
    tc = ROOT.TCanvas(upname+"_"+hname,upname+"_"+hname,700,600)
    tc.cd()
    if "alpha" in hname:
        tc.SetLogy()
        nom.GetYaxis().SetRangeUser(0.1,1)
        nom.Draw()
        up.Draw("same")
        #dwn.Draw("SAME")
        intlabel.Draw()
        intlabel2.Draw()
        tc.Update()
    outname = go.makeOutFile(upname+"_"+hname,"_debug_shapes_alphatest",".png","100","300","0","8E-1")
    tc.SaveAs(outname)
    canvases.append(tc)
        

###For multiple files
#for i,pair in enumerate(fsysts):
#    if i % 5 == 5:
#        print("on parameter ",i)
#    dwnname = pair[0].split('syst')[-1].split('wn_hem')[0]
#    upname  = pair[1].split('syst')[-1].split('up_hem')[0]

    #dwnname = pair[0].split('syst')[-1].split('wn_kf')[0]
    #upname  = pair[1].split('syst')[-1].split('up_kf')[0]
#    if upname not in dwnname:
#        print("UPS AND DOWNS DO NOT MATCH")
#        continue
#    else:
#        print(upname)
        #print(dwnname)
#        fup = ROOT.TFile(pair[1])
#        fdwn = ROOT.TFile(pair[0])
#        for key in keys:
#            hname = key.GetName()
#            print(hname)
#
#
#            if "extrap" in hname:
#                continue
#            
#            up = fup.Get(hname)
#            dwn = fdwn.Get(hname)
#            nom = fnom.Get(hname)
#            up.SetLineColor(ROOT.kRed)
#            up.SetMarkerColor(ROOT.kRed)
#            dwn.SetLineColor(ROOT.kBlue)
#            nom.SetLineColor(ROOT.kBlack)
#            nom.SetMarkerColor(ROOT.kBlack)
#
#            if "DYSR" in hname:
#                nom.GetYaxis().SetRangeUser(0,4)
#            if "dysrl" in hname:
#                nom.GetYaxis().SetRangeUser(0,5)
#
#            if "l" not in hname and ("fit" not in hname) and ("tot" not in hname) and ("fhsb" not in hname):
#                intlabel = ROOT.TPaveText(.6,.3,.9,.45,"NBNDC")
#                intlabel.AddText("Dwn integral: "+str(round(dwn.Integral(),2)))    
#                intlabel.AddText("Nom integral: "+str(round(nom.Integral(),2)))
#                intlabel.AddText("Up integral : "+str(round(up.Integral(),2)))
#                intlabel.SetFillColor(0)
#    
#                intlabel2 = ROOT.TPaveText(.6,.15,.9,.25,"NBNDC")
#                intlabel2.AddText("Dwn over nominal: "+str(round(dwn.Integral()/nom.Integral(),2)))    
#                intlabel2.AddText("Up over nominal : "+str(round(up.Integral()/nom.Integral(),2)))
#                intlabel2.SetFillColor(0)
#                
#            if ("l" in hname) or ("fit" in hname) or ("tot" in hname) or ("fhsb" in hname):
#                #print("found a fit")
#                intlabel = ROOT.TPaveText(.6,.35,.9,.50,"NBNDC")
#                intlabel.AddText("Dwn integral: "+str(round(dwn.Integral(limlims[hname][0],limlims[hname][1],0.0001),2)))    
#                intlabel.AddText("Nom integral: "+str(round(nom.Integral(limlims[hname][0],limlims[hname][1],0.0001),2)))
#                intlabel.AddText("Up integral : "+str(round(up.Integral(limlims[hname][0],limlims[hname][1],0.0001),2)))
#                intlabel.SetFillColor(0)
#    
#                intlabel2 = ROOT.TPaveText(.6,.2,.9,.3,"NBNDC")
#                intlabel2.AddText("Dwn over nominal: "+str(round(dwn.Integral(limlims[hname][0],limlims[hname][1],0.0001)/nom.Integral(limlims[hname][0],limlims[hname][1],0.0001),2)))    
#                intlabel2.AddText("Up over nominal : "+str(round(up.Integral(limlims[hname][0],limlims[hname][1],0.0001)/nom.Integral(limlims[hname][0],limlims[hname][1],0.0001),2)))
#                intlabel2.SetFillColor(0)
#
#                intlabel.GetLine(0).SetTextColor(ROOT.kBlue)
#                intlabel.GetLine(1).SetTextColor(ROOT.kBlack)
#                intlabel.GetLine(2).SetTextColor(ROOT.kRed)
#    
#
#                
#            tc = ROOT.TCanvas(upname+"_"+hname,upname+"_"+hname,700,600)
#            tc.cd()
#            if "alpha" in hname:
#                tc.SetLogy()
#                nom.GetYaxis().SetRangeUser(0.1,1)
#            nom.Draw()
#            up.Draw("same")
#            dwn.Draw("SAME")
#            intlabel.Draw()
#            intlabel2.Draw()
#            tc.Update()
#            outname = go.makeOutFile(upname+"_"+hname,"_debug_shapes_alphatest",".png","100","300","0","8E-1")
#            tc.SaveAs(outname)
#            canvases.append(tc)
#            
#

#outf = ROOT.TFile(outname,"recreate")
#for can in canvases:
#    can.Write()
#outf.Close()
#fup.Close()
#fdwn.Close()
#fnom.Close()
