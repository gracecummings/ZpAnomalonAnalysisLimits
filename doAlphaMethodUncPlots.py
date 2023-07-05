import ROOT
import gecorg_test as go
import numpy as np

def makeBinLowEdges(hist,lastnormbin):
    #this takes the histogram ou want to rebin, and you just give it the last normal bin
    nbins = hist.GetNbinsX()
    binw  = hist.GetBinWidth(1)
    iniedges = [hist.GetBinLowEdge(i) for i in range(nbins+2)]#+2 for the actual highedge
    if lastnormbin not in iniedges:
        print("desired bin edge not possible, try again")
    newedges = [edg for edg in iniedges if (edg <= lastnormbin)]
    newedges.append(iniedges[-1])
    newedges = newedges[1:]
    return np.array(newedges)


def getDeviatedOverNominal(hist,histnom):
    interest = hist.Integral()
    intnom = histnom.Integral()
    devonom = interest/intnom
    return devonom

def newNameAndStructure(hist,name,rebindiv,limrangelow,limrangehigh):
    hist.Rebin(rebindiv)
    nbins = hist.GetNbinsX()
    binw  = hist.GetBinWidth(1)
    newbins = [limrangelow+x*binw for x in range(int((limrangehigh-limrangelow)/binw))]
    nh = ROOT.TH1F(name,name,len(newbins),limrangelow,limrangehigh)
    for b,le in enumerate(newbins):
        bnum = hist.FindBin(le)
        bincontent = hist.GetBinContent(bnum)
        binerror   = hist.GetBinError(bnum)
        nh.SetBinContent(b+1,bincontent)
        nh.SetBinError(b+1,binerror)
    return nh

def getDeviatedOverNominalSummary(histup,hist,histdwn,name):
    upnum  = getDeviatedOverNominal(histup,hist)
    dwnnum = getDeviatedOverNominal(histdwn,hist)
    up = "Deviated over nominal  for "+name+" up variation: "+str(round(upnum,3))
    dwn = "Devoated over nominal for "+name+" dwn variation: "+str(round(dwnnum,3))
    return up,dwn

def makeDeviatedTPave(histup,hist,histdwn,name):
    updiv  = getDeviatedOverNominal(histup,hist)
    dwndiv = getDeviatedOverNominal(histdwn,hist)
    updiv = "up deviated over nominal: "+str(round(updiv,3))
    dwndiv = "dwn deviated over nominal: "+str(round(dwndiv,3))
    lab = ROOT.TPaveText(.5,.5,.85,.65,"NBNDC")
    lab.AddText(updiv)
    lab.AddText(dwndiv)
    lab.SetFillColor(0)
    return lab


if __name__=='__main__':

    limrangelow = 1800
    limrangehigh = 10000

    #dytf = ROOT.TFile('analysis_output_ZpAnomalon/2022-06-23/Run2_161718_dy_extraploationsystnominal_kfnom_btagnom_muidnom_elidnom_elreconom__Zptcut100.0_Hptcut300.0_metcut0.0_btagwp0.8.root')
    #dytf = ROOT.TFile('mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    #dytf = ROOT.TFile('analysis_output_ZpAnomalon/2023-01-23/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    #dytf = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-08/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    dytf = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-09/alpha_method_ttbar_likli/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    keys = dytf.GetListOfKeys()
    keys = [key.GetName() for key in keys]
    unckeys = [key for key in keys if "extrap_" in key]
    #print(unckeys)
    unckeysup = [key for key in unckeys if "Up" in key]
    unckeysdn = [key for key in unckeys if "Down" in key]
    unckeysup = sorted(unckeysup)
    unckeysdn = sorted(unckeysdn)
    
    nomdyr = dytf.Get('extrap')
    nomdy = newNameAndStructure(nomdyr,"DY",1,limrangelow,limrangehigh)
    newb = makeBinLowEdges(nomdy,2800)

    #print("The bin edges: ",newb)
    nomdy = nomdy.Rebin(len(newb)-1,"DY",newb)
    #nomdy.SetLineColor(ROOT.kBlack)
    #nomdy.GetYaxis().SetRangeUser(0,6)
    #nomdy.SetStats(0)

    rateuncdict = {}
    shapeuncdict = {}
    histlist = []
    
    for i in range(len(unckeysup)):
        applist = [0,1.0,0,0]#Setup to only apply to DY
        upkey = unckeysup[i]
        print(upkey)
        dnkey = unckeysdn[i]
        upr = dytf.Get(upkey)
        dnr = dytf.Get(dnkey)
        up = newNameAndStructure(upr,"DY_"+upkey,1,limrangelow,limrangehigh)
        dn = newNameAndStructure(dnr,"DY_"+dnkey,1,limrangelow,limrangehigh)
        up = up.Rebin(len(newb)-1,"DY_"+upkey,newb)
        dn = dn.Rebin(len(newb)-1,"DY_"+dnkey,newb)
        
        uprate = getDeviatedOverNominal(up,nomdy)
        dwnrate = getDeviatedOverNominal(dn,nomdy)
        genname = upkey.split("Up")[0]
        ratestr = "-"
            
        if abs(round(1-uprate,2))-abs(round(1-dwnrate,2)) == 0:
            ratestr = str(max(round(uprate,2),round(dwnrate,2)))
        elif "vv" not in genname:
             ratestr = str(round(dwnrate,2))+"/"+str(round(uprate,2))
        else:
            ratestr = str(round(dwnrate,3))+"/"+str(round(uprate,3))
        rateuncdict[genname] = {"type":"lnN","proc":["-",ratestr,"-","-"]}

        #now the shapes
        shapeuncdict[genname] = {"type":"shape","unc":1.0,"proc":[0,1,0,0]}
        histlist.append(up)
        histlist.append(dn)
        
        up.SetLineColor(ROOT.kRed)
        dn.SetLineColor(ROOT.kBlue)

        upOvnom,dnOvnom = getDeviatedOverNominalSummary(up,nomdy,dn,upkey.split("Up")[0])
        print(upOvnom)
        print(dnOvnom)
        devilab = makeDeviatedTPave(up,nomdy,dn,upkey.split("Up")[0])
        

        leg = ROOT.TLegend(0.45,0.7,0.85,0.85)
        leg.SetBorderSize(0)
        leg.AddEntry(nomdy,"nominal DY extrap","l")
        leg.AddEntry(up,upkey,"l")
        leg.AddEntry(dn,dnkey,"l")
        
        tc = ROOT.TCanvas("tc"+str(i),"tc"+str(i),500,400)
        tc.Draw()
        nomdy.Draw('hist')
        up.Draw('histsame')
        dn.Draw('histsame')
        leg.Draw()
        devilab.Draw()
        tc.Update()
        tc.SaveAs(go.makeOutFile('Run2_161718_ZllHbbMET_dySystematicsComp',upkey.split('Up')[0],'.png','100','300','75','8E-10'))
    
    print(rateuncdict)
    print(shapeuncdict)
    print(histlist)


