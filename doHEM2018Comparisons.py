import gecorg_test as go
import ROOT
import tdrstyle
import CMS_lumi

def applyStatsUncToSignal(hist,errseries):
    for ibin in range(hist.GetNbinsX()+1):
        if ibin == 0:
            continue
        else:
            binerr = errseries[ibin-1]
            hist.SetBinError(ibin,binerr)
    return hist


def gatherIndividualPlots(bkgs,tf,hname,bkgname,reg,years):
    hb = tf.Get(hname)
    hbempt = hb.Clone()
    hbempt.Reset("ICECM")
    hbkg = bkgs.getAddedHist(hbempt,bkgname,reg,hname,years=years)
    print(hbkg)
    return hbkg

def regionFormatter(regionstr):
    regdict = {'sideband':'sb','signalr':'sr','totalr':'tr'}
    formatted = regdict[regionstr]
    return formatted

def lumiFormatter(yearlist):
    lumidict = {16:35.9,17:41.53,18:59.74}
    lumi = 0
    for year in yearlist:
        lumi += lumidict[year]

    lumi = round(lumi,2)
    lumistr = str(lumi)+" fb^{-1}"
    return lumistr

def yearFormatter(yearlist):
    yearstr =''
    for year in yearlist:
        yearstr = yearstr+str(year)
    return yearstr

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

def makePlottableZprime(hzpnomsr):
    hzpnomsr.GetXaxis().SetTitle(titles['h_zp_jigm'][0])
    hzpnomsr.GetXaxis().SetTitleSize(0.05)
    hzpnomsr.GetXaxis().SetTitleOffset(1.1)
    hzpnomsr.GetXaxis().SetLabelSize(0.04)
    hzpnomsr.GetYaxis().SetTitle("Events")
    hzpnomsr.GetYaxis().SetTitleSize(0.05)
    hzpnomsr.GetYaxis().SetLabelSize(0.04)

def makePlottable(hzpnomsr,xaxtit):
    hzpnomsr.GetXaxis().SetTitle(xaxtit)
    hzpnomsr.GetXaxis().SetTitleSize(0.05)
    hzpnomsr.GetXaxis().SetTitleOffset(1.1)
    hzpnomsr.GetXaxis().SetLabelSize(0.04)
    hzpnomsr.GetYaxis().SetTitle("Events")
    hzpnomsr.GetYaxis().SetTitleSize(0.05)
    hzpnomsr.GetYaxis().SetLabelSize(0.04)

def getDeviatedOverNominal(hist,histnom):
    interest = hist.Integral()
    intnom = histnom.Integral()
    devonom = interest/intnom
    return devonom

def makeDeviatedTPave(hist,histup,title):
    updiv  = getDeviatedOverNominal(histup,hist)
    updiv = "Deviated over nominal: "+str(round(updiv,3))
    lab = ROOT.TPaveText(.4,.73,.9,.85,"NBNDC")
    lab.AddText(title)
    lab.AddText(updiv)
    lab.SetFillColor(0)
    return lab

def makeLegend(hist,histdev,samp):
    leg = ROOT.TLegend(0.6,0.4,0.9,0.6)
    leg.AddEntry(hist,samp+" nominal","l")
    leg.AddEntry(histdev,samp+" HEM","l")
    leg.SetFillColor(0)
    return leg

def getErrorOnIntegral(hist):
    sumsqrs = 0
    for i in range(hist.GetNbinsX()+1):
        err = hist.GetBinError(i)
        sumsqrs += err*err
    errint = sumsqrs**(1/2)
    return errint

if __name__=='__main__':
    sig_xsec  = '10'
    zptcut    = '100.0'
    hptcut    = '300.0'
    metcut    = '75.0'
    btagwp    = '0.8'
    year      = '18'
    regname   = 'sideband'
    systr     = 'systnominal_kfnom_btagnom_muidnom_mutrignom_elidnom_elreconom'
    systremu  = 'systnominal_kfnom_btagnom_muidnom_mutrignom_eltrignom_elidnom_elreconom'
    plot_data = 'True'
    sigdivsor = 5
    valid = False

    #Select Plotting years and region
    years = [16,17,18]
    if year:
        years = [int(year)]
    yearstr = yearFormatter(years)
    reg = regionFormatter(regname)
    if plot_data:
        print('Plotting the data!')
    else:
        print('You know you are not plotting data, right?')

    bkgsHEM  = go.backgrounds('mumu_2022-07-06_ProperMuIDTriSF_METXY_HEM_2018',zptcut,hptcut,metcut,btagwp,systr)
    bkgsnom  = go.backgrounds('mumu_2022-07-06_ProperMuIDTriSF_METXY',zptcut,hptcut,metcut,btagwp,systr)
    
    bkgsHEMeujec  = go.backgrounds('emu_2022-07-06_ProperSF_ProperEE_METXY-HEM_elec_rejec_jetHEM',zptcut,hptcut,metcut,btagwp,'systnominal_hemnom_kfnom_btagnom_muidnom_mutrignom_eltrignom_elidnom_elreconom')
    bkgsHEMeuveto  = go.backgrounds('emu_2022-07-06_ProperSF_ProperEE_METXY-HEM_elec_rejec',zptcut,hptcut,metcut,btagwp,'systnominal_hemnom_kfnom_btagnom_muidnom_mutrignom_eltrignom_elidnom_elreconom')
    bkgsnomeu  = go.backgrounds('emu_2022-07-06_ProperSF_ProperEE_METXY-NOHEM',zptcut,hptcut,metcut,btagwp,systremu)

    signom = go.signal('mumu_2022-07-06_ProperMuIDTriSF_METXY',zptcut,hptcut,metcut,btagwp,'10.0',[18],systr)
    sigHEM  = go.signal('mumu_2022-07-06_ProperMuIDTriSF_METXY_HEM_2018',zptcut,hptcut,metcut,btagwp,'10.0',[18],systr)

    signominfo = signom.getPreppedSig('sr',10.0,[18])
    sigHEMinfo = sigHEM.getPreppedSig('sr',10.0,[18])
    signominfo = sorted(signominfo,key = lambda sig: (sig["mzp"],sig["mnd"]))
    sigHEMinfo = sorted(sigHEMinfo,key = lambda sig: (sig["mzp"],sig["mnd"]))
    #bkgsnomAlphaSB  = go.backgrounds('mumu_2022-07-06_ProperMuIDTriSF_METXY_alphatest',zptcut,hptcut,metcut,btagwp,systr)
    #bkgsHEMAlphaSB  = go.backgrounds('mumu_2022-07-06_ProperMuIDTriSF_METXY_HEM_2018_alphatest',zptcut,hptcut,metcut,btagwp,systr)
    #data  = go.run2(pathplots,zptcut,hptcut,metcut,btagwp,'systnominal_btagnom_muidnom')
    #sigs =  go.signal(pathplots,zptcut,hptcut,metcut,btagwp,sig_xsec,years,systr)
    datareg = 'sb'

    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = lumiFormatter(years)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Simulation Preliminary"
    stkpadydims = [0.0,1.]
    ratpadydims = [0.0,0.0]
    tcanvasdims = [600,560]

    #Gather plots
    testyear = years[0]#picks first year in list, so desired year if only one
    testfile = bkgsnom.bkgs["DYJetsToLL"][testyear][reg][0][0]#stacked plots should always have DY
    testtfile = ROOT.TFile(testfile)
    #siginfo = sigs.getPreppedSig(reg,sig_xsec,years)
    #sigcolors = go.colsFromPalette(siginfo,ROOT.kCMYK)
    #siginfo = sorted(siginfo,key = lambda sig: (sig["mzp"],sig["mnd"])) 

    titles = {
        "h_z_pt":["Z p_{T} (GeV)",0,40,2],
        "h_z_eta":["\eta_{Z}",0,100,1],
        "h_z_phi":["\phi_{Z}",0,90,2],
        "h_z_phiw":["\phi_{Z}",0,90,2],
        "h_z_m":["m_{Z} (GeV)",0,100,1],
        "h_h_pt":["Higgs p_{T} (GeV)",0,60,1],
        "h_h_eta":["\eta_{Higss}",0,130,1],
        "h_h_phi":["\phi_{Higgs}",0,70,2],
        "h_h_phiw":["\phi_{Higgs}",0,70,2],
        "h_h_m":["m_{h} (GeV)",0,50,1],
        "h_h_sd":["Higgs Soft Drop Mass (GeV)",0,70,1],#45 normally max
        "h_metxy":["p_{T}^{miss} xy corrected (GeV)",0,100,1],
        "h_metxy_phi":["\phi p_{T}^{miss} xy corrected",0,80,2],
        "h_metxy_phiw":["\phi p_{T}^{miss} xy corrected",0,80,2],
        "h_met":["p_{T}^{miss} (GeV)",0,100,1],
        "h_met_phi":["\phi p_{T}^{miss}",0,80,2],
        "h_met_phiw":["\phi p_{T}^{miss}",0,80,2],
        "h_zp_jigm":["Jigsaw Mass Estimator Z'",0,60,2],
        "h_nd_jigm":["Jigsaw Mass Estimator ND",0,60,1],
        "h_ns_jigm":["Jigsaw Mass Estimator NS",0,100,1],
        "h_btag":["btag operating point",0,70,1],
        "h_dphi_zh":["\Delta\phi_{ZH}",0,80,2],
        "h_dphi_zmet":["\Delta\phi_{ZMET}",0,60,2],
        "h_dphi_hmet":["\Delta\phi_{HMET}",0,60,2],
        "h_dr_zh":["\Delta R(ZH)",0,170,1],
        "h_dr_lmuh":["\Delta R(lmu,H)",0,200,1],
        "h_dr_slmuh":["\Delta R(slmu,H)",0,200,1],
        "h_dr_slmulmu":["\Delta R(slmu,lmu)",0,60,1],
        "h_dr_gz_gh":["\Delta R(genZ,genH)",0,250,1],
        "h_dr_lmu_gh":["\Delta R(lmu,genH)",0,250,1],
        "h_dr_slmu_gh":["\Delta R(slmu,genH)",0,250,1],
        "h_LMu_pt":["leading \mu p_{T} (GeV)",0,40,1],
        "h_LMu_phi":["\phi_{leading \mu}",0,100,2],
        "h_LMu_eta":["\eta_{leading \mu} ",0,100,1],
        "h_sLMu_pt":["subleading \mu p_{T} (GeV)",0,40,1],
        "h_sLMu_phi":["\phi_{subleading \mu}",0,100,2],
        "h_sLMu_eta":["\eta_{subleading \mu}",0,100,1],
    }

    hnomtotplots = []
    hHEMtotplots = []
    hnomemuplots = []
    hHEMemvetplots = []
    hHEMemuplots = []
    hnomsrplots = []
    hnomalphasbplots = []
    hHEMsrplots = []
    hHEMalphasbplots = []
    hists = ['h_z_pt','h_z_phi','h_z_eta']
    bkg = "TT"
    divs = []
    #for bkg in bkgnames:
    #    print('Looking at bkg ',bkg)
    #For stright ttbar MC background
    hzpnomsr  = gatherIndividualPlots(bkgsnom,testtfile,'h_zp_jigm',bkg,'sr',years)
    hzpHEMsr  = gatherIndividualPlots(bkgsHEM,testtfile,'h_zp_jigm',bkg,'sr',years)
    
    for hist in hists:
        hdict  = {}
        #mumu
        hzpnom = gatherIndividualPlots(bkgsnom,testtfile,hist,bkg,'tr',years)
        hzpnom.Rebin(titles[hist][3])
        hzpnom.SetLineColor(ROOT.kBlack)
        makePlottable(hzpnom,titles[hist][0])
        #emu
        hzpnomeu = gatherIndividualPlots(bkgsnomeu,testtfile,hist,bkg,'tr',years)
        hzpnomeu.Rebin(titles[hist][3])
        hzpnomeu.SetLineColor(ROOT.kBlack)
        makePlottable(hzpnomeu,titles[hist][0])

        #mumu
        hzpHEM  = gatherIndividualPlots(bkgsHEM,testtfile,hist,bkg,'tr',years)
        hzpHEM.Rebin(titles[hist][3])
        hzpHEM.SetLineColor(ROOT.kRed)
        makePlottable(hzpHEM,titles[hist][0])
    
        #emu electron veto
        hzpHEMveto  = gatherIndividualPlots(bkgsHEMeuveto,testtfile,hist,bkg,'tr',years)
        hzpHEMveto.Rebin(titles[hist][3])
        hzpHEMveto.SetLineColor(ROOT.kRed)
        makePlottable(hzpHEMveto,titles[hist][0])
        #emu electron veto and jet scaleing
        hzpHEMemu =  gatherIndividualPlots(bkgsHEMeujec,testtfile,hist,bkg,'tr',years)
        hzpHEMemu.Rebin(titles[hist][3])
        hzpHEMemu.SetLineColor(ROOT.kRed)
        makePlottable(hzpHEMemu,titles[hist][0])
        #hzpHEMsr  = gatherIndividualPlots(bkgsHEM,testtfile,'h_zp_jigm',bkg,'sr',years)

        #halphanomsb = gatherIndividualPlots(bkgsnomAlphaSB,testtfile,'h_zp_jigm',bkg,'sb',years)
        #halphaHEMsb = gatherIndividualPlots(bkgsHEMAlphaSB,testtfile,'h_zp_jigm',bkg,'sb',years)

        #ttbar region divs
        #div = hzpHEM.Clone()
        #div.Divide(hzpHEM,hzpnom)
        #div.SetMarkerStyle(8)
        #div.SetLineColor(ROOT.kBlack)
        #div.SetMarkerColor(ROOT.kBlack)
        #divs.append(div)


        hnomtotplots.append(hzpnom)
        #hnomalphasbplots.append(halphanomsb)
        hHEMtotplots.append(hzpHEM)
        #hHEMalphasbplots.append(halphaHEMsb)
        hnomemuplots.append(hzpnomeu)
        hHEMemvetplots.append(hzpHEMveto)
        hHEMemuplots.append(hzpHEMemu)


    signoms = []
    sigHEMs = []
    sigdivs = []
    for s,sig in enumerate(signominfo):
        print(sig["name"])
        hemidx = -999
        for z,zig in enumerate(sigHEMinfo):
            if sig["name"] not in zig["name"]:
                continue
            if sig["name"] in zig["name"]:
                print(zig["name"])
                hemidx = z
        if hemidx < 0:
            continue
        hsig = sig["tfile"].Get('h_zp_jigm')
        hHEM = sigHEMinfo[hemidx]["tfile"].Get('h_zp_jigm')
        hsig = applyStatsUncToSignal(hsig,sig['errdf']['h_zp_jigm'])
        hHEM = applyStatsUncToSignal(hHEM,sigHEMinfo[hemidx]['errdf']['h_zp_jigm'])
        hsig = newNameAndStructure(hsig,sig['name'],titles['h_zp_jigm'][3],1400,3000)
        hHEM = newNameAndStructure(hHEM,sig['name']+'HEM',titles['h_zp_jigm'][3],1400,3000)
        hsig.SetLineColor(ROOT.kBlack)
        hsig.SetLineColor(ROOT.kRed)
        makePlottableZprime(hsig)
        makePlottableZprime(hHEM)

        div = hHEM.Clone()
        div.Divide(hHEM,hsig)
        div.SetMarkerStyle(8)
        div.SetLineColor(ROOT.kBlack)
        div.SetMarkerColor(ROOT.kBlack)

        tcmc = ROOT.TCanvas("tcmc",'HEM Check',600,800)
        pdmc = ROOT.TPad('pdmc','ttmc',0.0,0.3,1.0,1.0)
        ratmc = ROOT.TPad('ratmc','rat',0.0,0.0,1.0,0.3)
        
        tcmc.Draw()
        tcmc.cd()
        pdmc.Draw()
        pdmc.cd()
        hHEM.Draw('hist,e1')
        hsig.Draw('histsame,e1')
        lab = makeDeviatedTPave(hsig,hHEM,sig["name"].split("_Tune")[0])
        lab.Draw()
        leg = makeLegend(hsig,hHEM,sig["name"].split("_Tune")[0])
        leg.Draw()
        tcmc.cd()
        ratmc.Draw()
        ratmc.cd()
        div.Draw()
        tcmc.cd()
        tcmc.Update()
        
        tcmc.SaveAs("HEM_Check_"+sig["name"]+".png")


        #signoms.append(hsig)
        #sigHEMs.append(hHEM)
        #sigdivs.append(div)


    
    #for the emu, only care about the relative yields. The Z kinematics will not change
    print(" The nominal mumu yields: ")
    print("                      pT: {0} +- {1}".format(hnomtotplots[0].Integral(),getErrorOnIntegral(hnomtotplots[0])))
    print("    The HEM  mumu yields: ")
    print("                      pT: {0} +- {1}".format(hHEMtotplots[0].Integral(),getErrorOnIntegral(hHEMtotplots[0])))
    print(" The deviated over nominal mumu yields:")
    print("                                    pT: ",hHEMtotplots[0].Integral()/hnomtotplots[0].Integral())


    print("            The nominal emu yields: ")
    print("                                pT: {0} +- {1}".format(hnomemuplots[0].Integral(),getErrorOnIntegral(hnomemuplots[0])))
    print(" The HEM emu yields, electron veto: ")
    print("                                pT: {0} +- {1}".format(hHEMemvetplots[0].Integral(),getErrorOnIntegral(hHEMemvetplots[0])))#these are actually just weighting, not vetoed
    print("  The HEM emu yields, veto AND jet: ")
    print("                                pT: {0} +- {1}".format(hHEMemuplots[0].Integral(),getErrorOnIntegral(hHEMemuplots[0])))
    print("              el veto over nominal: ",hHEMemvetplots[0].Integral()/hnomemuplots[0].Integral())
    print("         veto and jet over nominal: ",hHEMemuplots[0].Integral()/hnomemuplots[0].Integral())
    print("            veto and jet over veto: ",hHEMemuplots[0].Integral()/hHEMemvetplots[0].Integral())
    


    #tc3 = ROOT.TCanvas('tc3','HEM Check ttbar emu',1800,1600)
    #pd1 = ROOT.TPad('pd1','mumutotpt',

    #Make plots ttbar MC background in mumu 
    #hzpnomsr  = newNameAndStructure(hzpnomsr,"TT_nominal",titles['h_zp_jigm'][3],1400,3000)
    #hzpHEMsr  = newNameAndStructure(hzpHEMsr,"TT_HEM",titles['h_zp_jigm'][3],1400,3000)
    #hzpnomsr.SetLineColor(ROOT.kBlack)
    #hzpHEMsr.SetLineColor(ROOT.kRed)
    #makePlottableZprime(hzpnomsr)
    #makePlottableZprime(hzpHEMsr)

    #div = hzpHEMsr.Clone()
    #div.Divide(hzpHEMsr,hzpnomsr)
    #div.SetMarkerStyle(8)
    #div.SetLineColor(ROOT.kBlack)
    #div.SetMarkerColor(ROOT.kBlack)

    #Single ttbar can be co-opted for sig
    
    #tcmc = ROOT.TCanvas("tcmc",'HEM Check',600,800)
    #pdmc = ROOT.TPad('pdmc','ttmc',0.0,0.3,1.0,1.0)
    #ratmc = ROOT.TPad('ratmc','rat',0.0,0.0,1.0,0.3)

    #tcmc.Draw()
    #tcmc.cd()
    #pdmc.Draw()
    #pdmc.cd()
    #hzpHEMsr.Draw('hist,e1')
    #hzpnomsr.Draw('histsame,e1')
    #lab = makeDeviatedTPave(hzpnomsr,hzpHEMsr,"TT MC SR")
    #lab.Draw()
    #leg = makeLegend(hzpnomsr,hzpHEMsr,"TT MC SR")
    #leg.Draw()
    #tcmc.cd()
    #ratmc.Draw()
    #ratmc.cd()
    #div.Draw()
    #tcmc.cd()
    #tcmc.Update()

    #tcmc.SaveAs("HEM_Check_TT_MC_SR.png")



    #Makre Plots - Alpha EMthod
#    hdysrnom = hnomsrplots[0].Rebin(titles['h_zp_jigm'][3])
#    hdysrHEM = hHEMsrplots[0].Rebin(titles['h_zp_jigm'][3])
#    hdyASBnom = hnomalphasbplots[0].Rebin(titles['h_zp_jigm'][3])
#    hdyASBHEM = hHEMalphasbplots[0].Rebin(titles['h_zp_jigm'][3])
#    httASBnom = hnomalphasbplots[1].Rebin(titles['h_zp_jigm'][3])
#    httASBHEM = hHEMalphasbplots[1].Rebin(titles['h_zp_jigm'][3])
#    hwzASBnom = hnomalphasbplots[2].Rebin(titles['h_zp_jigm'][3])
#    hwzASBHEM = hHEMalphasbplots[2].Rebin(titles['h_zp_jigm'][3])
#    hzzASBnom = hnomalphasbplots[3].Rebin(titles['h_zp_jigm'][3])
#    hzzASBHEM = hHEMalphasbplots[3].Rebin(titles['h_zp_jigm'][3])
    #hvvASBnom = hwzASBnom.Clone()
    #hvvASBnom.Add(hzzASBnom)
    #hvvASBHEM = hwzASBHEM.Clone()
    #hvvASBHEM.Add(hzzASBHEM)
    #noms = [hdysrnom,hdyASBnom,httASBnom,hvvASBnom]
    #hems = [hdysrHEM,hdyASBHEM,httASBHEM,hvvASBHEM]

    #Diboson Plots
#    noms = hnomsrplots
#    hems = hHEMsrplots
#
#    #divs = []
#    #stacks = []
#    for h in range(len(noms)):
#        noms[h] = newNameAndStructure(noms[h],bkgnames[h]+"_nominal",titles['h_zp_jigm'][3],1400,3000)
#        noms[h].SetLineColor(ROOT.kBlack)
#        hems[h] = newNameAndStructure(hems[h],bkgnames[h]+"_HEM",titles['h_zp_jigm'][3],1400,3000)
#        hems[h].SetLineColor(ROOT.kRed)
#        makePlottableZprime(noms[h])
#        makePlottableZprime(hems[h])
#
#        hmax = hems[h].GetMaximum()
#        hems[h].SetMaximum(hmax*1.4)
#
#        div = hems[h].Clone()
#        div.Divide(hems[h],noms[h])
#        div.SetMarkerStyle(8)
#        div.SetLineColor(ROOT.kBlack)
#        div.SetMarkerColor(ROOT.kBlack)
#        divs.append(div)
#
#        
#    tc = ROOT.TCanvas('tc','HEM Check',1200,1600)
#    pd1 = ROOT.TPad('pd1','sb',0.0,0.15,0.5,0.5)
#    pd2 = ROOT.TPad('pd2','sr',0.5,0.15,1.0,0.5)
#    pd3 = ROOT.TPad('pd1','sb',0.0,0.65,0.5,1)
#    pd4 = ROOT.TPad('pd2','sr',0.5,0.65,1.0,1.0)
#    prat1 = ROOT.TPad('prat1','sbrat',0.0,0.0,0.5,0.15)
#    prat2 = ROOT.TPad('prat2','srrat',0.5,0.0,1.0,0.15)
#    prat3 = ROOT.TPad('prat3','sbrat',0.0,0.5,0.5,0.65)
#    prat4 = ROOT.TPad('prat4','srrat',0.5,0.5,1.0,0.65)
#    pd1.SetLeftMargin(0.15)
#    pd1.SetRightMargin(0.05)
#    pd2.SetLeftMargin(0.15)
#    pd2.SetRightMargin(0.05)
#    pd3.SetLeftMargin(0.15)
#    pd3.SetRightMargin(0.05)
#    pd4.SetLeftMargin(0.15)
#    pd4.SetRightMargin(0.05)
#
#    prat1.SetRightMargin(.05)
#    prat1.SetLeftMargin(0.15)
#    prat1.SetBottomMargin(0.25)
#    prat1.SetTopMargin(0.05)
#    prat2.SetRightMargin(.05)
#    prat2.SetLeftMargin(0.15)
#    prat2.SetBottomMargin(0.25)
#    prat2.SetTopMargin(0.05)
#    prat3.SetRightMargin(.05)
#    prat3.SetLeftMargin(0.15)
#    prat3.SetBottomMargin(0.25)
#    prat3.SetTopMargin(0.05)
#    prat4.SetRightMargin(.05)
#    prat4.SetLeftMargin(0.15)
#    prat4.SetBottomMargin(0.25)
#    prat4.SetTopMargin(0.05)
#
#    #The Diboson Plots
#    tc2 = ROOT.TCanvas('tc2','tc2 HEM Check',1200,800)
#    pd11 = ROOT.TPad('pd1','sb',0.0,0.3,0.5,1.0)
#    pd22 = ROOT.TPad('pd2','sr',0.5,0.3,1.0,1.0)
#    prat11 = ROOT.TPad('prat1','sbrat',0.0,0.0,0.5,0.3)
#    prat22 = ROOT.TPad('prat2','srrat',0.5,0.0,1.0,0.3)
#    prat11.SetRightMargin(.05)
#    prat11.SetLeftMargin(0.15)
#    prat11.SetBottomMargin(0.25)
#    prat11.SetTopMargin(0.05)
#    prat22.SetRightMargin(.05)
#    prat22.SetLeftMargin(0.15)
#    prat22.SetBottomMargin(0.25)
#    prat22.SetTopMargin(0.05)
#    pd11.SetLeftMargin(0.15)
#    pd11.SetRightMargin(0.05)
#    pd22.SetLeftMargin(0.15)
#    pd22.SetRightMargin(0.05)
#
#    tc2.Draw()
#    tc2.cd()
#    pd11.Draw()
#    pd11.cd()
#    hems[2].Draw('hist')
#    noms[2].Draw('histsame')
#    lab1 = makeDeviatedTPave(noms[2],hems[2],"WZ SR")
#    lab1.Draw()
#    leg1 = makeLegend(noms[2],hems[2],"WZ SR")
#    leg1.Draw()
#    tc2.cd()
#    prat11.Draw()
#    prat11.cd()
#    divs[2].Draw('pe')
#    tc2.cd()
#    pd22.Draw()
#    pd22.cd()
#    hems[3].Draw('hist')
#    noms[3].Draw('histsame')
#    lab2 = makeDeviatedTPave(noms[3],hems[3],"ZZ SR")
#    lab2.Draw()
#    leg2 = makeLegend(noms[3],hems[3],"ZZ SR")
#    leg2.Draw()
#    tc2.cd()
#    prat22.Draw()
#    prat22.cd()
#    divs[3].Draw('pe')
#    tc2.cd()
#
#    tc2.SaveAs("HEM_Check_VV_Plots.png")
#    


    #The Alpha Ratio Ones
    #tc.Draw()
    #tc.cd()
    #pd1.Draw()
    #pd1.cd()
    ##hs.Draw('nostack')
    #hems[0].Draw('hist')
    #noms[0].Draw('histsame')
    #lab1 = makeDeviatedTPave(noms[0],hems[0],"DY SR")
    #lab1.Draw()
    #leg1 = makeLegend(noms[0],hems[0],"DY SR")
    #leg1.Draw()
    #tc.cd()
    #prat1.Draw()
    #prat1.cd()
    #divs[0].Draw('pe')
    #tc.cd()
    #pd2.Draw()
    #pd2.cd()
    #hems[1].Draw('hist')
    #noms[1].Draw('histsame')
    #lab2 = makeDeviatedTPave(noms[1],hems[1],"DY Alpha SB")
    #lab2.Draw()
    #leg2 = makeLegend(noms[1],hems[1],"DY Alpha SB")
    #leg2.Draw()
    #tc.cd()
    #prat2.Draw()
    #prat2.cd()
    #divs[1].Draw('pe')
    #tc.cd()
    #pd3.Draw()
    #pd3.cd()
    #hems[2].Draw('hist')
    #noms[2].Draw('histsame')
    #lab3 = makeDeviatedTPave(noms[2],hems[2],"TT Alpha SB")
    #lab3.Draw()
    #leg3 = makeLegend(noms[2],hems[2],"TT")
    #leg3.Draw()
    #tc.cd()
    #prat3.Draw()
    #prat3.cd()
    #divs[2].Draw('pe')
    #tc.cd()
    #pd4.Draw()
    #pd4.cd()
    #hems[3].Draw('hist')
    #noms[3].Draw('histsame')
    #lab4 = makeDeviatedTPave(noms[3],hems[3],"VV Alpha SB")
    #lab4.Draw()
    #leg4 = makeLegend(noms[3],hems[3],"VV")
    #leg4.Draw()
    #tc.cd()
    #prat4.Draw()
    #prat4.cd()
    #divs[3].Draw('pe')
    #tc.cd()
    #
    #tc.SaveAs("HEM_Check_AlphaRatio_Plots.png")

    
    
        
        
        #hdict['h_h_pt'] = [hhptnom,hhptHEM]
        #hdict['h_metxy'] = [hmetnom,hmetHEM]
        #hdict['h_zp_jigm'] = [hzpnom,hzpHEM]

        #tc = ROOT.TCanvas('tc',bkg,1800,800)

        #for key in hdict.keys():
        #    rebins = [h.Rebin(titles[key][3]) for h in hdict[key]]
        #    print("   ",key)
        #    print("       Nominal: ",rebins[0].Integral())
        #    print("      HEMshift: ",rebins[1].Integral())
        #    print("      devi/nom: ",rebins[1].Integral()/rebins[0].Integral())

