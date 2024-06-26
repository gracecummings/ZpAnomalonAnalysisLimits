import argparse
import tdrstyle
import CMS_lumi
import ROOT
import gecorg_test as go
import configparser
import glob

ROOT.gROOT.SetBatch(ROOT.kTRUE)

#conditions with ratio
stkpadydims = [0.3,1.]
ratpadydims = [0.0,0.3]
tcanvasdims = [600,800]

bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
bkgcols = go.colsFromPalette(bkgnames,ROOT.kLake)

def castFitIntoHistogram(empty,dicfits,fitlow):
    histlist = []
    for key in dicfits.keys():
        name = key
        #print("****************************",key)
        if "up" in key:
            name = key.replace("up","Up")
        if "dwn" in key:
            name = key.replace("dwn","Down")
        hnew = empty.Clone()
        hnew.SetName(name)
        hnew.SetTitle(name)
        fit = dicfits[key]
        for b in range(hnew.GetNbinsX()+1):
            if hnew.GetBinCenter(b) < fitlow:
                hnew.SetBinContent(b,-999.0)
            else:
                if ((hnew.GetBinCenter(b-1) >= fitlow) and (fit.Eval(hnew.GetBinCenter(b-1)) < fit.Eval(hnew.GetBinCenter(b)))):
                    #print("fit turned over - bin ",b)
                    #print("   b-1 fit val  ",fit.Eval(hnew.GetBinCenter(b-1)))
                    #print("   b-1 bin cont ",hnew.GetBinContent(b-1))
                    #print("   b   fit val  ",fit.Eval(hnew.GetBinCenter(b)))
                    
                    hnew.SetBinContent(b,hnew.GetBinContent(b-1))
                    hnew.SetBinError(b,0.0)
                    #print("   b bin after set ",hnew.GetBinContent(b))
                else:
                    #print("Filiing bin with center: ",hnew.GetBinCenter(b))
                    #print("  Filling with: ",fit.Eval(hnew.GetBinCenter(b)))
                    hnew.SetBinContent(b,fit.Eval(hnew.GetBinCenter(b)))
                    hnew.SetBinError(b,0.0)
        histlist.append(hnew)
    return histlist


def makeTPadOfFitGOF(fit,normfits = False):
    lims = [0.55,0.3,.9,0.45]
    if normfits:
        lims = [0.17,0.8,0.52,0.9]
    chi2 = fit.GetChisquare()
    ndof = fit.GetNDF()
    fitlabel = ROOT.TPaveText(lims[0],lims[1],lims[2],lims[3],"NBNDC")
    #print("chi 2 ",chi2)
    #print("\Chi^2 = {0} , ndof = {1}".format(round(chi2,2),ndof))
    fitlabel.AddText("Simple exp fit to estimation")
    fitlabel.AddText("likelihood fit")
    fitlabel.AddText("\Chi^2 = {0} , ndof = {1}".format(round(chi2,2),ndof))
    fitlabel.AddText("\Chi^2 / ndof = {0}".format(round(chi2/ndof,4)))
    fitlabel.SetFillColor(0)
    return chi2,ndof,fitlabel

def makePlotandPng(hmc,hdata,name,**kwargs):
    fitdict = kwargs.get('fitdict',None)
    empty = hdata.Clone()
    empty.SetDirectory(0)
    empty.Reset("ICESM")

    #Plotting sytle
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = go.lumiFormatter(years)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Internal"

    #Do some settings, make some plots
    hmc.SetFillColor(bkgcols[1])
    hmc.SetLineColor(bkgcols[1])
    hdata.SetBinErrorOption(1)
    hdata.SetMarkerStyle(8)
    hdata.SetMarkerSize(0.7)
    hdata.SetMarkerColor(ROOT.kBlack)

    #make ratio
    hdiv = hdata.Clone()
    hdiv.Divide(hmc)
    ratline = ROOT.TLine(hdiv.GetBinLowEdge(1),1,hdiv.GetBinLowEdge(hdiv.GetNbinsX())+hdiv.GetBinWidth(hdiv.GetNbinsX()),1)
    
    leg = ROOT.TLegend(0.5,0.45,0.92,0.80)
    leg.SetBorderSize(0)
    leg.AddEntry(hmc,"ttbar MC in SR","f")
    leg.AddEntry(hdata,"data-driven estimation")

    #Plotting itself
    tc = ROOT.TCanvas("tc","tc",tcanvasdims[0],tcanvasdims[1])
    p1 = ROOT.TPad("p1","ttbar_comparison",0,stkpadydims[0],1.0,stkpadydims[1])
    p1.SetLeftMargin(0.15)
    p1.SetRightMargin(0.05)
    p2 = ROOT.TPad("p2","ttbar_estimation_ratio",0,ratpadydims[0],1.0,ratpadydims[1])
    p2.SetRightMargin(.05)
    p2.SetLeftMargin(0.15)
    p2.SetBottomMargin(0.25)
    p2.SetTopMargin(0.05)

    #Prepare first pad for stack
    p1.Draw()
    p1.cd()
    #p1.SetLogy()
        
    #Draw the stack
    hmc.Draw("HIST")#add PFC for palette drawing
    hmc.GetXaxis().SetTitle("RJR Estimator")
    hmc.GetXaxis().SetTitleSize(0.05)
    hmc.GetXaxis().SetTitleOffset(1.1)
    hmc.GetXaxis().SetLabelSize(0.04)
    hmc.GetYaxis().SetTitle("Events")
    hmc.GetYaxis().SetTitleSize(0.05)
    hmc.GetYaxis().SetLabelSize(0.04)
    CMS_lumi.CMS_lumi(p1,4,13)
    hdata.Draw("histsame,pe")
    p1.Update()

    if fitdict:#should be the fits
        keys = fitdict.keys()
        if len(keys) > 2:
            #make new forms of stuff
            hshifts = castFitIntoHistogram(empty,{"h"+fitdict["name"]+"up":fitdict["up"],"h"+fitdict["name"]+"dwn":fitdict["dwn"],"hnominal":fitdict["nom"]},1300)
            hshiftsdvs = []
            for h in hshifts:
                hshiftsdvs.append(h.Clone())
                hshiftsdvs[-1].Divide(hshifts[-1])

            #plot
            fitdict["nom"].SetLineColor(ROOT.kBlack)
            fitdict["up"].SetLineColor(2)
            fitdict["dwn"].SetLineColor(4)
            fitdict["nom"].Draw("same")
            fitdict["up"].Draw("same")
            fitdict["dwn"].Draw("same")

            leg.AddEntry(fitdict["nom"],"Nominal Exp Fit")
            leg.AddEntry(fitdict["up"],fitdict["name"]+" up")
            leg.AddEntry(fitdict["dwn"],fitdict["name"]+" dwn")
            leg.Draw()

            p1.Update()
            tc.cd()
            p2.Draw()
            p2.cd()

            hshiftsdvs[-1].SetLineColor(ROOT.kBlack)
            hshiftsdvs[-1].Draw()
            hshiftsdvs[-1].GetXaxis().SetTitle("bin center")
            hshiftsdvs[-1].GetXaxis().SetTitleSize(0.11)
            hshiftsdvs[-1].GetXaxis().SetTitleOffset(0.65)
            hshiftsdvs[-1].GetXaxis().SetLabelSize(0.075)
            hshiftsdvs[-1].GetYaxis().SetTitle("syst/fit")
            hshiftsdvs[-1].GetYaxis().SetTitleSize(0.11)
            hshiftsdvs[-1].GetYaxis().SetTitleOffset(.45)
            hshiftsdvs[-1].GetYaxis().SetLabelSize(0.08)
            hshiftsdvs[-1].GetYaxis().SetLabelOffset(0.02)
            hshiftsdvs[-1].GetYaxis().SetNdivisions(503)
            hshiftsdvs[-1].SetMinimum(-0.5)
            hshiftsdvs[-1].SetMaximum(5.)
            hshiftsdvs[0].SetLineColor(2)
            hshiftsdvs[1].SetLineColor(4)
            hshiftsdvs[0].Draw("same")
            hshiftsdvs[1].Draw("same")
            hshiftsdvs[-1].Draw("same")

            legfit = ROOT.TLegend(0.25,0.6,0.55,0.92)
            legfit.SetBorderSize(0)
            legfit.AddEntry(hshiftsdvs[-1],"Nominal Fit Value","l")
            legfit.AddEntry(hshiftsdvs[0],fitdict["name"]+" up","l")
            legfit.AddEntry(hshiftsdvs[1],fitdict["name"]+" dwn","l")
            legfit.Draw()
            tc.cd()

        else:
            #divisions
            fitdict["nom"].SetLineColor(ROOT.kBlack)
            fitdict["nom"].Draw("same")
            leg.AddEntry(fitdict["nom"],"Nominal Exp Fit")
            chi2,ndof,fitinfo = makeTPadOfFitGOF(fitdict["nom"])
            fitinfo.Draw()
            leg.Draw()
            p1.Update()
            tc.cd()
            p2.Draw()
            p2.cd()
            hdiv.Draw("pe")
            hdiv.GetXaxis().SetTitle("bin center")
            hdiv.GetXaxis().SetTitleSize(0.11)
            hdiv.GetXaxis().SetTitleOffset(0.65)
            hdiv.GetXaxis().SetLabelSize(0.075)
            hdiv.GetYaxis().SetTitle("est/MC")
            hdiv.GetYaxis().SetTitleSize(0.11)
            hdiv.GetYaxis().SetTitleOffset(.45)
            hdiv.GetYaxis().SetLabelSize(0.08)
            hdiv.GetYaxis().SetLabelOffset(0.02)
            hdiv.GetYaxis().SetNdivisions(503)
            hdiv.SetMinimum(-0.5)
            hdiv.SetMaximum(2.5)
            ratline.Draw()
            tc.cd()

    else:
        leg.Draw()
        p1.Update()
        tc.cd()
        p2.Draw()
        p2.cd()
        hdiv.Draw("pe")
        hdiv.GetXaxis().SetTitle("bin center")
        hdiv.GetXaxis().SetTitleSize(0.11)
        hdiv.GetXaxis().SetTitleOffset(0.65)
        hdiv.GetXaxis().SetLabelSize(0.075)
        hdiv.GetYaxis().SetTitle("est/MC")
        hdiv.GetYaxis().SetTitleSize(0.11)
        hdiv.GetYaxis().SetTitleOffset(.45)
        hdiv.GetYaxis().SetLabelSize(0.08)
        hdiv.GetYaxis().SetLabelOffset(0.02)
        hdiv.GetYaxis().SetNdivisions(503)
        hdiv.SetMinimum(0.)
        hdiv.SetMaximum(2.)
        ratline.Draw()
        tc.cd()

    pngname = go.makeOutFile('ttbar_'+name,'ratio_161718_signalr','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(pngname)

def makeCombineHistandPng(hdata,hfitnom,hfitup,hfitdn,name):
    #can only take un-rebinned plots
    empty = hdata.Clone()
    empty.SetDirectory(0)
    empty.Reset("ICESM")

    #Rebin the emu data 
    hemu = hdata.Clone()
    hemu = go.newNameAndStructure(hemu,"TT_emu_data",1,limrangelow,limrangehigh)
    hemu = hemu.Rebin(len(newbinedges)-1,"TT_emu_data",newbinedges)#The 6 bin that go to combine
    int_emu = hemu.Integral()

    #gather the fits, make them hists, and rebin
    fitnames = ["TT_"+name+"_up","TT_"+name+"_dwn","TT_nom_fit"]
    hfits = []
    hfits.append(castFitIntoHistogram(empty,{fitnames[0]:hfitup},1300)[0])
    hfits.append(castFitIntoHistogram(empty,{fitnames[1]:hfitdn},1300)[0])  
    hfits.append(castFitIntoHistogram(empty,{fitnames[2]:hfitnom},1300)[0])
    hfits_use = []
    
    for i,hfit in enumerate(hfits):
        hfits_use.append(go.newNameAndStructure(hfit,fitnames[i],1,limrangelow,limrangehigh))
        hfits_use[-1] = hfits_use[-1].Rebin(len(newbinedges)-1,fitnames[i],newbinedges)

    print("Checking Yields")
    print("    ",name)
    int_nom = hfits_use[-1].Integral()
    int_up  = hfits_use[0].Integral()
    int_dwn  = hfits_use[1].Integral()
    print("    nom: ",int_nom)
    print("     up: ",int_up)
    print("    dwn: ",int_dwn)


    rat_up = int_up/int_nom
    rat_dwn = int_dwn/int_nom

    #make divs of the fits -- dividing each by nominal fit hist
    hshiftsdvs = []
    for h in hfits_use:
        #h.SetBinErrorOption(1)
        hshiftsdvs.append(h.Clone())
        hshiftsdvs[-1].Divide(hfits_use[-1])

    #Plotting sytle
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = go.lumiFormatter(years)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Internal"

    #Do some settings, make some plots
    hdata.SetBinErrorOption(1)
    hdata.SetMarkerStyle(8)
    hdata.SetMarkerSize(0.7)
    hdata.SetMarkerColor(ROOT.kBlack)
    hfits_use[-1].SetLineColor(ROOT.kBlack)
    hfits_use[0].SetLineColor(2)
    hfits_use[1].SetLineColor(4)

    leg = ROOT.TLegend(0.5,0.45,0.92,0.80)
    leg.SetBorderSize(0)
    leg.AddEntry(hdata,"data-driven estimation")
    leg.AddEntry(hfits_use[-1],"Nominal Exp Fit")
    leg.AddEntry(hfits_use[0],fitnames[0])
    leg.AddEntry(hfits_use[1],fitnames[1])

    fitlabel = ROOT.TPaveText(0.55,0.3,0.9,0.45,"NBNDC")
    fitlabel.AddText("emu yield: "+str(round(int_emu,4)))
    fitlabel.AddText("fit yield: "+str(round(int_nom,4)))
    fitlabel.AddText(name+" up effect: "+str(round(rat_up,4)))
    fitlabel.AddText(name+" dwn effect: "+str(round(rat_dwn,4)))
    fitlabel.SetFillColor(0)

    #Plotting itself
    tc = ROOT.TCanvas("tc","tc",tcanvasdims[0],tcanvasdims[1])
    p1 = ROOT.TPad("p1","ttbar_comparison",0,stkpadydims[0],1.0,stkpadydims[1])
    p1.SetLeftMargin(0.15)
    p1.SetRightMargin(0.05)
    p2 = ROOT.TPad("p2","ttbar_estimation_ratio",0,ratpadydims[0],1.0,ratpadydims[1])
    p2.SetRightMargin(.05)
    p2.SetLeftMargin(0.15)
    p2.SetBottomMargin(0.25)
    p2.SetTopMargin(0.05)

    #Prepare first pad for stack
    p1.Draw()
    p1.cd()
    hfits_use[0].Draw("hist")
    hfits_use[0].GetXaxis().SetTitle("RJR Estimator")
    hfits_use[0].GetXaxis().SetTitleSize(0.05)
    hfits_use[0].GetXaxis().SetTitleOffset(1.1)
    hfits_use[0].GetXaxis().SetLabelSize(0.04)
    hfits_use[0].GetYaxis().SetTitle("Events")
    hfits_use[0].GetYaxis().SetTitleSize(0.05)
    hfits_use[0].GetYaxis().SetLabelSize(0.04)
    CMS_lumi.CMS_lumi(p1,4,13)
    hfits_use[1].Draw("hist,same")
    hfits_use[-1].Draw("hist,same")
    hdata.Draw("histsame,pe")
    leg.Draw()
    fitlabel.Draw()
    p1.Update()
    tc.cd()


    p2.Draw()
    p2.cd()
    hshiftsdvs[-1].SetLineColor(ROOT.kBlack)
    hshiftsdvs[-1].Draw()
    hshiftsdvs[-1].GetXaxis().SetTitle("bin center")
    hshiftsdvs[-1].GetXaxis().SetTitleSize(0.11)
    hshiftsdvs[-1].GetXaxis().SetTitleOffset(0.65)
    hshiftsdvs[-1].GetXaxis().SetLabelSize(0.075)
    hshiftsdvs[-1].GetYaxis().SetTitle("syst/fit")
    hshiftsdvs[-1].GetYaxis().SetTitleSize(0.11)
    hshiftsdvs[-1].GetYaxis().SetTitleOffset(.45)
    hshiftsdvs[-1].GetYaxis().SetLabelSize(0.08)
    hshiftsdvs[-1].GetYaxis().SetLabelOffset(0.02)
    hshiftsdvs[-1].GetYaxis().SetNdivisions(503)
    hshiftsdvs[-1].SetMinimum(-0.5)
    hshiftsdvs[-1].SetMaximum(5.)
    hshiftsdvs[0].SetLineColor(2)
    hshiftsdvs[1].SetLineColor(4)
    hshiftsdvs[0].Draw("same")
    hshiftsdvs[1].Draw("same")
    hshiftsdvs[-1].Draw("same")
    
    legfit = ROOT.TLegend(0.25,0.6,0.55,0.92)
    legfit.SetBorderSize(0)
    legfit.AddEntry(hshiftsdvs[-1],"Nominal Fit Value","l")
    legfit.AddEntry(hshiftsdvs[0],fitnames[0],"l")
    legfit.AddEntry(hshiftsdvs[1],fitnames[1],"l")
    legfit.Draw()
    tc.cd()

    pngname = go.makeOutFile('ttbar_'+name,'ratio_161718_combineplots','.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    tc.SaveAs(pngname)
    
    print("recapping the histogram integral")
    for h in hfits_use:
        print("    ",h.Integral())

    return hfits_use


def doExpShifts(nomfit,parnum,name,lowr,highr,shiftedparamsin = ROOT.TVector()):
    #print("doing the shifted fits!")
    fitpars = []
    fiterrs = []
    shiftedfits = []
    for i in range(parnum):
        parvecedit = ROOT.TVector(parnum)
        for j in range(parnum):
            parvecedit[j] = shiftedparamsin[i][j]
        fit = ROOT.expFitSetParsAndErrs(name+"par"+str(i),parvecedit,lowr,highr)
        shiftedfits.append(fit)
    return fitpars,fiterrs,shiftedfits

if __name__=='__main__':
    #Years
    years   = [16,17,18]
    yearstr = go.yearFormatter(years)

    #Get systematic info
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)

    #Starting parameters
    zptcut  = '100.0'
    hptcut  = '300.0'
    metcut  = '75.0'
    btagwp  = '0.8'
    rebindiv = 2
    limrangelow = 1800
    limrangehigh = 10000

    #import files
    print("Importing root files")
    combf = glob.glob('unblindedDatacardHolder/Run2_161718_ZllHbbMET_unblind_mumu*.root')
    estf  = ROOT.TFile.Open('Run2_161718_emu_extrapolation_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    bkgs = go.backgrounds(config.get('nominal','pathnom'),zptcut,hptcut,metcut,btagwp,config.get('nominal','strnom'))
    #data = go.run2(config.get('nominal','pathdata'),zptcut,hptcut,metcut,btagwp,config.get('nominal','strdata'))

    #Get some plots
    #mc
    tf1 = ROOT.TFile.Open(bkgs.f17dyjetsr[0])
    htest = tf1.Get('h_zp_jigm')
    empty = htest.Clone()
    empty.SetDirectory(0)
    tf1.Close()
    empty.Reset("ICESM")
    hsrttori = bkgs.getAddedHist(empty,"TT","sr","h_zp_jigm")
    hsrttori.SetDirectory(0)
    hsrtt = go.newNameAndStructure(hsrttori.Clone(),"TTmc",rebindiv,limrangelow,limrangehigh)
    newbinedges = go.makeBinLowEdges(hsrtt,2800)#last normal bin
    print("Getting bin edges for what is passed to Combine")
    print("    the edges: ",newbinedges)
    hsrtt = hsrtt.Rebin(len(newbinedges)-1,"TTmcLimitRange",newbinedges)#The 6 bin that go to combine

    hsrtt6TeV = go.newNameAndStructure(hsrttori.Clone(),"TTmc6TeV",2,0,6000)

    #What is in the combine file
    tf2 = ROOT.TFile.Open(combf[0])#does not matter which, all have same bkg
    hesttt = tf2.Get('TT')
    hesttt.SetDirectory(0)
    tf2.Close()

    #the full emu estation
    hemufull = estf.Get('TT_scalenominal')
    hemu6TeV = go.newNameAndStructure(hemufull.Clone(),"TTdata6TeV",1,0,6000)

    #Lets do some nominal fits
    print("Compiling fits")
    ROOT.gSystem.CompileMacro("../ZpAnomalonAnalysisUproot/cfunctions/alphafits.C","kfc")
    ROOT.gSystem.Load("../ZpAnomalonAnalysisUproot/cfunctions/alphafits_C")
    print("=================       Doing base fit        =================")
    emuexpfit = ROOT.expFitTH1F(hemu6TeV,"emuexp","LRE0+",1300,4000)
    hemuexp = castFitIntoHistogram(empty.Clone(),{"hemuexp":emuexpfit},1300)[0]

        #Lets do an alternate fit
    print("=================       Doing Alt fit         =================")
    #print("****offset***")
    #emuexpoffsetfit = ROOT.expOffsetFitTH1F(hemu6TeV,"emuexpoffset","LRE0+",1300,4000)
    print("****sqr***")
    emuexpsqrfit = ROOT.expSqFitTH1F(hemu6TeV,"emuexpsqr","LRE0+",1300,4000)
    #print("****pwr***")
    #emuexppwrfit = ROOT.expPwrFitTH1F(hemu6TeV,"emuepwr","LR0+",1300,4000)
    #print("****poly3***")
    #hemu6TeV.Fit("pol3","L+")
    #emupoly3fit = hemu6TeV.GetFunction("pol3")
    #emupoly3fit = ROOT.poly3FitTH1F(hemu6TeV,"poly3","LRE0+",1300,3000)

    #print("****poly2***")
    #emupoly2fit = ROOT.poly2FitTH1F(hemu6TeV,"poly2","LRE0+",1100,3000)
    #print("*******line******")
    #emulinefit = ROOT.lineFitTH1F(hemu6TeV,"line","LRE0+",1300,4000)

    print("****xExp***")
    emuxexpfit = ROOT.xexpFitTH1F(hemu6TeV,"xexp","LRE0+",1100,4000)
    #print("****x2Exp***")
    #emux2expfit = ROOT.x2expFitTH1F(hemu6TeV,"x2exp","LRE0+",1100,4000)

    #altfit_syst = {"nom":emuexpfit,"up":emuexpoffsetfit,"dwn":emuexpsqrfit,"name":"AltFunc"}
    #altfit_syst1 = {"nom":emuexpfit,"up":emupoly2fit,"dwn":emuexpsqrfit,"name":"AltFunc"}
    #altfit_syst2 = {"nom":emuexpfit,"up":emuexppwrfit,"dwn":emuexpsqrfit,"name":"AltFunc"}
    altfit_syst3 = {"nom":emuexpfit,"up":emuxexpfit,"dwn":emuexpsqrfit,"name":"AltFunc"}
    #altfit_syst4 = {"nom":emuexpfit,"up":emux2expfit,"dwn":emuexpsqrfit,"name":"AltFunc"}

    #Doing the systematic shifted fits
    print("================= Doing fits to gather errors =================")
    emuexpfitparamsup = ROOT.expFitTH1DecorrParamsShiftedUp(hemu6TeV,"emushiftup","LQRE0+",1300,4000)
    emuexpfitparamsdn = ROOT.expFitTH1DecorrParamsShiftedDown(hemu6TeV,"emushiftdn","LQRE0+",1300,4000)
    print("=================     Getting Shifted Fits    =================")
    emuupparams,emuuperrs,emuupfits = doExpShifts(emuexpfit,2,"ttup",1300,4000,emuexpfitparamsup)
    emudnparams,emudnerrs,emudnfits = doExpShifts(emuexpfit,2,"ttdwn",1300,4000,emuexpfitparamsdn)
    par0_syst = {"nom":emuexpfit,"up":emuupfits[0],"dwn":emudnfits[0],"name":"Par0"}
    par1_syst = {"nom":emuexpfit,"up":emuupfits[1],"dwn":emudnfits[1],"name":"Par1"}
    
    #Make Plots
    makePlotandPng(hsrtt,hesttt,"estimation_over_mc")
    makePlotandPng(hsrttori.Rebin(2),hemufull,"fullemuregion_mccomp")
    makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_mccomp",fitdict = {"nom":emuexpfit,"name":"nominal exp"})
    makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncexpsqr",fitdict = {"nom":emuexpsqrfit,"name":"alt exp sqr"})
    #makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncexppwr",fitdict = {"nom":emuexppwrfit,"name":"alt exp pwr"})
    #makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncexpoff",fitdict = {"nom":emuexpoffsetfit,"name":"alt exp offset"})
    ##makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncpoly3",fitdict = {"nom":emupoly3fit,"name":"alt ploy3"})
    #makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncpoly2",fitdict = {"nom":emupoly2fit,"name":"alt ploy2"})
    #makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncline",fitdict = {"nom":emulinefit,"name":"alt line"})
    makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncxexp",fitdict = {"nom":emuxexpfit,"name":"alt xexp"})
    makePlotandPng(hsrtt6TeV,hemu6TeV,"6Tevemu_par0syst",fitdict = par0_syst)
    makePlotandPng(hsrtt6TeV,hemu6TeV,"6Tevemu_par1syst",fitdict = par1_syst)
    makePlotandPng(hsrttori,hemufull,"fullemuregion_par0syst",fitdict = par0_syst)
    makePlotandPng(hsrttori,hemufull,"fullemuregion_par1syst",fitdict = par1_syst)
    #makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncts_exp",fitdict =altfit_syst)
    #makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncts_poly2exp",fitdict =altfit_syst1)
    #makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncts_pwrexp",fitdict =altfit_syst2)
    makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncts_xexp",fitdict =altfit_syst3)
    #makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_altfuncts_x2exp",fitdict =altfit_syst4)

    #Do Combine Plots 
    #do check of recovered ttbar file and combine submission dists
    hemufromfullspec = hemufull.Clone()
    hemufromfullspec.SetDirectory(0)
    hemufromfullspec = go.newNameAndStructure(hemufromfullspec,"TTemuFromFullinLimitRange",1,limrangelow,limrangehigh)
    hemufromfullspec = hemufromfullspec.Rebin(len(newbinedges)-1,"TTemuFromFullinLimitRange",newbinedges)#The 6 bin that go to combine
    makePlotandPng(hemufromfullspec,hesttt,"combineest_over_estfromrecovdredfile")
    par0hists = makeCombineHistandPng(hemufull.Clone(),emuexpfit,emuupfits[0],emudnfits[0],"Par0")
    print(par0hists[0].Integral())
    print(par0hists[1].Integral())
    print(par0hists[2].Integral())
    par1hists = makeCombineHistandPng(hemufull.Clone(),emuexpfit,emuupfits[1],emudnfits[1],"Par1")
    #altfunchists = makeCombineHistandPng(hemufull.Clone(),emuexpfit,emuexpoffsetfit,emuexpsqrfit,"AltFunctOffset")
    #altfunchists1 = makeCombineHistandPng(hemufull.Clone(),emuexpfit,emupoly2fit,emuexpsqrfit,"AltFunct")
    #altfunchists2 = makeCombineHistandPng(hemufull.Clone(),emuexpfit,emuexppwrfit,emuexpsqrfit,"AltFunctPwr")
    altfunchists3 = makeCombineHistandPng(hemufull.Clone(),emuexpfit,emuxexpfit,emuexpsqrfit,"AltFunctxExp")
    #altfunchists4 = makeCombineHistandPng(hemufull.Clone(),emuexpfit,emux2expfit,emuexpsqrfit,"AltFunctx2Exp")

    outname = go.makeOutFile('ttbar_shenanigan_fits_to_data','ratio_161718_signalr','.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    rootout = ROOT.TFile(outname,'recreate')
    hemuexp.Write()

    for h in par0hists:
        h.Write()
    for h in par1hists[:-1]:
        h.Write()
    #for h in altfunchists:
    #    h.Write()
    #for h in altfunchists1:
    #    h.Write()
    #for h in altfunchists2:
    #    h.Write()
    for h in altfunchists3[:-1]:
        h.Write()
    #for h in altfunchists4:
    #    h.Write()

    rootout.Close()

