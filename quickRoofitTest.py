import tdrstyle
import CMS_lumi
import ROOT
import gecorg_test as go
import configparser

ROOT.gROOT.SetBatch(ROOT.kTRUE)

#conditions with ratio
stkpadydims = [0.3,1.]
ratpadydims = [0.0,0.3]
tcanvasdims = [600,800]

bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
bkgcols = go.colsFromPalette(bkgnames,ROOT.kLake)

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
    #fitlabel.AddText("\Chi^2 / ndof = {0}".format(round(chi2/ndof,4)))
    fitlabel.SetFillColor(0)
    return chi2,ndof,fitlabel

def makePlotandPng(hmc,hdata,name,**kwargs):
    print("In the plotter")
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
            #chi2,ndof,fitinfo = makeTPadOfFitGOF(fitdict["nom"])
            #fitinfo.Draw()
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
    estf  = ROOT.TFile.Open('Run2_161718_emu_extrapolation_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')

    #Compile RooFit stuff
    print("Compiling fits")
    ROOT.gSystem.CompileMacro("generic_roofits.C","kfc")
    ROOT.gSystem.Load("generic_roofits_C")
    
    #the full emu estation
    hemufull = estf.Get('TT_scalenominal')
    hemu6TeV = go.newNameAndStructure(hemufull.Clone(),"TTdata6TeV",1,0,6000)

    bkgs = go.backgrounds(config.get('nominal','pathnom'),zptcut,hptcut,metcut,btagwp,config.get('nominal','strnom'))

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
    hsrtt = hsrtt.Rebin(len(newbinedges)-1,"TTmcLimitRange",newbinedges)#The 6 bin that go to combine
    hsrtt6TeV = go.newNameAndStructure(hsrttori.Clone(),"TTmc6TeV",2,0,6000)

    #do a fit
    #fit = ROOT.SimpExpFit(hemu6TeV,1300,4000)
    #print(fit)
    #print(fit.GetParameter(0))

    print("Trying Roofiting in pyroot")

    #do fit in pyroot
    x = ROOT.RooRealVar("x","RJR",1300,4000)
    dh = ROOT.RooDataHist("dh","emu RJR dist",[x],Import=hemu6TeV)
    x.setRange(1300,4000)
    initau = ROOT.guessDecayConstantTH1F(hemu6TeV,hemu6TeV.GetMaximum())
    tau = ROOT.RooRealVar("tau","tau",initau,-5,5)
    model = ROOT.RooExponential("model","model",x,tau)
    model.fitTo(dh)
    print(tau)
    fit = model.asTF(x,[tau],x)
    print(fit.GetParameter(0))

    #do another fit
    modelext = ROOT.RooExtendPdf("modelext","modelext",model,hemu6TeV.Integral(1200,4000))

    #distnorm = hemu6TeV.Integral(1200,4000)
    #fit2 = distnorm*fit
    
    
    tc = ROOT.TCanvas('tc','tc',600,600)
    tc.Draw()
    tc.cd()
    fit.Draw()
    tc.SaveAs("test_roofit.png")
    #make plot
    makePlotandPng(hsrtt6TeV,hemu6TeV,"6TeVemu_roofitexp",fitdict = {"nom":fit,"name":"RooExp"})
