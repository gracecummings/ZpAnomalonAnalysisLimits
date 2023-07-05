import ROOT
import gecorg_test as go
import tdrstyle
import CMS_lumi
import numpy as np

def makeStackedPlot(empty,bkgs,reg,years,bkgcols,dyEst,dynorm,leg,limrangelow,limrangehigh):
    #empty1 = empty.Clone()
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    empty4 = empty.Clone()
    empty5 = empty.Clone()
    empty6 = empty.Clone()
    
    #Gather histograms
    #hdy  = bkgs.getAddedHist(empty2,"DYJetsToLL",reg,hname,years = years)
    if "sr" in reg:
        hdy  = dyEst.Get("extrphistnoerrs").Clone()
    else:
        hdy  = bkgs.getAddedHist(empty2,"DYJetsToLL",reg,hname,years = years)
        hdy.Scale(dynorm)
    htt  = bkgs.getAddedHist(empty3,"TT",reg,hname,years = years)
    hzz  = bkgs.getAddedHist(empty4,"ZZTo2L2Q",reg,hname,years = years)
    hwz  = bkgs.getAddedHist(empty5,"WZTo2L2Q",reg,hname,years = years)
    
    hvv = hzz.Clone()
    hvv.Add(hwz)

    #do rebining
    htt = go.newNameAndStructure(htt,"TT",2,limrangelow,limrangehigh)
    hdy = go.newNameAndStructure(hdy,"DY",1,limrangelow,limrangehigh)
    hvv = go.newNameAndStructure(hvv,"VV",2,limrangelow,limrangehigh)

    if "sr" not in reg:
        hdy = go.newNameAndStructure(hdy,"DY",2,limrangelow,limrangehigh)
    
    newbinedges = go.makeBinLowEdges(hvv,2800)
    hvv = hvv.Rebin(len(newbinedges)-1,"VV",newbinedges)
    htt = htt.Rebin(len(newbinedges)-1,"TT",newbinedges)
    hdy = hdy.Rebin(len(newbinedges)-1,"DY",newbinedges)

    hdy.SetFillColor(bkgcols[0])
    hdy.SetLineColor(bkgcols[0])
    htt.SetFillColor(bkgcols[1])
    htt.SetLineColor(bkgcols[1])
    hvv.SetFillColor(bkgcols[2])
    
    #add lines to legend
    leg.AddEntry(htt,"TT","f")
    leg.AddEntry(hdy,"DYJetsToLL","f")
    #leg.AddEntry(hwz,"WZTo2L2Q","f")
    #leg.AddEntry(hzz,"ZZTo2L2Q","f")
    leg.AddEntry(hvv,"VV","f")

    #Bkg Stack for Plotting
    hsbkg = ROOT.THStack("hsbkg","")
    #hsbkg.Add(hzz)
    #hsbkg.Add(hwz)
    hsbkg.Add(hvv)
    hsbkg.Add(hdy)
    hsbkg.Add(htt)
    hsbkg.SetMinimum(0)
    hsbkg.SetMaximum(50)

    #Make added hist for ratio plotting
    hbkg = hvv.Clone()
    #hbkg.Add(hwz)
    #hbkg.Add(hzz)
    #hbkg.Add(hvv)
    hbkg.Add(htt)
    hbkg.Add(hdy)

    return hsbkg,hbkg


def regionFormatter(regionstr):
    regdict = {'sideband':'sb','signalr':'sr','totalr':'tr'}
    formatted = regdict[regionstr]
    return formatted

def lumiFormatter(yearlist):
    lumidict = {16:36.31,17:41.53,18:59.74}
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


if __name__=='__main__':

    #Get parameters
    sig_xsec  = "1.0" #args.xsec
    zptcut    = "100.0"#args.zptcut
    hptcut    = "300.0"#args.hptcut
    metcut    = "75.0"#args.metcut
    btagwp    = "0.8"#args.btagwp
    year      = False
    regname   = "sideband"#args.region
    #lumi      = 0
    pathplots = "mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref"#args.directory
    systr     = "systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco"#args.syst
    plot_data = True#args.data
    sigdivsor = 5
    hname     = "h_zp_jigm"

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


    #Gather files
    bkgs  = go.backgrounds(pathplots,zptcut,hptcut,metcut,btagwp,systr)
    #data  = go.run2(pathplots,zptcut,hptcut,metcut,btagwp,systr.replace("_mutrignom",""))
    data  = go.run2(pathplots,zptcut,hptcut,metcut,btagwp,'systnominal_btagnom_muidnom')
    sigs =  go.signal(pathplots,zptcut,hptcut,metcut,btagwp,sig_xsec,years,systr)
    dynorm = np.load(pathplots+'/Run2_161718_dynormalization_alphat_'+systr+'_signalblind_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.npy')[0]
    dyEst = ROOT.TFile('analysis_output_ZpAnomalon/2023-05-09/alpha_method_ttbar_likli/Run2_161718_dy_extraploationalphat_systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco__Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')

    #Colors, Naming, general style
    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = lumiFormatter(years)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Simulation Preliminary"

    #Gather plots
    testyear = years[0]#picks first year in list, so desired year if only one
    testfile = bkgs.bkgs["DYJetsToLL"][testyear][reg][0][0]#stacked plots should always have DY
    testtfile = ROOT.TFile(testfile)
    
    #Make holder histograms
    h = testtfile.Get(hname)
    empty = h.Clone()
    empty.Reset("ICESM")#creates an empty hist with same structure

    #Make legends
    legsb = ROOT.TLegend(0.55,0.5,0.9,0.8)
    legsb.SetBorderSize(0)
    legsr = ROOT.TLegend(0.55,0.5,0.9,0.8)
    legsr.SetBorderSize(0)
    legsbs = ROOT.TLegend(0.55,0.5,0.9,0.8)
    legsbs.SetBorderSize(0)

    #Make the sideband MC plots
    hsbkgsb, hbkgsb = makeStackedPlot(empty,bkgs,'sb',years,bkgcols,dyEst,dynorm,legsb,1800,10000)
    print("doing signal region")
    hsbkgsr, hbkgsr = makeStackedPlot(empty,bkgs,'sr',years,bkgcols,dyEst,dynorm,legsr,1800,10000)

    #Get the Data plot
    hdat = data.getAddedHist(empty.Clone(),'sb',hname,years=years)
    hdat = go.newNameAndStructure(hdat,"data_obs",2,1800,10000)
    newbinedges = go.makeBinLowEdges(hdat,2800)
    hdat = hdat.Rebin(len(newbinedges)-1,"data_obs",newbinedges)
 
    #hdat.SetBinErrorOption(1)
    hdat.SetMarkerStyle(8)
    hdat.SetMarkerSize(0.7)
    hdat.SetMarkerColor(ROOT.kBlack)
    legsb.AddEntry(hdat,"Data SB","pe")
    legsr.AddEntry(hdat,"Data SB","pe")

    #Scale the CR
    hdats = hdat.Clone()
    srscale = hbkgsr.Integral()/hdat.Integral()
    hdats.Scale(srscale)
    legsbs.AddEntry(hdats,"SB Data scaled to SR MC Yield","pe")
    

    #Get ratio plots
    hdivsb = hdat.Clone()
    hdivsr = hdat.Clone()
    hdivsbs = hdats.Clone()
    hdivsb.Divide(hdat,hbkgsb)
    hdivsr.Divide(hdat,hbkgsr)
    hdivsbs.Divide(hdats,hbkgsr)
    ratline = ROOT.TLine(hbkgsb.GetBinLowEdge(1),1,hbkgsb.GetBinLowEdge(hbkgsb.GetNbinsX())+hbkgsb.GetBinWidth(hbkgsb.GetNbinsX()),1)

    #Plots
    tc = ROOT.TCanvas("tc",hname+"GoF",1800,800)
    pdsbonly = ROOT.TPad("pdsbonly","SB MC SB Data",0,.25,.33,1)
    pdsbassr = ROOT.TPad("pdsbassr","SR MC SB Data",.33,.25,.67,1)
    pdsbscaled = ROOT.TPad("pdpdsbscaled","SR MC Scaled SB Data",.67,.25,1,1)
    pdratiosb = ROOT.TPad("pdratiosb","SB ratio",0,0,0.33,.25)
    pdratiosr = ROOT.TPad("pdratiosr","SR SB ratio",0.33,0,0.67,.25)
    pdratiosbs = ROOT.TPad("pdratiosbs","scaled ratio",0.67,0,1,.25)
    label1 = ROOT.TPaveText(.45,.35,.85,.45,"NBNDC")
    label2 = ROOT.TPaveText(.45,.35,.85,.45,"NBNDC")
    label3 = ROOT.TPaveText(.35,.35,.9,.45,"NBNDC")
    label1.AddText("SB MC, SB Data")
    label2.AddText("SR MC, SB Data")
    label3.AddText("SR MC, SB Data")
    label3.AddText("scaled to SR MC Yield")
    label1.SetFillColor(0)
    label2.SetFillColor(0)
    label3.SetFillColor(0)
    

    tc.cd()
    pdsbonly.Draw()
    pdsbonly.cd()
    hsbkgsb.Draw("HIST")
    xaxsb = hsbkgsb.GetXaxis()
    yaxsb = hsbkgsb.GetYaxis()
    xaxsb.SetTitle("M_{Z'}")
    xaxsb.SetTitleSize(0.05)
    xaxsb.SetLabelSize(0.035)
    yaxsb.SetTitle("Events / 200 GeV")
    yaxsb.SetTitleSize(0.05)
    yaxsb.SetLabelSize(0.04)
    yaxsb.SetLabelOffset(0.025)
    CMS_lumi.CMS_lumi(pdsbonly,4,13)
    hdat.Draw("same,e1")
    legsb.Draw()
    label1.Draw()
    pdsbonly.Update()

    tc.cd()
    pdratiosb.Draw()
    pdratiosb.cd()
    hdivsb.Draw()
    ratline.Draw()

    tc.cd()
    pdsbassr.Draw()
    pdsbassr.cd()
    hsbkgsr.Draw("HIST")
    xaxsr = hsbkgsr.GetXaxis()
    yaxsr = hsbkgsr.GetYaxis()
    xaxsr.SetTitle("M_{Z'}")
    xaxsr.SetTitleSize(0.05)
    xaxsr.SetLabelSize(0.035)
    yaxsr.SetTitle("Events / 200 GeV")
    yaxsr.SetTitleSize(0.05)
    yaxsr.SetLabelSize(0.04)
    yaxsr.SetLabelOffset(0.025)
    CMS_lumi.CMS_lumi(pdsbassr,4,13)
    hdat.Draw("same,e1")
    legsr.Draw()
    label2.Draw()
    pdsbassr.Update()

    tc.cd()
    pdratiosr.Draw()
    pdratiosr.cd()
    hdivsr.Draw()
    ratline.Draw()

    tc.cd()
    pdsbscaled.Draw()
    pdsbscaled.cd()
    hsbkgsr.Draw("HIST")
    CMS_lumi.CMS_lumi(pdsbscaled,4,13)
    hdats.Draw("same,el")
    legsr.Draw()
    label3.Draw()
    pdsbscaled.Update()

    tc.cd()
    pdratiosbs.Draw()
    pdratiosbs.cd()
    hdivsbs.Draw()
    ratline.Draw()

    gofplots = go.makeOutFile('Run2_'+yearstr,'gofplots_'+systr,'.png',"100","300","75","8E-10")
    tc.SaveAs(gofplots)

    gofroot = go.makeOutFile('Run2_'+yearstr,'gofplots_'+systr,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))

    outf = ROOT.TFile(gofroot,"recreate")
    hdats.Write()
    outf.Close()
    
    
