import glob
import ROOT
import tdrstyle
import CMS_lumi
import configparser
import gecorg_test as go

def orderDY(histFile):
    s0 = histFile.split("topWithSummaryPlots/")[-1]
    s1 = s0.split("to")[0]
    s2 = s1.split("HT-")[1]
    return int(s2)       

def getDatahist(datadict,t0,hname,years=[16,17,18]):
    h0 = t0.Get(hname).Clone()
    h0.Reset("ICESM")
    for year in years:
        files = datadict[year]
        for f in files:
            #print(f)
            tf = ROOT.TFile(f)
            h = tf.Get(hname)
            hdat = h.Clone()
            h0.Add(hdat)
    return h0
        

def getAddedHist(samp,bkgdict,t0,hname,config,years = [16,17,18]):
    xspairs = config.items(samp)
    h0 = t0.Get(hname).Clone()
    h0.Reset("ICESM")
    for year in years:
        if year == 16:
            lumi = 35.9
        if year == 17:
            lumi = 41.53
        if year == 18:
            lumi = 59.74
                
        files = bkgdict[samp][year]
        if "DYJetsToLL" in samp:
            files.sort(key = orderDY)
        if "TT" in samp:
            files.sort()

        for i,f in enumerate(files):
            fparts = f.split("/")
            name = fparts[-1]
            #print("     adding ",name)
            tf = ROOT.TFile(f)
            numevents = tf.Get('hnorigevnts').GetBinContent(1)
            xs = float(xspairs[i][1].split()[0])*1000#Into Femtobarn
            scale = go.findScale(numevents,xs,lumi)
            h = tf.Get(hname)
            hscaled = h.Clone()
            hscaled.Scale(scale)
            h0.Add(hscaled)
                
    return h0

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
    dy16 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Summer16v3.DYJets*topiary*")
    dy17 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Fall17.DYJets*topiary*")
    dy18 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Autumn18.DYJets*topiary*")
    tt16 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Summer16v3.TTTo*topiary*")
    tt17 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Fall17.TTTo*topiary*")
    tt18 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Autumn18.TTTo*topiary*")
    wz16 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Summer16v3.WZ*topiary*")
    wz17 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Fall17.WZ*topiary*")
    wz18 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Autumn18.WZ*topiary*")
    zz16 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Summer16v3.ZZ*topiary*")
    zz17 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Fall17.ZZ*topiary*")
    zz18 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Autumn18.ZZ*topiary*")
    dt16 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Run2016*topiary*")
    dt17 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Run2017*topiary*")
    dt18 = glob.glob("analysis_output_ZpAnomalon/topWithSummaryPlots/Run2018*topiary*")

    bkgs = {"DYJetsToLL":{16:dy16,17:dy17,18:dy18},
            "TT":{16:tt16,17:tt17,18:tt18},
            "WZTo2L2Q":{16:wz16,17:wz17,18:wz18},
            "ZZTo2L2Q":{16:zz16,17:zz17,18:zz18}}
    data = {16:dt16,17:dt17,18:dt18}

    #Select Plotting years and region
    years = [16,17,18]
    yearstr = yearFormatter(years)

    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('xsects_2017.ini')#2017 and 2018dy+jets xs same
    config.read_file(fp)


    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = lumiFormatter(years)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Simulation Preliminary"

    
    hists = ['hjetpt','hjeteta','hjetphi']
    titles = {'hjetpt':"all fat jet pT, no DY k-factors",'hjeteta':"all fat jet eta",'hjetphi':"all fat jet phi"}

    testfile = ROOT.TFile(bkgs["DYJetsToLL"][18][0])
    for hist in hists:
        print(hist)
        print("  DY")
        hdy = getAddedHist("DYJetsToLL",bkgs,testfile,hist,config,years)
        print("  TT")
        htt = getAddedHist("TT",bkgs,testfile,hist,config,years)
        print("  WZ")
        hwz = getAddedHist("WZTo2L2Q",bkgs,testfile,hist,config,years)
        print("  ZZ")
        hzz = getAddedHist("ZZTo2L2Q",bkgs,testfile,hist,config,years)

        hbkgs = [hdy,htt,hwz,hzz]
        hsbkg = ROOT.THStack()
        leg = ROOT.TLegend(0.65,0.55,0.9,0.8)
        for i,h in enumerate(hbkgs):
            h.SetFillColor(bkgcols[i])
            h.SetLineColor(bkgcols[i])
            hsbkg.Add(h)
            leg.AddEntry(h,bkgnames[i],"f")

        htotbkg = hdy.Clone()
        htotbkg.Add(htt)
        htotbkg.Add(hwz)
        htotbkg.Add(hzz)

        print("  data")
        hdat = getDatahist(data,testfile,hist,years)
        hdat.SetBinErrorOption(1)
        hdat.SetMarkerStyle(8)
        hdat.SetMarkerSize(0.7)
        hdat.SetMarkerColor(ROOT.kBlack)

        leg.AddEntry(hdat,'Data',"pe")

        hdiv = hdat.Clone()
        hdiv.Divide(hdat,htotbkg)
        hdiv.SetMarkerStyle(8)
        hdiv.SetMarkerSize(0.7)
        hdiv.SetMarkerColor(ROOT.kBlack)
        


        #Plotting itself
        print("  Plotting")
        ratline = ROOT.TLine(htotbkg.GetBinLowEdge(1),1,htotbkg.GetBinLowEdge(htotbkg.GetNbinsX())+htotbkg.GetBinWidth(htotbkg.GetNbinsX()),1)
        tc = ROOT.TCanvas("tc",hist,600,800)
        p1 = ROOT.TPad("p1","stack_"+hist,0,0.3,1.0,1.0)
        p1.SetLeftMargin(0.15)
        p1.SetRightMargin(0.05)
        p2 = ROOT.TPad("p2","signif_"+hist,0,0.0,1.0,0.3)
        p2.SetRightMargin(.05)
        p2.SetLeftMargin(0.15)
        p2.SetBottomMargin(0.25)
        p2.SetTopMargin(0.05)

        #Prepare first pad for stack
        tc.Draw()
        tc.cd()
        
        p1.Draw()
        p1.cd()
        p1.SetLogy()
        #Draw the stack
        hsbkg.SetMaximum(htotbkg.GetMaximum()*100)
        hsbkg.Draw("HIST")#add PFC for palette drawing
        hsbkg.GetXaxis().SetTitle(titles[hist])
        hsbkg.GetXaxis().SetTitleSize(0.05)
        hsbkg.GetXaxis().SetTitleOffset(1.1)
        hsbkg.GetXaxis().SetLabelSize(0.04)
        hsbkg.GetYaxis().SetTitle("Events")
        hsbkg.GetYaxis().SetTitleSize(0.05)
        hsbkg.GetYaxis().SetLabelSize(0.04)
        CMS_lumi.CMS_lumi(p1,4,13)
        hdat.Draw("histsame,pe")
        leg.Draw()

        tc.cd()
        p2.Draw()
        p2.cd()
        hdiv.Draw()
        hdiv.GetXaxis().SetTitle("bin center")
        hdiv.GetXaxis().SetTitleSize(0.11)
        hdiv.GetXaxis().SetTitleOffset(0.7)
        hdiv.GetXaxis().SetLabelSize(0.1)
        hdiv.GetYaxis().SetTitle("data/MC")
        hdiv.GetYaxis().SetTitleSize(0.11)
        hdiv.GetYaxis().SetTitleOffset(.45)
        hdiv.GetYaxis().SetLabelSize(0.1)
        hdiv.GetYaxis().SetLabelOffset(0.02)
        hdiv.GetYaxis().SetNdivisions(503)
        hdiv.SetMinimum(0.)
        hdiv.SetMaximum(2.)
        ratline.Draw()
        
        tc.cd()

        pngname = go.makeOutFile(hist,'alljets','.png','ThereIsAZ','ThereIsAJet','0','0')
        tc.SaveAs(pngname)
        

        
