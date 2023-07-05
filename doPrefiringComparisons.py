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

    #Select Plotting years and region
    years = [16,17]
    if year:
        years = [int(year)]
    yearstr = yearFormatter(years)
    reg = regionFormatter(regname)
    if plot_data:
        print('Plotting the data!')
    else:
        print('You know you are not plotting data, right?')

    bkgsnom = go.backgrounds('mumu_2022-07-06_ProperMuIDTriSF_METXY',zptcut,hptcut,metcut,btagwp,'systnominal_kfnom_btagnom_muidnom_mutrignom_elidnom_elreconom')
    bkgspre = go.backgrounds('mumu_2022-07-06_ProperMuIDTriSF_METXY_prefire',zptcut,hptcut,metcut,btagwp,'systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco')

    #Gather plots
    testyear = years[0]#picks first year in list, so desired year if only one
    testfile = bkgsnom.bkgs["DYJetsToLL"][testyear][reg][0][0]#stacked plots should always have DY
    testtfile = ROOT.TFile(testfile)

    bkgnames = ["DYJetsToLL","TT","WZTo2L2Q","ZZTo2L2Q"]
    bkgcols  = go.colsFromPalette(bkgnames,ROOT.kLake)
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = lumiFormatter(years)
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Simulation Preliminary"

    nomplotssr = []
    nomplotssb = []
    nomplotstr = []
    pfplots = []

    for bkg in bkgnames:
        nomplotssr.append(gatherIndividualPlots(bkgsnom,testfile,'h_zp_jignom',bkg,'sr',years))
        nomplotssb.append(gatherIndividualPlots(bkgsnom,testfile,'h_zp_jignom',bkg,'sb',years))
        nomplotstr.append(gatherIndividualPlots(bkgsnom,testfile,'h_zp_jignom',bkg,'tr',years))
        pfplotssr.append(gatherIndividualPlots(bkgspre,testfile,'h_zp_jignom',bkg,'sr',years))
        pfplotssb.append(gatherIndividualPlots(bkgspre,testfile,'h_zp_jignom',bkg,'sb',years))
        pfplotstr.append(gatherIndividualPlots(bkgspre,testfile,'h_zp_jignom',bkg,'tr',years))
