import glob
import ROOT as ROOT
import gecorg_test as go

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

def applyStatsUncToSignal(hist,errseries):
    for ibin in range(hist.GetNbinsX()+1):
        if ibin == 0:
            continue
        else:
            binerr = errseries[ibin-1]
            hist.SetBinError(ibin,binerr)
    return hist

def makeSignalInfoDict(sigclass, region,sigxs):
    sigs = sigclass.getPreppedSig(region,sigxs)
    sigdict = {}
    for sig in sigs:
        sigdict[sig["name"]] = sig
    return sigdict


if __name__=='__main__':
    #Gather the signals
    syststr = 'systnominal_hem_kf_btag_muid_mutrig_eltrig_elid_elreco'
    sigs = go.signal('mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref/ns1sigs',100.0,300.0,75.0,0.8,1,[16,17,18],syststr)
    siginfo = sigs.getPreppedSig('sr',1)
    signom = makeSignalInfoDict(sigs,'sr',1)
    sigcolors = go.colsFromPalette(siginfo,ROOT.kCMYK)
    nomsigs = [s["name"] for s in siginfo]

    #make the outfile
    outfname = go.makeOutFile('Run2_161718_ZllHbbMET',syststr,'.root','100','300','75','8E-10')
    outf = ROOT.TFile(outfname,"recreate")
    
    for s in nomsigs:
        #gathering names
        name = signom[s]["name"]
        print(name)
        signame = "holder"
        #if "NS1" in name:
        #    continue
        if "Tune" in name:
            strippedname = name.split("_Tune")[0]
            signame = strippedname.replace("-","")
        else:
            signame = name.replace("-","")

        #The histogtrams
        hsigori = signom[s]["tfile"].Get("h_zp_jigm")
        #The lumiu scaled
        hsig = hsigori.Clone()
        hsig.Scale(signom[s]["scale"])
        hsig = applyStatsUncToSignal(hsig,signom[s]["errdf"]["h_zp_jigm"]*signom[s]["scale"])
        hsigf = hsig.Clone()
        hsig = newNameAndStructure(hsig,signame+"_lumixs_combiner",2,1400,3000)#combine ready
        hsigf.Rebin(2)
        hsigf.SetName(signame+"_lumixs")#lumi scaled whole distribution
        #The straight
        hog = hsigori.Clone()
        hog = applyStatsUncToSignal(hog,signom[s]["errdf"]["h_zp_jigm"])
        hog.Rebin(2)
        hog.SetName(signame)
        #Write it
        outf.cd()
        hsig.Write()
        hsigf.Write()
        hog.Write()
