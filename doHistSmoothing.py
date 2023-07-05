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

if __name__=='__main__':

    bkgs = go.backgrounds('mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref_unclmet','100.0','300.0','75.0','0.8','systunclup_hem_kf_btag_muid_mutrig_eltrig_elid_elreco')

    tf1 = ROOT.TFile(bkgs.f17wzsr[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")
    empty1 = empty.Clone()
    empty2 = empty.Clone()
    hzz  = bkgs.getAddedHist(empty1,"ZZTo2L2Q","sr","h_zp_jigm")
    hwz  = bkgs.getAddedHist(empty2,"WZTo2L2Q","sr","h_zp_jigm")
    hvv  = hzz.Clone()
    hvv.Add(hwz)
    #hvv.Rebin(2)

    hvv = newNameAndStructure(hvv,"VV",2,1800,10000)
    newbinedges = makeBinLowEdges(hvv,2800)
    hvv = hvv.Rebin(len(newbinedges)-1,"VV",newbinedges)

    nbins = hvv.GetNbinsX()
    binlow = hvv.GetBinLowEdge(1)
    binhigh = hvv.GetBinLowEdge(nbins+1)+hvv.GetBinWidth(1)

    smooth = hvv.Clone()
    smooth.Reset("ICESM")
    
    binconts = []
    for b in range(nbins):
        binconts.append(hvv.GetBinContent(b+1))
    arrbinconts = np.array(binconts)
    spec = ROOT.TSpectrum()
    spec.SmoothMarkov(arrbinconts,nbins,1)
    for b in range(nbins):
        smooth.SetBinContent(b+1,arrbinconts[b])


    hvv.SetLineColor(ROOT.kRed)
    smooth.SetLineColor(ROOT.kBlue)

    #hvv.Rebin(2)
    #smooth.Rebin(2)

    #hvv.GetXaxis().SetRangeUser(1400,10000)
    hvv.GetYaxis().SetRangeUser(0,0.6)
    #smooth.GetXaxis().SetRangeUser(1400,10000)
    tc = ROOT.TCanvas("tc","tc",800,600)
    tc.cd()
    hvv.Draw('hist')
    smooth.Draw('histsame')
    tc.SaveAs("markovchain_test_uncup_newrange.png")
