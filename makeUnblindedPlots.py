import ROOT
import tdrstyle
import CMS_lumi
import gecorg_test as go
import numpy as np

def convertToPhysicalHist(hemp,hbad):
    h = hemp.Clone()
    if hemp.GetNbinsX() == hbad.GetNbinsX():
        for i in range(hbad.GetNbinsX()+1):
            h.SetBinContent(i,hbad.GetBinContent(i))
            h.SetBinError(i,hbad.GetBinError(i))
    else:
        print("Bad histograms, needs attention")
    return h

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


    #Select Plotting years and region
    years = [16,17,18]
    yearstr = yearFormatter(years)

    postfit = False
    isLog   = False

    signame = "Zp5500ND1800NS200"
    f = ROOT.TFile('unblindedDatacardHolder/fitDiagnosticsRun2_161718_ZllHbbMET_datacard_unblind_mumu_'+signame+'.root')
    f2 = ROOT.TFile('unblindedDatacardHolder/Run2_161718_ZllHbbMET_unblind_mumu_'+signame+'_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    htest = f2.Get("TT")
    empty = htest.Clone()
    empty.SetDirectory(0)
    f2.Close()
    empty.Reset('ICESM')

    hprett = convertToPhysicalHist(empty,f.Get('shapes_prefit/'+signame+'_mumu/TT'))
    hprevv = convertToPhysicalHist(empty,f.Get('shapes_prefit/'+signame+'_mumu/VV'))
    hpredy = convertToPhysicalHist(empty,f.Get('shapes_prefit/'+signame+'_mumu/DY'))
    hpretot = convertToPhysicalHist(empty,f.Get('shapes_prefit/total_background'))
    gdat = f.Get('shapes_prefit/total_data')

    hposttt = convertToPhysicalHist(empty,f.Get('shapes_fit_s/'+signame+'_mumu/TT'))
    hpostvv = convertToPhysicalHist(empty,f.Get('shapes_fit_s/'+signame+'_mumu/VV'))
    hpostdy = convertToPhysicalHist(empty,f.Get('shapes_fit_s/'+signame+'_mumu/DY'))
    hposttot = convertToPhysicalHist(empty,f.Get('shapes_fit_s/total_background'))

    hpost0sigtt = convertToPhysicalHist(empty,f.Get('shapes_fit_b/'+signame+'_mumu/TT'))
    hpost0sigvv = convertToPhysicalHist(empty,f.Get('shapes_fit_b/'+signame+'_mumu/VV'))
    hpost0sigdy = convertToPhysicalHist(empty,f.Get('shapes_fit_b/'+signame+'_mumu/DY'))
    hpost0sigtot = convertToPhysicalHist(empty,f.Get('shapes_fit_b/total_background'))

    bkgs = {"prefit":[hpredy,hprett,hprevv,hpretot],
            "postfit":[hpostdy,hposttt,hpostvv,hposttot],
            "posfitbkgonly":[hpost0sigdy,hpost0sigtt,hpost0sigvv,hpost0sigtot],
}

    #Print some debug info
    #ROOT.gROOT.SetBatch(True)
    #for i in range(hpost0sigtot.GetNbinsX()+1):
    #    print("Comparing bin ",i)
    #    print("  combine bin low edge: ",hpost0sigtot.GetBinLowEdge(i))
    #    print("   manual bin low edge: ",empty.GetBinLowEdge(i))
    #    print("     combine bin width: ",hpost0sigtot.GetBinWidth(i))
    #    print("      manual bin width: ",empty.GetBinWidth(i))

    #make an hdat
    hdat = hpost0sigtot.Clone()
    hdat.Reset("ICESM")
    hdat.SetMarkerStyle(8)
    hdat.SetMarkerSize(0.7)
    hdat.SetMarkerColor(ROOT.kBlack)

    datx = ROOT.TVector(gdat.GetN())
    daty = ROOT.TVector(gdat.GetN())
    daterrxl = ROOT.TVector(gdat.GetN())
    daterrxh = ROOT.TVector(gdat.GetN())
    daterryl = ROOT.TVector(gdat.GetN())
    daterryh = ROOT.TVector(gdat.GetN())
    for i in range(gdat.GetN()):
        gx = gdat.GetPointX(i)
        gy = gdat.GetPointY(i)
        hdat.SetBinContent(i+1,gy)
        daty[i] = gy
        daterrxl[i] = gdat.GetErrorXlow(i)
        daterrxh[i] = gdat.GetErrorXhigh(i)
        daterryl[i] = gdat.GetErrorYlow(i)
        daterryh[i] = gdat.GetErrorYhigh(i)
        datx[i] = hdat.GetBinCenter(i+1)

    #print(datx)
    #print(daty)
    #print(daterryl)
    #print(daterryh)

    gdat2 = ROOT.TGraphAsymmErrors(datx,daty,daterrxl,daterrxh,daterryl,daterryh)
        
    #Make the plots
    #general setup
    bkgcols = go.colsFromPalette(bkgs.keys(),ROOT.kLake)
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = lumiFormatter(years)
    CMS_lumi.writeExtraText = 1
    for key in bkgs.keys():
        CMS_lumi.extraText = "Internal, "+key+" bkg"
        if key == "postfit":
            CMS_lumi.extraText = "Internal, "+signame+" "+key+" bkg"
        if "bkg" in key:
            CMS_lumi.extraText = "Internal, bkg only postfit"
        stkpadydims = [0.3,1.]
        ratpadydims = [0.0,0.3]
        #tcanvasdims = [600,560]#no ratio pad
        tcanvasdims = [600,800]
        
        #ratio plot
        hdiv = hdat.Clone()
        hdiv.Divide(bkgs[key][3])
        

        #Colors
        bkgs[key][0].SetFillColor(bkgcols[0])
        bkgs[key][0].SetLineColor(bkgcols[0])
        bkgs[key][1].SetFillColor(bkgcols[1])
        bkgs[key][1].SetLineColor(bkgcols[1])
        bkgs[key][2].SetFillColor(bkgcols[2])
        bkgs[key][2].SetLineColor(bkgcols[2])
        bkgs[key][3].SetLineColor(bkgcols[0])#make the total one invisible
        bkgs[key][3].SetFillColor(ROOT.kBlack-3)
        bkgs[key][3].SetFillStyle(3245)

        #for when this was a graph
        #hdat.SetBinErrorOption(1)
        gdat2.SetMarkerStyle(8)
        gdat2.SetMarkerSize(0.7)
        gdat2.SetMarkerColor(ROOT.kBlack)

        #Bkg Stack for Plotting
        hsbkg = ROOT.THStack("hsbkg","")
        hsbkg.Add(bkgs[key][2])
        hsbkg.Add(bkgs[key][1])
        hsbkg.Add(bkgs[key][0])
        hsbkg.SetMinimum(0)
        hsbkg.SetMaximum(30)
        if isLog:
            hsbkg.SetMinimum(0.01)
            hsbkg.SetMaximum(100)

        #Legend
        leg = ROOT.TLegend(0.5,0.45,0.88,0.80)
        leg.SetBorderSize(0)
        leg.AddEntry(bkgs[key][0],"DYJetsToLL","f")
        leg.AddEntry(bkgs[key][1],"TT","f")
        leg.AddEntry(bkgs[key][2],"VV","f")
        leg.AddEntry(gdat2,"Data")

        #Label
        #lab = ROOT.TPaveText(.15,.80,.4,.90,"NBNDC")
        #lab.AddText("Signal Region")
        #lab.AddText(key+" background")
        #lab.SetFillColor(0)

        #ratio line
        ratline = ROOT.TLine(hdat.GetBinLowEdge(1),1,hdat.GetBinLowEdge(hdat.GetNbinsX())+hdat.GetBinWidth(hdat.GetNbinsX()),1)

        #tCanvas
        tc = ROOT.TCanvas("tc"+key,"tc"+key,tcanvasdims[0],tcanvasdims[1])
        p1 = ROOT.TPad("p1","stack_"+key,0,stkpadydims[0],1.0,stkpadydims[1])
        p1.SetLeftMargin(0.15)
        p1.SetRightMargin(0.05)
        p2 = ROOT.TPad("p2","ratio_"+key,0,ratpadydims[0],1.0,ratpadydims[1])
        p2.SetRightMargin(.05)
        p2.SetLeftMargin(0.15)
        p2.SetBottomMargin(0.25)
        p2.SetTopMargin(0.05)
    
        #Prepare first pad for stack
        tc.Draw()
        p1.Draw()
        p1.cd()
        hsbkg.Draw("HIST")#add PFC for palette drawing
        hsbkg.GetXaxis().SetTitle("Jigsaw Mass Estimator Z'")
        hsbkg.GetXaxis().SetTitleSize(0.05)
        hsbkg.GetXaxis().SetTitleOffset(1.1)
        hsbkg.GetXaxis().SetLabelSize(0.04)
        hsbkg.GetYaxis().SetTitle("Events")
        hsbkg.GetYaxis().SetTitleSize(0.05)
        hsbkg.GetYaxis().SetLabelSize(0.04)
        CMS_lumi.CMS_lumi(p1,4,13)
        ROOT.gStyle.SetErrorX(0.5)#binwidth for the weird combine stuff
        bkgs[key][3].Draw("same,E2")
        gdat2.Draw('P')
        leg.Draw()
        #lab.Draw()
        p1.Update()
        tc.cd()
        p2.Draw()
        p2.cd()
        hdiv.Draw()
        hdiv.GetXaxis().SetTitle("bin center")
        hdiv.GetXaxis().SetTitleSize(0.11)
        hdiv.GetXaxis().SetTitleOffset(0.65)
        hdiv.GetXaxis().SetLabelSize(0.075)
        hdiv.GetYaxis().SetTitle("data/bkg est")
        hdiv.GetYaxis().SetTitleSize(0.11)
        hdiv.GetYaxis().SetTitleOffset(.45)
        hdiv.GetYaxis().SetLabelSize(0.08)
        hdiv.GetYaxis().SetLabelOffset(0.02)
        hdiv.GetYaxis().SetNdivisions(503)
        if "bkg" in key:
            hdiv.SetMinimum(0)
            hdiv.SetMaximum(2)
        else:
            hdiv.SetMinimum(-4)
            hdiv.SetMaximum(5)
        ratline.Draw()
        tc.cd()

        pngname = go.makeOutFile("Run2_ZllHbbMET_"+key,signame+'_unblinded_'+yearstr,'.png',"100","300","75","8E-10")
        tc.SaveAs(pngname)


