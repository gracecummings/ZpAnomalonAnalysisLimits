import ROOT
import tdrstyle
import CMS_lumi
import gecorg_test as go

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

    signame = "Zp2000ND400NS200"
    f = ROOT.TFile('unblindedDatacardHolder/fitDiagnosticsRun2_161718_ZllHbbMET_datacard_unblind_mumu_'+signame+'.root')


    hprett = f.Get('shapes_prefit/'+signame+'_mumu/TT')
    hprevv = f.Get('shapes_prefit/'+signame+'_mumu/VV')
    hpredy = f.Get('shapes_prefit/'+signame+'_mumu/DY')
    hpretot = f.Get('shapes_prefit/total_background')
    gdat = f.Get('shapes_prefit/total_data')

    hposttt = f.Get('shapes_fit_s/'+signame+'_mumu/TT')
    hpostvv = f.Get('shapes_fit_s/'+signame+'_mumu/VV')
    hpostdy = f.Get('shapes_fit_s/'+signame+'_mumu/DY')
    hposttot = f.Get('shapes_fit_s/total_background')

    hpost0sigtt = f.Get('shapes_fit_b/'+signame+'_mumu/TT')
    hpost0sigvv = f.Get('shapes_fit_b/'+signame+'_mumu/VV')
    hpost0sigdy = f.Get('shapes_fit_b/'+signame+'_mumu/DY')
    hpost0sigtot = f.Get('shapes_fit_b/total_background')

    bkgs = {"prefit":[hpredy,hprett,hprevv,hpretot],
            "postfit":[hpostdy,hposttt,hpostvv,hposttot],
            "posfitbkgonly":[hpost0sigdy,hpost0sigtt,hpost0sigvv,hpost0sigtot],
}

    #Print some debug info
    #ROOT.gROOT.SetBatch(True)
    #make an hdat
    hdat = hpost0sigtot.Clone()
    hdat.Reset("ICESM")
    hdat.SetMarkerStyle(8)
    hdat.SetMarkerSize(0.7)
    hdat.SetMarkerColor(ROOT.kBlack)

    for i in range(gdat.GetN()):
        gx = gdat.GetPointX(i)
        gy = gdat.GetPointY(i)
        hdat.SetBinContent(i+1,gy)

    #Make the plots
    #general setup
    bkgcols = go.colsFromPalette(bkgs.keys(),ROOT.kLake)
    tdrstyle.setTDRStyle()
    CMS_lumi.lumi_13TeV = lumiFormatter(years)
    CMS_lumi.writeExtraText = 1

    for key in bkgs.keys():
        CMS_lumi.extraText = "Internal, "+key+" bkg"
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
        gdat.SetMarkerStyle(8)
        gdat.SetMarkerSize(0.7)
        gdat.SetMarkerColor(ROOT.kBlack)

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
        leg.AddEntry(gdat,"Data")

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
        gdat.Draw('P')
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

        pngname = go.makeOutFile("Run2_ZllHbbMET_"+key,'unblinded_'+yearstr,'.png',"100","300","75","8E-10")
        tc.SaveAs(pngname)


