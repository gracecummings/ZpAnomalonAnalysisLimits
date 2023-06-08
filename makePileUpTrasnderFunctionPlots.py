import ROOT
import gecorg_test as go

tf = ROOT.TFile("analysis_output_ZpAnomalon/2023-05-30/pileup_tranfers_recoquantToTrue_Zptcut_Hptcut_metcut_btagwp.root")

keys = tf.GetListOfKeys()
keynames = [key.GetName() for key in keys]

nvtxprofs = sorted([key for key in keynames if "h_nvtx_truenum_tprof" in key])
nvtxvtruh = sorted([key for key in keynames if "h_nvtx_truenum_hist" in key])

y16nvtx = sorted([key for key in nvtxprofs if "Summer16" in key])
y17nvtx = sorted([key for key in nvtxprofs if "Fall17" in key])
y18nvtx = sorted([key for key in nvtxprofs if "Autumn18" in key])

y16th2 = sorted([key for key in nvtxvtruh if "Summer16" in key])
y17th2 = sorted([key for key in nvtxvtruh if "Fall17" in key])
y18th2 = sorted([key for key in nvtxvtruh if "Autumn18" in key])

byyear = {"2016":[y16nvtx,y16th2],"2017":[y17nvtx,y17th2],"2018":[y18nvtx,y18th2]}

pairs = list(zip(nvtxvtruh,nvtxprofs))

tcs = []
legs = []
averages = []
ROOT.gStyle.SetOptStat(0)
#for pair in pairs:
#    nameh = pair[0].split("_h")[0]
#    namep = pair[1].split("_h")[0]
#    if nameh != namep:
#        continue
#    tcs.append(ROOT.TCanvas("tc"+nameh,"tc"+nameh,600,400))
#    h = tf.Get(pair[0])
#    p = tf.Get(pair[1])
#    h.GetXaxis().SetTitle("True Number of Vtx")
#    h.GetYaxis().SetTitle("NVtx Reco")
#    h.SetTitle(nameh)
#    p.SetLineColor(ROOT.kRed)
#    p.SetLineWidth(2)
#    h.Draw("COLZ")
#    p.Draw("E, same")

#    outname = go.makeOutFile("PileUpTransferFunction_",nameh,".png","","","","")
#    tcs[-1].SaveAs(outname)
    

for year in byyear.keys():
    print("Plotting the TProfiles for ",year)
    ls = byyear[year]
    tcs.append(ROOT.TCanvas("tc"+year,"tc"+year,1200,800))
    legs.append(ROOT.TLegend(0.4,0.15,0.85,0.45))
    cols = go.colsFromPalette(ls[1],ROOT.kCMYK)
    for i,plot in enumerate(ls[0]):
        h = tf.Get(plot)
        h.SetLineColor(cols[i])
        h.GetXaxis().SetTitle("True Number of Vtx")
        h.GetYaxis().SetTitle("NVtx Reco")
        h.SetTitle(year+" transfer functions")
        legs[-1].AddEntry(h,plot.split("_h")[0],"l")
        if i == 0:
            av = h.Clone()
            h.Draw()
        else:
            av.Add(h)
            h.Draw("SAME")

    av.SetLineColor(ROOT.kRed)
    legs[-1].AddEntry(av,"average relation","l")
    av.Draw("hist,SAME")
    legs[-1].Draw()
    av.SetDirectory(0)
    av.SetName("tprofile_allmc_"+year)
    averages.append(av)
    
    outname = go.makeOutFile("PileUpTransferFunctions_ALL_",year,".png","","","","")
    tcs[-1].SaveAs(outname)

    outrootname = go.makeOutFile("Run2_161718","pileUpTransferFunctions_Averages",".root","","","","")
    outf = ROOT.TFile(outrootname,"RECREATE")
    for p in averages:
        p.Write()
    outf.Close()

