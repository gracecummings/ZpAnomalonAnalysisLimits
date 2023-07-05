import ROOT
import uproot as up


pathrecosf = 'leptonsf/Run2018_electron_SF.root'
pathidsf   = 'leptonsf/Run2018_electron_SFID.root'

f1 = ROOT.TFile(pathrecosf,"read")
f2 = ROOT.TFile(pathidsf,"read")
d1 = f1.Get("EGamma_SF2D")
d2 = f2.Get("EGamma_SF2D")

combosf = d1.Clone()
#combosf.Add(d2)

recobinedgesf = up.open(pathrecosf)
idbinegdesf  = up.open(pathidsf)

recosfup = recobinedgesf["EGamma_SF2D"]
idsfup   = idbinegdesf["EGamma_SF2D"]

npreco = recosfup.to_numpy()
npid   = idsfup.to_numpy()

print("reco hist info")
print(npreco)
print("id sf hist info")
print(npid)

#for bx in range(d1.GetNbinsX()+2):
#    print("The reco sf hist x axis low edge is ",d1.GetXaxis().GetBinLowEdge(bx))
#    print("The ID   sf hist x axis low edge is ",d2.GetXaxis().GetBinLowEdge(bx))
#    for by in range(d1.GetNbinsY()+2):
#        print("The reco sf hist y axis low edge is ",d1.GetYaxis().GetBinLowEdge(by))
#        print("The ID   sf hist y axis low edge is ",d2.GetYaxis().GetBinLowEdge(by))
#        print("Bin number {0} {1}".format(bx,by))
#        print("   Bin Error in reco sf hist ",d1.GetBinError(bx,by))
#        print("   Bin Error in id sf hist   ",d2.GetBinError(bx,by))

