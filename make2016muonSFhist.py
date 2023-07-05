import ROOT
import uproot as up

#pathBCDEF = 'leptonsf/Run2016BCDEF_muon_SF_ID.root'
#pathGH    = 'leptonsf/Run2016GH_muon_SF_ID.root'

#pathBCDEF = 'leptonsf/RunBCDEF_SF_ID.root'
#pathGH    = 'leptonsf/RunGH_SF_ID.root'

#f1 = ROOT.TFile(pathBCDEF,"read")
#f2 = ROOT.TFile(pathGH,"read")
#d1 = f1.Get("NUM_TightID_DEN_genTracks_eta_pt")
#d2 = f1.Get("NUM_TightID_DEN_genTracks_eta_pt")


#d1 = f1.Get("NUM_TightID_DEN_genTracks_pt_abseta")
#d2 = f2.Get("NUM_TightID_DEN_genTracks_pt_abseta")
#lumibcdef = 20.
#lumigh = 16.

#scaledbcdef = d1.Clone()
#scaledbcdef.Scale(lumibcdef)
#scaledgh = d2.Clone()
#scaledgh.Scale(lumigh)
#scaledbcdef.Add(scaledgh)
#scaledbcdef.Scale(1/(lumibcdef+lumigh))

#Just to flip generic hists
pathBCDEF = 'leptonsf/Ele115orEleIso32orPho200_SF_2018.root'

fforchar = up.open(pathBCDEF)
fog = ROOT.TFile(pathBCDEF)

#need the SF hist and the efficieny hist
#listofhists = ['NUM_Mu50_TkMu100_DEN_TightID_abseta_pt','NUM_Mu50_TkMu100_DEN_TightID_abseta_pt_efficiencyMC']
listofhists = ['SF_TH2F']
histstowrite = []
for hname in listofhists:

    hforchar = fforchar[hname]
    npchar = hforchar.to_numpy()
    print(npchar)


    htoflip = fog.Get(hname)
    xbins = htoflip.GetNbinsX()
    ybins = htoflip.GetNbinsY()

    hflip = ROOT.TH2D(hname,hname,ybins,npchar[2],xbins,npchar[1])

    #print("The orignal X bins: ",scaledbcdef.GetNbinsX())
    print("The orignal X bins: ",htoflip.GetNbinsX())
    print("The new X bins:     ",hflip.GetNbinsX())
    print("The orignal Y bins: ",htoflip.GetNbinsY())
    print("The new Y bins:     ",hflip.GetNbinsY())
    
    for xbin in range(htoflip.GetNbinsX()+1):
        for ybin in range(htoflip.GetNbinsY()+1):
            sf = htoflip.GetBinContent(xbin,ybin)
            sferr = htoflip.GetBinError(xbin,ybin)
            hflip.SetBinContent(ybin,xbin,sf)
            hflip.SetBinError(ybin,xbin,sferr)

    histstowrite.append(hflip)

fsave = ROOT.TFile('leptonsf/Run2018_electron_SF_TRIGGER.root','recreate')

for hist in histstowrite:
    hist.Write()



