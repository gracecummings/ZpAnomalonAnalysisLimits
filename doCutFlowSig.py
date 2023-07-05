import subprocess
import gecorg_test as go
import argparse
import ROOT
import numpy as np

def gatherYields(topiaryf,bkgbinf):
    topiary       = ROOT.TFile(topiaryf)        
    scale         = bkgbinf["scale"]
    evntstot      = float(str(bkgbinf["tfile"].Get('hnevents').GetString()))*scale
    evntspassmet  = float(str(bkgbinf["tfile"].Get('hnevents_pMET').GetString()))*scale#sdmass too
    evntspasszpt  = float(str(bkgbinf["tfile"].Get('hnevents_pZ').GetString()))*scale
    evntspasshpt  = float(str(bkgbinf["tfile"].Get('hnevents_ph').GetString()))*scale
    evntspassbtag = float(str(bkgbinf["tfile"].Get('hnevents_btag').GetString()))*scale
    evntspasssb   = float(str(bkgbinf["tfile"].Get('hnevents_sb').GetString()))*scale
    evntsinskim   = topiary.Get('hnskimed').GetBinContent(1)*scale
    evntspasstrig = topiary.Get('htrigpass').GetBinContent(1)*scale
    evntspasszrec = topiary.Get('hZpass').GetBinContent(1)*scale
    evntspasshrec = topiary.Get('hHpass').GetBinContent(1)*scale

    return evntstot,evntsinskim,evntspasstrig,evntspasszrec,evntspasshrec,evntspassmet,evntspasszpt,evntspasshpt,evntspassbtag,evntspasssb

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-x","--xsec", type=float,help = "desired siganl cross section in fb")
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    parser.add_argument("-date","--date",help="date folder with plots to stack")
    args = parser.parse_args()

    #Get command line parameters
    sig_xsec      = args.xsec
    zptcut        = args.zptcut
    hptcut        = args.hptcut
    metcut        = args.metcut
    btagwp        = args.btagwp

    zp1500nd600ns200top = "pfMET_nominal/Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8_topiary_mumu_systnominal_Zptcut0.0_Hptcut250.0_metcut0.0_btagwp0.0.root"
    zp1500nd600ns200up = "post_muonsf_nominal/Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8_upout_signalr_systnominal_btagnom_muidnom_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root"
    zp2000nd800ns200up = "post_muonsf_nominal/Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8_upout_signalr_systnominal_btagnom_muidnom_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root"
    zp2000nd800ns200top = "pfMET_nominal/Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8_topiary_mumu_systnominal_Zptcut0.0_Hptcut250.0_metcut0.0_btagwp0.0.root"
    zp3000nd400ns200up = "post_muonsf_nominal/Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8_upout_signalr_systnominal_btagnom_muidnom_DeepMassDecorrelTagZHbbvsQCD_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root"
    zp3000nd400ns200top = "pfMET_nominal/Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8_topiary_mumu_systnominal_Zptcut0.0_Hptcut250.0_metcut0.0_btagwp0.0.root"
    sigsamps = [zp1500nd600ns200up,zp2000nd800ns200up,zp3000nd400ns200up]
    sigtops  = [zp1500nd600ns200up,zp2000nd800ns200top,zp3000nd400ns200top]

    sigcols  = go.colsFromPalette(sigsamps,ROOT.kLake)
    siginfo  = go.prepSig(sigsamps,sigcols,sig_xsec,101.27)#gathers xs scales
    skimstr  = "Events with \( > 0\) fat jets and \( Z(\ell^{+}\ell^{-}) cand \)"
    trigstr  = "Events Pass trigger"
    zrecostr = "Events with \( 70 < m_{\ell\ell} < 110 \)"
    hrecostr = "Events Pass $H$ reco "
    metstr = "Events Passing \(MET > "+str(metcut)+"\)"
    zptstr = "Events Passing \(Z p_{T} > "+str(zptcut)+"\)"
    hptstr = "Events Passing Fat Jet \(p_{T} > "+str(hptcut)+"\)"
    btgstr = "Events Passing {0} btag WP".format(str(btagwp))
    sbstr  = "Events in side band"
    
    
    cfdict = {metstr:{},zptstr:{},hptstr:{},btgstr:{},sbstr:{},skimstr:{},trigstr:{},zrecostr:{},hrecostr:{}}
    
    totorig = 0
    totskim = 0
    tottrig = 0
    totzrec = 0
    tothrec = 0
    totymet = 0
    totyzpt = 0
    totyhpt = 0
    totybtg = 0
    totysb  = 0
    
    for sig in siginfo:
        topi = 'str'
        for top in sigtops:
            if sig["name"] in top:
                topi= top

        origevnts = 0
        yieldskim = 0
        yieldtrig = 0
        yieldzrec = 0
        yieldhrec = 0
        yieldmet  = 0
        yieldzpt  = 0
        yieldhpt  = 0
        yieldbtag = 0
        yieldsb   = 0
        #cfdict

        topiary       = ROOT.TFile(topi)
        print(sig["name"])
        
        scale         = sig["scale"]
        print(scale)
        evntstot17      = float(str(sig["tfile"].Get('hnevents').GetString()))*scale
        evntspassmet17  = float(str(sig["tfile"].Get('hnevents_pMET').GetString()))*scale#sdmass too
        evntspasszpt17  = float(str(sig["tfile"].Get('hnevents_pZ').GetString()))*scale
        evntspasshpt17  = float(str(sig["tfile"].Get('hnevents_ph').GetString()))*scale
        evntspassbtag17 = float(str(sig["tfile"].Get('hnevents_btag').GetString()))*scale
        evntspasssb17   = float(str(sig["tfile"].Get('hnevents_sr').GetString()))*scale##sr
        #evntsinskim17   = topiary.Get('hnskimed').GetBinContent(1)*scale
        #evntspasstrig17 = topiary.Get('htrigpass').GetBinContent(1)*scale
        #evntspasszrec17 = topiary.Get('hZpass').GetBinContent(1)*scale
        #evntspasshrec17 = topiary.Get('hHpass').GetBinContent(1)*scale

        #cfdict[skimstr][sig["name"]] = round(evntsinskim17,2)
        #cfdict[trigstr][sig["name"]] = round(evntspasstrig17,2)
        #cfdict[zrecostr][sig["name"]] = round(evntspasszrec17,2)
        #cfdict[hrecostr][sig["name"]] = round(evntspasshrec17,2)
        cfdict[metstr][sig["name"]] = round(evntspassmet17,2)
        cfdict[zptstr][sig["name"]] = round(evntspasszpt17,2)
        cfdict[hptstr][sig["name"]] = round(evntspasshpt17,2)
        cfdict[btgstr][sig["name"]] = round(evntspassbtag17,2)
        cfdict[sbstr][sig["name"]] = round(evntspasssb17,2)

        print(sig["name"])

    cutFlowTableTex = go.makeOutFile("ZpAnomalon_mumu_sig","cutflow",'.tex',str(int(zptcut)),str(int(hptcut)),str(int(metcut)),str(btagwp).split('.')[-1]+'E-10')
    cftab = open(cutFlowTableTex,"w")
    cftab.write(r'\clearpage')
    cftab.write('\n')
    cftab.write(r'\begin{sidewaystable}[h!]\centering')
    cftab.write('\n')
    cftab.write(r'\begin{tabular}{l | c | c | c | c | c }')#need to make dynamic
    cftab.write('\n')
    cftab.write("\hline\hline\n")
    cftab.write(r'cut description & mzp1500    & mzp2000      & mzp3000     \\')
    cftab.write('\n')
    cftab.write(r'                & mnd600mns200 & mnd800mns200 & mnd400mns200 \\')
    cftab.write('\n')
    cftab.write("\hline\n")
    cftab.write(metstr+' & '+str(cfdict[metstr]['Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'])+' & '+str(cfdict[metstr]['Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'])+' & '+str(cfdict[metstr]['Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write(zptstr+' & '+str(cfdict[zptstr]['Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'])+' & '+str(cfdict[zptstr]['Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'])+' & '+str(cfdict[zptstr]['Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write(hptstr+' & '+str(cfdict[hptstr]['Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'])+' & '+str(cfdict[hptstr]['Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'])+' & '+str(cfdict[hptstr]['Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write(btgstr+' & '+str(cfdict[btgstr]['Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'])+' & '+str(cfdict[btgstr]['Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'])+' & '+str(cfdict[btgstr]['Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write("\hline\n")
    cftab.write('Events in Signal Region & '+str(cfdict[sbstr]['Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'])+' & '+str(cfdict[sbstr]['Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'])+' & '+str(cfdict[sbstr]['Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8']))
    cftab.write(r'\\')
    cftab.write('\n')
    cftab.write("\end{tabular}\n")
    cftab.write("\caption{Signal yields with luminosity scaling for 2017 and 2018. XS 10 fb}\n")
    cftab.write("\label{tab:sigeff}\n")
    cftab.write("\end{sidewaystable}\n")
    cftab.write(r'\clearpage')
    cftab.close()
        
