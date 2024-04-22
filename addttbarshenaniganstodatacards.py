import ROOT
import gecorg_test as go
import glob
import numpy as np

def makeTxtWithExtraDCardLines(info):
    prepCardName = go.makeOutFile('Run2_161718_ZllHbbMET','datacard_extrattshenanigan_mumu','.txt',str(100.0),str(300.0),str(75.0),str(0.8))
    card = open(prepCardName,"w")

    systs = info.keys()
    maxnamelength = max([len(x) for x in systs])
    for syst in systs:
        uprate = info[syst][0]
        dwnrate = info[syst][1]
        ratestr = "-"
        if abs(round(1-uprate,2))-abs(round(1-dwnrate,2)) == 0:
            ratestr = str(max(round(uprate,2),round(dwnrate,2)))
        else:
            ratestr = str(round(dwnrate,3))+"/"+str(round(uprate,3))
        
        nameoffset = maxnamelength - len(syst)
        
        ratestr = "{0} {1} {2} {3} {4} {5}\n".format(syst+" "*nameoffset,'lnN','-','-',ratestr,'-')
        shapestr = "{0} {1} {2} {3} {4} {5}\n".format(syst+" "*nameoffset,'shape','-','-','1','-')
        card.write(ratestr)
        card.write(shapestr)

    card.close()


if __name__=='__main__':
    #Years
    years   = [16,17,18]
    yearstr = go.yearFormatter(years)

    #Starting parameters
    zptcut  = '100.0'
    hptcut  = '300.0'
    metcut  = '75.0'
    btagwp  = '0.8'
    rebindiv = 2
    limrangelow = 1800
    limrangehigh = 10000
    newbinedges = np.array([1800.0,2000.0,2200.0,2400.0,2600.0,2800.0,10000.0])

    #import files
    print("Importing root files")
    combf = glob.glob('unblindedDatacardHolder/Run2_161718_ZllHbbMET_unblind_mumu*.root')
    #newttf = 'analysis_output_ZpAnomalon/2024-02-08/ttbar_shenanigan_fits_to_data_ratio_161718_signalr_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root'
    #newttf = 'analysis_output_ZpAnomalon/2024-03-11/ttbar_shenanigan_fits_to_data_ratio_161718_signalr_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root'
    newttf = 'analysis_output_ZpAnomalon/2024-04-01/ttbar_shenanigan_fits_to_data_ratio_161718_signalr_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root'
    
    #get new ttbar
    esttf = ROOT.TFile.Open(newttf,'read')
    hest = esttf.Get('TT_nom_fit')
    hest.SetDirectory(0)
    #hestlimit = go.newNameAndStructure(hest,"TT",1,limrangelow,limrangehigh)
    #hestlimit = hestlimit.Rebin(len(newbinedges)-1,"TT",newbinedges)
    print('Total new ttbar yield: ',hest.Integral())

    #get new ttbar systematics
    hpar0up = esttf.Get('TT_Par0_up')
    hpar0dn = esttf.Get('TT_Par0_dwn')
    hpar1up = esttf.Get('TT_Par1_up')
    hpar1dn = esttf.Get('TT_Par1_dwn')
    haltup  = esttf.Get('TT_AltFunctxExp_up')
    haltdn  = esttf.Get('TT_AltFunctxExp_dwn')

    hpar0up.SetDirectory(0)
    hpar0dn.SetDirectory(0)
    hpar1up.SetDirectory(0)
    hpar1dn.SetDirectory(0)
    haltup.SetDirectory(0)
    haltdn.SetDirectory(0)

    hest.SetName('TT')
    hpar0up.SetName('TT_TT_par0Up')
    hpar0dn.SetName('TT_TT_par0Down')
    hpar1up.SetName('TT_TT_par1Up')
    hpar1dn.SetName('TT_TT_par1Down')
    haltup.SetName('TT_TT_altfuncUp')
    haltdn.SetName('TT_TT_altfuncDown')
    
    par0up = go.getDeviatedOverNominal(hpar0up,hest)
    print(hpar0up.Integral())
    print(hest.Integral())
    par0dn = go.getDeviatedOverNominal(hpar0dn,hest)
    par1up = go.getDeviatedOverNominal(hpar1up,hest)
    par1dn = go.getDeviatedOverNominal(hpar1dn,hest)
    altup = go.getDeviatedOverNominal(haltup,hest)
    altdn = go.getDeviatedOverNominal(haltdn,hest)
    
    print("making lines to add to the datacards")
    info = {'TT_par0':[par0up,par0dn],'TT_par1':[par1up,par1dn],'TT_altfunc':[altup,altdn]}
    makeTxtWithExtraDCardLines(info)

    print("making new root files")
    for oldf in combf:
        sigpnt = oldf.split("_")[5]
        print(sigpnt)
        tf = ROOT.TFile.Open(oldf,'read')
        keys = tf.GetListOfKeys()
        hl = []
        for key in keys:
            hname = key.GetName()
            if "TT" in hname:
                continue
            else:
                hl.append(tf.Get(hname))
                hl[-1].SetDirectory(0)
        
        outname = go.makeOutFile('Run2_161718_ZllHbbMET','unblindttshenanigans_withsysts_mumu_'+sigpnt,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        outf = ROOT.TFile(outname,"recreate")
        hest.Write()
        hpar0up.Write()
        hpar0dn.Write()
        hpar1up.Write()
        hpar1dn.Write()
        haltup.Write()
        haltdn.Write()
        
        for h in hl:
            h.Write()
        outf.Write()
        outf.Close()
        tf.Close()
                
