import ROOT
import gecorg_test as go
import glob
import numpy as np

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
    newttf = 'analysis_output_ZpAnomalon/2024-02-08/ttbar_shenanigan_fits_to_data_ratio_161718_signalr_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root'
    
    #get new ttbar
    esttf = ROOT.TFile.Open(newttf,'read')
    hest = esttf.Get('hemuexp')
    hest.SetDirectory(0)
    hestlimit = go.newNameAndStructure(hest,"TT",1,limrangelow,limrangehigh)
    hestlimit = hestlimit.Rebin(len(newbinedges)-1,"TT",newbinedges)
    print('Total new ttbar yield: ',hestlimit.Integral())
    

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
        tf.Close()
        outname = go.makeOutFile('Run2_161718_ZllHbbMET','unblindttshenanigans_mumu_'+sigpnt,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        outf = ROOT.TFile(outname,"recreate")
        hestlimit.Write()
        for h in hl:
            h.Write()
        outf.Close()
