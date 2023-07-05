#this was a start of adding the years together, bit it is not clear if that is the better option.
import glob
import ROOT as ROOT
import gecorg_test as go
import configparser

if __name__=='__main__':
    #Gather the types of signals
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()

    signal_run2  = go.signal_run2(config.get('nominal','pathsignom'),'100.0','300.0','75.0','0.8','1.0',[16,17,18],config.get('nominal','strnom'))
    tf1 = ROOT.TFile(signal_run2.sig16sr[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")#creates an empty hist with same structure


    #make the outfile
    outfname = go.makeOutFile('Run2_161718_ZllHbbMET',"inputsForInterp_"+config.get('nominal','strnom'),'.root','100','300','75','8E-10')
    outf = ROOT.TFile(outfname,"recreate")
    signalsamps = signal_run2.sigsbyname
    hs = []
    for signalname in signalsamps.keys():
        signame = "holder"
        signame = signalname.replace("-","")
        #make the hist combining all three years scaled to 1 fb @ full Run 2 lumi
        hs.append(signal_run2.getAddedHist(signalname,1.0,empty.Clone(),"sr","h_zp_jigm"))
        hs[-1].SetName(signalname+"_lumixs")#lumi and xs scaled dist
        outf.cd()
        hs[-1].Write()

    outf.Close()
                  

    systdirec = ["up","dwn"]
    for syst in systs:
        #make a file for the  nominal, then each up and down
        if syst == 'nominal':
            continue
        for direc in systdirec:
            signal_systfluc = go.signal_run2(config.get(syst,'pathsig'+direc),'100.0','300.0','75.0','0.8','1.0',[16,17,18],config.get(syst,'str'+direc))
            
            #make the outfile
            outfname = go.makeOutFile('Run2_161718_ZllHbbMET','inputsForInterp_'+config.get(syst,'str'+direc),'.root','100','300','75','8E-10')
            outf = ROOT.TFile(outfname,"recreate")
            signalsamps = signal_run2.sigsbyname
            hs = []
            for signalname in signalsamps.keys():
                signame = "holder"
                signame = signalname.replace("-","")
                #make the hist combining all three years scaled to 1 fb @ full Run 2 lumi
                hs.append(signal_run2.getAddedHist(signalname,1.0,empty.Clone(),"sr","h_zp_jigm"))
                hs[-1].SetName(signalname+"_lumixs")#lumi and xs scaled dist
                outf.cd()
                hs[-1].Write()
                
            outf.Close()

