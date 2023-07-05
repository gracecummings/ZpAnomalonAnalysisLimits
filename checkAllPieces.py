import argparse
import glob
import ROOT

samp_mumu =['Run2016B-17Jul2018_ver2-v1.SingleMuon',
            'Run2016C-17Jul2018-v1.SingleMuon',
            'Run2016D-17Jul2018-v1.SingleMuon',
            'Run2016E-17Jul2018-v1.SingleMuon',
            'Run2016F-17Jul2018-v1.SingleMuon',
            'Run2016G-17Jul2018-v1.SingleMuon',
            'Run2016H-17Jul2018-v1.SingleMuon',
            'Summer16v3.DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'Summer16v3.DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'Summer16v3.DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'Summer16v3.DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'Summer16v3.DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'Summer16v3.DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'Summer16v3.DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8',
            'Summer16v3.TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',
            'Summer16v3.TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8',
            'Summer16v3.TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8',
            'Summer16v3.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',#ok
            'Summer16v3.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',#ok
            'Run2017B-31Mar2018-v1.SingleMuon',
            'Run2017C-31Mar2018-v1.SingleMuon',
            'Run2017D-31Mar2018-v1.SingleMuon',
            'Run2017E-31Mar2018-v1.SingleMuon',
            'Run2017F-31Mar2018-v1.SingleMuon',
            'Fall17.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_new_pmx',
            'Fall17.TTToHadronic_TuneCP5_13TeV-powheg-pythia8_new_pmx',##,ttbar bg
            'Fall17.TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_new_pmx',##ttbar bg
            'Fall17.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',
            'Fall17.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',
            'Fall17.DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8',
            'Fall17.DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8',
            'Fall17.DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8',
            'Fall17.DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8',
            'Fall17.DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8',
            'Fall17.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8',
            'Fall17.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
            'Run2018C-17Sep2018-v1.SingleMuon',
            'Run2018B-17Sep2018-v1.SingleMuon',
            'Run2018A-17Sep2018-v1.SingleMuon',
            'Run2018D-22Jan2019-v2.SingleMuon',
            'Autumn18.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',##ttbar bg
            'Autumn18.TTToHadronic_TuneCP5_13TeV-powheg-pythia8',##ttbar bg
            'Autumn18.TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',##ttbar bg
            'Autumn18.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',
            'Autumn18.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8',
            'Autumn18.DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
            'Autumn18.DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
            'Autumn18.DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
            'Autumn18.DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
            'Autumn18.DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
            'Autumn18.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
            'Autumn18.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1600-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND1600-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8',
            'Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp1200-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp1300-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp1300-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND800-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND1000-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND800-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1000-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1400-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1200-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1600-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND800-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1000-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1400-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND1600-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8',
#            'Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8']
]
samp_ee = [ 'Run2016B-17Jul2018_ver2-v1.SingleElectron',
            'Run2016C-17Jul2018-v1.SingleElectron',
            'Run2016D-17Jul2018-v1.SingleElectron',
            'Run2016E-17Jul2018-v1.SingleElectron',
            'Run2016F-17Jul2018-v1.SingleElectron',###This is bum!!!
            'Run2016G-17Jul2018-v1.SingleElectron',
            'Run2016H-17Jul2018-v1.SingleElectron',
]
samp_emu = ['Run2016B-17Jul2018_ver2-v1.SingleElectron',
            'Run2016C-17Jul2018-v1.SingleElectron',
            'Run2016D-17Jul2018-v1.SingleElectron',
            'Run2016E-17Jul2018-v1.SingleElectron',
            'Run2016F-17Jul2018-v1.SingleElectron',###This is bum!!!
            'Run2016G-17Jul2018-v1.SingleElectron',
            'Run2016H-17Jul2018-v1.SingleElectron',
            'Run2016B-17Jul2018_ver2-v1.SingleMuon',
            'Run2016C-17Jul2018-v1.SingleMuon',
            'Run2016D-17Jul2018-v1.SingleMuon',
            'Run2016E-17Jul2018-v1.SingleMuon',
            'Run2016F-17Jul2018-v1.SingleMuon',
            'Run2016G-17Jul2018-v1.SingleMuon',
            'Run2016H-17Jul2018-v1.SingleMuon',
            'Summer16v3.TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8',
            'Summer16v3.TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8',
            'Summer16v3.TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8',
            'Run2017B-31Mar2018-v1.SingleElectron',
            'Run2017C-31Mar2018-v1.SingleElectron',
            'Run2017D-31Mar2018-v1.SingleElectron',
            'Run2017E-31Mar2018-v1.SingleElectron',
            'Run2017F-31Mar2018-v1.SingleElectron',
            'Fall17.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_new_pmx',
            'Fall17.TTToHadronic_TuneCP5_13TeV-powheg-pythia8_new_pmx',##,ttbar bg
            'Fall17.TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_new_pmx',##ttbar bg
            'Run2018C-17Sep2018-v1.EGamma',
            'Run2018B-17Sep2018-v1.EGamma',
            'Run2018A-17Sep2018-v1.EGamma',
            'Run2018D-22Jan2019-v2.EGamma',
            'Autumn18.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',##ttbar bg
            'Autumn18.TTToHadronic_TuneCP5_13TeV-powheg-pythia8',##ttbar bg
            'Autumn18.TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',##ttbar bg
    ]

samp_emu_mix = ['Run2016B',
                'Run2016C',
                'Run2016D',
                'Run2016E',
                'Run2016F',
                'Run2016G',
                'Run2016H',
                'Run2017B',
                'Run2017C',
                'Run2017D',
                'Run2017E',
                'Run2017F',
                'Run2018C',
                'Run2018B',
                'Run2018A',
                'Run2018D',
                
]


def nameFormatter(name):
    new = "newname"
    if "Run201" not in name:
        new = name.split("pythia8")[0]+"pythia8"
        if "Fall17.TTTo" in new:
            new = new+"_new_pmx"
    else:
        new = name.split("_")[0]

    return new

def generalChecks(fs,chan):
    samps = []
    if "mumu" in chan:
        samps = samp_mumu
    if "emu" in chan:
        samps = samp_emu
    if "ee" in chan:
        samps = samp_ee

    sampset = set(samps)
    fset    = set(fs)

    missingsamps = sampset-fset

    return missingsamps

def comparableSample(chan):
    samps = []
    if "mumu" in chan:
        samps = samp_mumu
    if "emu" in chan:
        samps = samp_emu
    if "ee" in chan:
        samps = samp_ee

    return samps

regions = ["signalr","sideband","totalr"]
filetypes = ["topiary","upout"]

parser = argparse.ArgumentParser()
if __name__=='__main__':
    parser.add_argument("-dir","--directory", type=str,help = "date folder with output")
    parser.add_argument("-c","--chan", type=str,help = "channel: emu, mumu, ee")
    args = parser.parse_args()

    fdir = glob.glob(args.directory+"/*")

    rootfs = [f.split(args.directory+"/")[-1] for f in fdir if ".root" in f]
    rootnames = [nameFormatter(f) for f in rootfs]
    #print(rootnames)

    chan = args.chan
    compsamps = comparableSample(chan)

    print("Checking the parts of the analysis needed for the {0} channel, generally speaking.".format(chan))
    mostmiss = generalChecks(rootnames,chan)#If there is no version of anything, this will find it
    if len(mostmiss) > 0:
        print("The following samples have no version of processing:")
        missing = list(mostmiss)
        missing = sorted(missing)
        for miss in missing:
            print("         ",miss)

    print("Beyond the above, now checking for subsquent files that are missed.")
    sampset = set(compsamps)
    for ftype in filetypes:
        types = [f for f in rootfs if ftype in f]
        typenames = [nameFormatter(f) for f in types]
        typeset = set(typenames)

        misst = sampset - typeset
        newmiss = misst - (mostmiss & misst)#missing types minus the intersection of the already missing
        missingtypes = list(newmiss)
        missingtypes = sorted(missingtypes)
        if len(missingtypes) > 0:
            print("Found {0} that are missing, and not covered in above list (i.e, questionable version(ing) exists)".format(ftype))
            print("These new missing {0} are:".format(ftype))
            for m in missingtypes:
                print("         ",m)



    print("Checking the validity of the root files. Just because it is there does not mean it is OK.")
    badfs = []
    for f in rootfs:
        tf = ROOT.TFile(args.directory+"/"+f,"read")
        keys = tf.GetListOfKeys()
        if (len(keys) == 0):
            print("Error, tfile has not keys".format(f))
            badfs.append(f)
    print("If no errors were printed, all root files are OK.")
    print("If errors were reported, these are the bad files: ")
    for f in badfs:
        print("         ",f)
    

    #print("Now checking that all regions are there.")
    #sampset = set(compsamps)
    #for ftype in filetypes:
    #    types = [f for f in rootfs if ftype in f]
    #    typenames = [nameFormatter(f) for f in types]
    #    typeset = set(typenames)

    #    misst = sampset - typeset
    #    newmiss = misst - (mostmiss & misst)#missing types minus the intersection of the already missing
    #    missingtypes = list(newmiss)
    #    missingtypes = sorted(missingtypes)
    #    if len(missingtypes) > 0:
    #        print("Found {0} that are missing, and not covered in above list (i.e, questionable version(ing) exists)".format(ftype))
    #        print("These new missing {0} are:".format(ftype))
    #        for m in missingtypes:
    #            print("         ",m)#
