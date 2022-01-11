import subprocess
from datetime import date

if __name__=='__main__':

    #steps to run
    #assumes you have run the whole thing at the start of the day
    #steps = {"selections":True,"uncs":True,"ratios":True,"opts":True}
    steps = {"topiary":False,"selections":False,"uncs":False,"ratios":False,"opts":False,"cutflow":False,"alphaNorm":False,"alphaR":False,"datacards":False,"limitplots":False}
    
    #cut list, Zpt, Hpt, met,btagger,btagwp
    cutlist = [#['150.0','300.0','50.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['150.0','300.0','100.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['100.0','300.0','50.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['100.0','300.0','100.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['200.0','300.0','50.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               ######['100.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['125.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['150.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['175.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['200.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['125.0','300.0','50.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['125.0','300.0','100.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['175.0','300.0','50.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['175.0','300.0','100.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['200.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               ['100.0','0.0','0.0','DeepMassDecorrelTagZHbbvsQCD','0.0'],#ttbar background
               ]

    lumi = "41.53"

    #year specifics, year as 2 digit end, integrated luminosity
    eras = [["16","35.9"],
            ["17","41.53"],
            ["18","59.74"]
            ]

    plots = ['h_z_pt']#['h_h_pt','h_z_pt','h_met','h_nd_jigm','h_zp_jigm','h_h_sd','h_btag']

    #topiary sample list: dateforfolder, samplename
    samplelist = [['','Run2016B-17Jul2018_ver2-v1.SingleMuon'],
                  [''.'Run2016C-17Jul2018-v1.SingleMuon'],
                  [''.'Run2016D-17Jul2018-v1.SingleMuon'],
                  [''.'Run2016E-17Jul2018-v1.SingleMuon'],
                  [''.'Run2016F-17Jul2018-v1.SingleMuon'],
                  [''.'Run2016G-17Jul2018-v1.SingleMuon'],
                  ['','Run2016H-17Jul2018-v1.SingleMuon'],
                  ['','Summer16v3.DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                  ['','Summer16v3.DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                  ['','Summer16v3.DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                  ['','Summer16v3.DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                  ['','Summer16v3.DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                  ['','Summer16v3.DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                  ['','Summer16v3.DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                  ['','Summer16v3.TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8_'],
                  ['','Summer16v3.TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8'],
                  ['','Summer16v3.TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8'],
                  ['','Summer16v3.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  ['','Summer16v3.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  ##['pfMET_nominal','Run2017B-31Mar2018-v1.SingleMuon'],
                  ##['pfMET_nominal','Run2017C-31Mar2018-v1.SingleMuon'],
                  ##['pfMET_nominal','Run2017D-31Mar2018-v1.SingleMuon'],
                  ##['pfMET_nominal','Run2017E-31Mar2018-v1.SingleMuon'],
                  ##['pfMET_nominal','Run2017F-31Mar2018-v1.SingleMuon'],
                  #['2021-06-17','ZpAnomalonHZ_UFO-Zp1200-ND175-NS1'],
                  #['2021-06-17','ZpAnomalonHZ_UFO-Zp2000-ND300-NS1'],
                  #['2021-06-17','ZpAnomalonHZ_UFO-Zp2000-ND500-NS200'],
                  #['2021-06-17','ZpAnomalonHZ_UFO-Zp2000-ND800-NS200'],
                  #['2021-06-17','ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1'],
                  ##['pfMET_nominal','Fall17.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_new_pmx'],
                  ##['pfMET_nominal','Fall17.TTToHadronic_TuneCP5_13TeV-powheg-pythia8_new_pmx'],##ttbar bg
                  ##['analysis_output_ZpAnomalon/2021-12-09','Fall17.TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_new_pmx'],##ttbar bg
                  ##['analysis_output_ZpAnomalon/2021-12-09','Fall17.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  ##['analysis_output_ZpAnomalon/2021-12-09','Fall17.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                  ##['pfMET_nominal','Fall17.DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ##['pfMET_nominal','Fall17.DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ##['pfMET_nominal','Fall17.DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ##['pfMET_nominal','Fall17.DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ##['pfMET_nominal','Fall17.DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ##['pfMET_nominal','Fall17.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ##['pfMET_nominal','Fall17.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8'],
                  ##['pfMET_nominal','Run2018C-17Sep2018-v1.SingleMuon'],
                  ##['pfMET_nominal','Run2018B-17Sep2018-v1.SingleMuon'],
                  ##['pfMET_nominal','Run2018A-17Sep2018-v1.SingleMuon'],
                  ##['pfMET_nominal','Run2018D-22Jan2019-v2.SingleMuon'],#ttbar background
                  ##['pfMET_nominal','Autumn18.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8'],##ttbar bg
                  ##['pfMET_nominal','Autumn18.TTToHadronic_TuneCP5_13TeV-powheg-pythia8'],##ttbar bg
                  ##['pfMET_nominal','Autumn18.TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'],##ttbar bg
                  ##['analysis_output_ZpAnomalon/2021-12-09','Autumn18.DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                  ##['analysis_output_ZpAnomalon/2021-12-09','Autumn18.DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                  ##['analysis_output_ZpAnomalon/2021-12-09','Autumn18.DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                  ##['analysis_output_ZpAnomalon/2021-12-09','Autumn18.DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                  ##['analysis_output_ZpAnomalon/2021-12-09','Autumn18.DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                  ##['analysis_output_ZpAnomalon/2021-12-09','Autumn18.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                  ##['analysis_output_ZpAnomalon/2021-12-09','Autumn18.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp1200-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp1300-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND800-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND1000-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND800-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1000-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND800-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1000-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND1200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND1600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND1600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  #['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  ##['pfMET_nominal','Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                  ]

    if steps["topiary"]:
        for samp in samplelist:
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-y","2018"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-y","2017"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu","-upjec","1"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu","-dwnjec","1"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu","-upuncl","1"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu","-dwnuncl","1"])
            subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu"])
            
    for cut in cutlist:
        print("Doing ZpT cut {0}, HpT cut {1}, MET cut {2}, btag wp {3}".format(cut[0],cut[1],cut[2],cut[4]))

        #do selections
        if steps["selections"]:
            for samp in samplelist:
                #subprocess.run(["python","quickSelectionMinimal.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0]])
                
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0]])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-sr","True"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-c","True"])

                #####ttbar background
                subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","0.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-c","True","-channel","mumu"])

                ####btagging sf
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-syst","btagdwn"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-sr","True","-syst","btagdwn"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-c","True","-syst","btagdwn"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-syst","btagup"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-sr","True","-syst","btagup"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-date",samp[0],"-c","True","-syst","btagup"])
                
        if steps["alphaNorm"]:
            print("    Doing alpha method normalization")
            subprocess.run(["python","doAlphaRatioNormalization_systematics.py","-m", cut[2], "-z", cut[0],"-j",cut[1],"-wp",cut[4],"-dir","analysis_output_ZpAnomalon/2021-11-05","-v","False"])

        if steps["alphaR"]:
            print("    Doing alpha ratio extrapolation")
            subprocess.run(["python","doAlphaRatioFits_systematics.py","-m", cut[2], "-z", cut[0],"-j",cut[1],"-wp",cut[4]])

        for era in eras:
            print("   Beginning plottng and analysis for year {0}, with a luminosity of {1}".format(era[0],era[1]))
            if steps["uncs"]:
                #subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0],"-r","totalr"])
                #subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0],"-r","sideband"])
                subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date","2021-11-05","-y",era[0],"-r","sideband"])
                #subprocess.run(["python","doStackedUncertainty.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0],"-r","signalr"])

        #stack all  
        if steps["ratios"]:
            #subprocess.run(["python","stackAll.py","-L",era[1],"-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0]])
            subprocess.run(["python","stackAll.py","-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-r","sideband"])#add -y for a specific year
            #subprocess.run(["python","stackAll.py","-L","41.53","-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y","17")
            #subprocess.run(["python","stackAll.py","-L",lumi,"-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today())+'/signalregion_only/'])

        #Optimization Plots
        if steps["opts"]:
            for plot in plots:
                subprocess.run(["python","stackForOptimization.py","-L",era[1],"-x","10.0","-p",plot,"-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0]])
                #subprocess.run(["python","stackForOptimization.py","-L",lumi,"-x","100.0","-p",plot,"-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",'2021-02-03'])

        if steps["cutflow"]:
            print("Creating cutflow table")
            subprocess.run(["python","doCutFlow.py","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-r","sideband"])
            #subprocess.run(["python","doCutFlowSig.py","-L",era[1],"-x","10.0","-m",cut[2],"-z",cut[0],"-j",cut[1],"-wp",cut[4],"-date",str(date.today()),"-y",era[0]])

        if steps["datacards"]:
            subprocess.run(["python", "makeDataCardsAndRootFiles.py","-z",cut[0],"-j",cut[1],"-m",cut[2],"-wp",cut[4]])

        if steps["limitplots"]:
            subprocess.run(["python", "makeLimiPlot.py","-z",cut[0],"-j",cut[1],"-m",cut[2],"-wp",cut[4]])
            
