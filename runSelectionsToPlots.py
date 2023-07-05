import subprocess
from datetime import date

if __name__=='__main__':

    #steps to run
    #assumes you have run the whole thing at the start of the day
    #steps = {"selections":True,"uncs":True,"ratios":True,"opts":True}
    steps = {"topiary":True,"selections":False,"uncs":False,"ratios":False,"opts":False,"cutflow":False,"alphaNorm":False,"alphaR":False,"datacards":False,"limitplots":False}
    
    #cut list, Zpt, Hpt, met,btagger,btagwp
    cutlist = [#['150.0','300.0','50.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['150.0','300.0','100.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['100.0','300.0','50.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['100.0','300.0','100.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['200.0','300.0','50.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               ['100.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['100.0','300.0','0.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],#lose alphar
               #['100.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.0'],
               #['125.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['150.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['175.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['200.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['125.0','300.0','50.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['125.0','300.0','100.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['175.0','300.0','50.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['175.0','300.0','100.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['200.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8'],
               #['100.0','-1.0','-1.0','DeepMassDecorrelTagZHbbvsQCD','-1.0'],#ttbar background
               ]

    lumi = "41.53"

    #year specifics, year as 2 digit end, integrated luminosity
    eras = [["16","35.9"],
            ["17","41.53"],
            ["18","59.74"]
            ]

    plots = ['h_z_pt']#['h_h_pt','h_z_pt','h_met','h_nd_jigm','h_zp_jigm','h_h_sd','h_btag']
    
    #topiary sample list: dateforfolder, samplename
    samplelist = [#['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2016B-17Jul2018_ver2-v1.SingleMuon'],
                  #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2016C-17Jul2018-v1.SingleMuon'],
                  #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2016D-17Jul2018-v1.SingleMuon'],
                  #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2016E-17Jul2018-v1.SingleMuon'],
                  #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2016F-17Jul2018-v1.SingleMuon'],
                  #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2016G-17Jul2018-v1.SingleMuon'],
                  #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2016H-17Jul2018-v1.SingleMuon'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016B-17Jul2018_ver2-v1.SingleElectron'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016C-17Jul2018-v1.SingleElectron'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016D-17Jul2018-v1.SingleElectron'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016E-17Jul2018-v1.SingleElectron'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016F-17Jul2018-v1.SingleElectron'],###This is bum!!!
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016G-17Jul2018-v1.SingleElectron'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016H-17Jul2018-v1.SingleElectron'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016B'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016C'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016D'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016E'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016F'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016G'],
                  #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2016H'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8'],
                 #['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],#ok
                 #['analysis_output_ZpAnomalon/2023-05-30','Summer16v3.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],#ok
                 #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2017B-31Mar2018-v1.SingleMuon'],
                 #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2017C-31Mar2018-v1.SingleMuon'],
                 #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2017D-31Mar2018-v1.SingleMuon'],
                 #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2017E-31Mar2018-v1.SingleMuon'],
                 #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2017F-31Mar2018-v1.SingleMuon'],
                 #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2017B-31Mar2018-v1.SingleElectron'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2017C-31Mar2018-v1.SingleElectron'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2017D-31Mar2018-v1.SingleElectron'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2017E-31Mar2018-v1.SingleElectron'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2017F-31Mar2018-v1.SingleElectron'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2017B'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2017C'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2017D'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2017E'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2017F'],
                  #['2021-06-17','ZpAnomalonHZ_UFO-Zp1200-ND175-NS1'],
                  #['2021-06-17','ZpAnomalonHZ_UFO-Zp2000-ND300-NS1'],
                  #['2021-06-17','ZpAnomalonHZ_UFO-Zp2000-ND500-NS200'],
                  #['2021-06-17','ZpAnomalonHZ_UFO-Zp2000-ND800-NS200'],
                  #['2021-06-17','ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Fall17.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_new_pmx'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Fall17.TTToHadronic_TuneCP5_13TeV-powheg-pythia8_new_pmx'],##,ttbar bg
                 ['analysis_output_ZpAnomalon/2023-05-30','Fall17.TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_new_pmx'],##ttbar bg
                 #['analysis_output_ZpAnomalon/2023-05-30','Fall17.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                 ##['analysis_output_ZpAnomalon/2023-05-30','Fall17.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Fall17.DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Fall17.DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Fall17.DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Fall17.DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Fall17.DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Fall17.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Fall17.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8'],
                 #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2018C-17Sep2018-v1.SingleMuon'],
                 #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2018B-17Sep2018-v1.SingleMuon'],
                 #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2018A-17Sep2018-v1.SingleMuon'],
                 #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Run2018D-22Jan2019-v2.SingleMuon'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2018C-17Sep2018-v1.EGamma'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2018B-17Sep2018-v1.EGamma'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2018A-17Sep2018-v1.EGamma'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2018D-22Jan2019-v2.EGamma'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2018C'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2018B'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2018A'],
                 #['emu_2022-07-17_ProperSF_EE_METXY_HEMveto_Pref','Run2018D'],
                 ['analysis_output_ZpAnomalon/2023-05-30','Autumn18.TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8'],##ttbar bg
                 ['analysis_output_ZpAnomalon/2023-05-30','Autumn18.TTToHadronic_TuneCP5_13TeV-powheg-pythia8'],##ttb#ar bg
                 ['analysis_output_ZpAnomalon/2023-05-30','Autumn18.TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'],##ttbar bg
                ####['mumu_2022-03-31_ProperREOIDSF','Autumn18.WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8'],
                #['analysis_output_ZpAnomalon/2023-05-30','Autumn18.WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],####################################2018_new
                #['analysis_output_ZpAnomalon/2023-05-30','Autumn18.ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'],####################################2018_new
                ['analysis_output_ZpAnomalon/2023-05-30','Autumn18.DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                ['analysis_output_ZpAnomalon/2023-05-30','Autumn18.DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                ['analysis_output_ZpAnomalon/2023-05-30','Autumn18.DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                ['analysis_output_ZpAnomalon/2023-05-30','Autumn18.DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                ['analysis_output_ZpAnomalon/2023-05-30','Autumn18.DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                ['analysis_output_ZpAnomalon/2023-05-30','Autumn18.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
               ['analysis_output_ZpAnomalon/2023-05-30','Autumn18.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND1600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
                #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Autumn18.ZpAnomalonHZ_UFO-Zp5500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp1200-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp1300-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp1300-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp1500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp2000-ND800-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND1000-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp2500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND1200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp3000-ND800-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1000-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND1400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp3500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND1600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp4000-ND800-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1000-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND1400-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp4500-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND1600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND200-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['analysis_output_ZpAnomalon/2022-10-10','Autumn18.ZpAnomalonHZ_UFO-Zp5000-ND600-NS1_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp2000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp2500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp2500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp3000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp3000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp3500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp3500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp3500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp4000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp4000-ND1600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp4000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp4000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp4500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp4500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp4500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp4500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp5000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp5000-ND1600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp5000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp5000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp5500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp5500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Summer16.ZpAnomalonHZ_UFO-Zp5500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Summer16.ZpAnomalonHZ_UFO-Zp5500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp1500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp2000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp2000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp2500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp2500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp3000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp3000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp3000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp3500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp3500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp3500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp4000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp4000-ND1600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp4000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp4000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp4500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp4500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp4500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp4500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp5000-ND1200-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp5000-ND1600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp5000-ND400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp5000-ND800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp5500-ND1000-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp5500-ND1400-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2023-05-03_ProperSF_EE_METXY_HEMveto_Pref_qcdscale','Fall17.ZpAnomalonHZ_UFO-Zp5500-ND1800-NS200_TuneCP5_13TeV-madgraph-pythia8'],
        #['mumu_2022-07-18_ProperSF_EE_METXY_HEMveto_Pref','Fall17.ZpAnomalonHZ_UFO-Zp5500-ND600-NS200_TuneCP5_13TeV-madgraph-pythia8'],
            
                  ]

    if steps["topiary"]:
        for samp in samplelist:
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-y","2018"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-y","2017"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu","-upjec","1"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu","-dwnjec","1"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu","-upjer","1"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu","-dwnjer","1"])
            subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu","-upuncl","1"])
            subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu","-dwnuncl","1"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","mumu"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","ee"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","emu"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","emu","-upjec","1"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","emu","-dwnjec","1"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","emu","-upuncl","1"])
            #subprocess.run(["python","runTopiary.py","-s",samp[1],"-c","emu","-dwnuncl","1"])

            
    for cut in cutlist:
        print("Doing ZpT cut {0}, HpT cut {1}, MET cut {2}, btag wp {3}".format(cut[0],cut[1],cut[2],cut[4]))

        #do selections
        if steps["selections"]:
            for samp in samplelist:
                #unblind the signal region
                subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu","-unblind","True"])
                
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True"])#alpha sideband, no syst

                #pileup checking
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True","-puup","True"])#alpha sideband
                ##subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True","-pumap","True"])#alpha sideband
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-sr","True","-puup","True"])#signal region
                ##subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-sr","True","-pumap","True"])#signal region
                ##subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True","-pumap","True","-tot","True"])#alpha totalr
                
                #alpha region validation
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-b",cut[3],"-wp",cut[4],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-d",samp[0],"-v","True","-a","True","-sr","True","-c","mumu"])
                
                if "Run" in samp[1]:
                    continue
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu","-a","True"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu"])

                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu","-a","True"])#alpha total r, no weight syst
                
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True","-syst","btagup"])#alpha sb up
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True","-syst","btagdwn"])#alpha sb dwn

                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu","-syst","btagdwn","-a","True"])#total down for alpah
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu","-syst","btagup","-a","True"])#total up for alpah
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu","-syst","btagdwn"])#sr dwn
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu","-syst","btagup"])#sr up


                #emu
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","emu","-syst","btagup"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","emu","-syst","btagup"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","emu","-syst","btagdwn"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","emu","-syst","btagdwn"])

                #prefire uncs
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","--prefireup","True","-a","True"])#alpha r sb
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","--prefiredwn","True","-a","True"])#alpha r sb
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu","--prefireup","True","-a","True"])#alphar tot
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu","--prefireup","True"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu","--prefiredwn","True","-a","True"])#alpha r tot
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu","--prefiredwn","True"])
                
                #Alpha ratio norm specific
                #if "Run" in samp[1]:#her for sf uncs
                #    continue
                #if "ZpA" not in samp[1]:#If not signal, do other regions
                #    subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True","-qcdup","True"])#alpha sideband
                #    subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True","-qcddwn","True"])#alpha sideband
                #    subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True","-tot","True","-qcdup","True"])#alpha total, for norm
                #    subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True","-tot","True","-qcddwn","True"])#alpha total, for norm
                #if "TT" not in samp[1]:#DY need SR for Alpha, signal needs it for signal, VV needs it

                ###subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-sr","True","-qcdup","True"])#alpha signal region, will have different tag. Only needed for DY
               ### subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-sr","True"])
                ##subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-sr","True","-qcddwn","True"])#alpha signal region, will have different tag. Only needed for DY


                #Alpha Ratio
                #sideband with no met cut
                
                #if "Run" in samp[1]:
                #    continue
                #total region witn no met cut
                
                
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu","-a","True","-syst","muidup"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu","-a","True","-syst","muidup"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu","-a","True","-syst","muiddwn"])

                #Alpha method rvalidation region stuff
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu","-a","True","-v","True"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-a","True","-v","True"])

                #####ttbar background
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","emu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","emu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","emu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu","--prefireup","True"])
                #ubprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu","--prefiredwn","True"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","emu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","0.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","emu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","0.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","emu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","0.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","emu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","0.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-c","mumu"])

                ####sf uncs
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-syst","btagdwn","-c","mumu"])
                #####subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-syst","btagdwn","-c","mumu"])
                ##subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-syst","btagdwn","-c","mumu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-syst","btagup","-c","mumu"])
                #####subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-syst","btagup","-c","mumu"])
                ##subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-syst","btagup","-c","mumu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-syst","muidup","-c","mumu"])
                #####subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-syst","muidup","-c","mumu"])
                ##subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-syst","muidup","-c","mumu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","mumu","-syst","muiddwn"])
                #####subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-syst","muiddwn","-c","mumu"])
                ##subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-tot","True","-syst","muiddwn","-c","mumu"])

                ####subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-syst","mutrigup","-c","mumu"])
                ####subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-syst","mutrigdwn","-c","mumu"])
                #####subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu","-pdfup","True"])
                #####subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-c","mumu","-pdfdwn","True"])
                

                ###Electron or emu channel only
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","emu","-syst","eliddwn"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-syst","eliddwn","-c","emu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","emu","-syst","elidup"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-syst","elidup","-c","emu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","emu","-syst","elrecodwn"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-syst","elrecodwn","-c","emu"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-c","emu","-syst","elrecoup"])
                #subprocess.run(["python","doSelections.py","-f",samp[1],"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-d",samp[0],"-sr","True","-syst","elrecoup","-c","emu"])


                
                
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
