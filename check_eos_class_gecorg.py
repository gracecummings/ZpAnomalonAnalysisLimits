import gecorg_test as go
import configparser
import ROOT 
import pandas

#Get systematic info                                                                       
config = configparser.RawConfigParser()
config.optionxform = str
fp = open('systematics.ini')
config.read_file(fp)
systs = config.sections()

full_check = True

if full_check:
    for syst in systs:
        if "jec" in syst:
            bkgs = go.backgrounds(config.get(syst,'pathup'),100.0,300.0,75.0,0.8,"alphat_"+config.get(syst,'strup'))
            #bkgs = go.backgrounds(config.get(syst,'pathup'),100.0,300.0,75.0,0.8,config.get(syst,'strup'))
            print("Have to manually change to check the opposite direction of systematic")
            #continue
        elif "jer" in syst:
            bkgs = go.backgrounds(config.get(syst,'pathup'),100.0,300.0,75.0,0.8,"alphat_"+config.get(syst,'strup'))
            #bkgs = go.backgrounds(config.get(syst,'pathup'),100.0,300.0,75.0,0.8,config.get(syst,'strup'))
            print("Have to manually change to check the opposite direction of systematic")
            #continue
        elif "unclmet" in syst:

            bkgs = go.backgrounds(config.get(syst,'pathup'),100.0,300.0,75.0,0.8,"alphat_"+config.get(syst,'strup'))
            #bkgs = go.backgrounds(config.get(syst,'pathup'),100.0,300.0,75.0,0.8,config.get(syst,'strup'))
            print("Have to manually change to check the opposite direction of systematic")
            #continue
        elif "nominal" in syst:
            bkgs = go.backgrounds(config.get('nominal','pathnom'),100.0,300.0,75.0,0.8,"alphat_"+config.get('nominal','strnom'))
            #bkgs = go.backgrounds(config.get('nominal','pathnom'),100.0,300.0,75.0,0.8,config.get('nominal','strnom'))
        else:
            #bkgs = go.backgrounds(config.get(syst,'pathup'),100.0,300.0,75.0,0.8,"alphat_"+config.get(syst,'strup'))
            #bkgs = go.backgrounds(config.get(syst,'pathup'),100.0,300.0,75.0,0.8,config.get(syst,'strup'))
            #print("Have to manually change to check the opposite direction of systematic")
            continue

        print('Checking files for ',syst)        
        for sample in bkgs.bkgs.keys():
            print('    Checking ',sample)
            for year in bkgs.bkgs[sample].keys():
                print('      Checking ',year)
                if "prefire" in syst and year == 18:
                    continue
                for region in bkgs.bkgs[sample][year]:
                    print('        Checking ',region)
                    for f in bkgs.bkgs[sample][year][region][0]:
                        g = ROOT.TFile.Open(f)
                        g.Close()
                    for q in bkgs.bkgs[sample][year][region][1]:
                        p = pandas.read_pickle(q)

else:
    print("Not checking if all files exist/can be accessed. Toggle bool if that is the goal")

check_stack = False

if check_stack:
    print("Checking if the stacking works")
    bkgs = go.backgrounds(config.get('nominal','pathnom'),100.0,300.0,75.0,0.8,config.get('nominal','strnom'))
    
    #Make cleared histograms
    testfile = bkgs.bkgs["DYJetsToLL"][16]['sb'][0][0]#stacked plots should always have DY
    testtfile = ROOT.TFile.Open(testfile)
    h = testtfile.Get('h_z_pt')
    h.SetDirectory(0)
    testtfile.Close()
    
    empty = h.Clone()
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty1 = empty.Clone()
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    empty4 = empty.Clone()
    empty5 = empty.Clone()
    
    
    
    reg = 'sb'
    hname = 'h_z_pt'
    years = [16,17,18]
    
    #Gather histograms
    hdy  = bkgs.getAddedHist(empty2,"DYJetsToLL",reg,hname,years = years)
    print('DY hist integral: ',hdy.Integral())
    htt  = bkgs.getAddedHist(empty3,"TT",reg,hname,years = years)
    print('TT hist integral: ',htt.Integral())
    hzz  = bkgs.getAddedHist(empty4,"ZZTo2L2Q",reg,hname,years = years)
    print('ZZ hist integral: ',hzz.Integral())
    hwz  = bkgs.getAddedHist(empty5,"WZTo2L2Q",reg,hname,years = years)
    print('WZ hist integral: ',hwz.Integral())

else:
    print("Not checking if stacking works in class. Toggle bool if that is the goal")
