import subprocess
import argparse

parser = argparse.ArgumentParser()

configdict = {
    "none":[],
    "unblind":["-unblind","True"],#override unblinding protections
    "alpham-basic":["-a","True"],
}

cmdlinesystdict = {
    "pdf_up":["-pdfup","True"],
    "pdf_dwn":["-pdfdwn","True"],
    "qcd_up":["-qcdup","True"],
    "qcd_dwn":["-qcddwn","True"],
    "prefire_up":["-prefireup","True"],
    "prefire_dwn":["-prefiredwn","True"],
    "pumapWeight_up":["-puupder","True"],#up weight taken from mapped PU vertice num
    "pumapWeight_dwn":["-pudwnder","True"],#dwn weight taken from mapped PU vertice num
    "pumapPUnum_up":["-puupmap","True"],#pu weight taken when the maped PU number is one std up
    "pumapPUnum_dwn":["-pudwnmap","True"],#pu weight taken when the maped PU number is one std down
    "btag_up":["-syst","btagup"],
    "btag_dwn":["-syst","btagdwn"],
    "muid_up":["-syst","muidup"],
    "muid_dwn":["-syst","muiddwn"],
    "mutrg_dwn":["-syst","mutrigdwn"],
    "mutrg_up":["-syst","mutrigup"],
}

if __name__=='__main__':
    parser.add_argument("-s","--samples",type=str,help="samples to do selections on")
    parser.add_argument("-c","--channel",type=str,help="channel")
    parser.add_argument("-conf","--config",type=str,help="common configuration")
    parser.add_argument("-evntwsyst","--evntwsyst",type=str,help="event weight style systematic")
    
    args = parser.parse_args()

    #Prepare the input files to run over
    inputsraw = args.samples.strip("[]").split(",")
    inputs = list(map(lambda x:x.strip("''"),inputsraw))
    
    #Input cuts: Zpt, Hpt, MET, Tagger, Tagger wp
    cut= ['100.0','300.0','75.0','DeepMassDecorrelTagZHbbvsQCD','0.8']

    #Check addiitonal configuratiosn
    extrun = []
    if args.config:
        if args.config in configdict.keys():
            extrun = configdict[args.config]
        else:
            print("Your extra configuration has not been added, output may not be what you want")

    #Check for event weight systematic flags
    ewts = []
    wsysts = [x.split("_")[0] for x in cmdlinesystdict.keys()]
    if args.evntwsyst:
        ewtsraw = args.evntwsyst.strip("[]").split(",")
        ewts = list(map(lambda x:x.strip("''"),ewtsraw))

    #Do the thing
    print("Doing ZpT cut {0}, HpT cut {1}, MET cut {2}, btag wp {3}".format(cut[0],cut[1],cut[2],cut[4]))

    runstrings = []
    for f in inputs:
        #print("only doing one - make sure you remove this for real running")


        runstring = ["python","doSelections.py","-f",f,"-zpt",cut[0],"-hpt",cut[1],"-met",cut[2],"-sdm","30.0","-b",cut[3],"-wp",cut[4],"-c","mumu","-pumap","True"]#sideband -- default

        if args.evntwsyst and ("jec" not in f) and ("uncl" not in f) and ("jer" not in f):
            #Check if doing an object level syst --- will not combine with event weight syst
            #makes sure to only run event weight syst on nominal topiaries

            for evntw in ewts:
                #Make the event weight systematic stuff
                upflag = []
                dwnflag = []
                
                if evntw in wsysts:
                    upflag = cmdlinesystdict[evntw+"_up"]
                    dwnflag = cmdlinesystdict[evntw+"_dwn"]
                else:
                    print("Your systematic was not found to match an already implemented systematic -- output may not be what you want.")

                runstrings.append(runstring+extrun+upflag)
                runstrings.append(runstring+["-sr","True"]+extrun+upflag)#signalregion
                runstrings.append(runstring+["-tot","True"]+extrun+upflag)#totalregion
                runstrings.append(runstring+extrun+dwnflag)#sideband
                runstrings.append(runstring+["-sr","True"]+extrun+dwnflag)#signalregion
                runstrings.append(runstring+["-tot","True"]+extrun+dwnflag)#totalregion
        else:
            runstrings.append(runstring+extrun)
            runstrings.append(runstring+["-sr","True"]+extrun)#signalregion
            runstrings.append(runstring+["-tot","True"]+extrun)#totalregion
        
        #if args.evntwsyst and (("jec" in  f) or ("uncl" in f) or ("jer" in f)):
        #    print("Do not have a nominal topiary to use for the event weight systematic: ",args.evntwsyst)

    if len(runstrings) == 0:
        print("No commands to run -- check combination of topiaries and systematic flags! Possibly combined object level systematic topiary with event weight shifts")
    for cmd in runstrings:
        print(cmd)
        #print("running commented out for debug")
        subprocess.run(cmd)
