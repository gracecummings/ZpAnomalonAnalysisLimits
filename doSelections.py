import uproot3 as up3
import pandas as pd
import numpy as np
import boost_histogram as bh
import argparse
import glob
import gecorg_test as go
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()

def makeEventWeightSeries(df,syststr):
    colnames = list(df.columns)
    wnames = [x for x in colnames if ("event_weight" in x) and ("unc" not in x)]
    sfnames = [x.split("event_weight_")[1] for x in wnames]

    sfname = "_".join(sfnames)
    
    wdf = df[wnames]
    ewdf = wdf.prod(axis=1)
    return ewdf,sfname

def makeEventWeightSeriesWithUncertainties(df,syststr):
    colnames = list(df.columns)
    wnames = [x for x in colnames if ("event_weight" in x) and ("unc" not in x)]
        
    sfnames = [x.split("event_weight_")[1] for x in wnames]
    sfname = "_".join(sfnames)

    #print(sfnames)

    if "up" in syststr:
        uncname = syststr.split("up")[0]
        unctag = syststr.split("up")[0]+"uncup"
        uncsign = 1.0
    if "dwn" in syststr:
        uncname = syststr.split("dwn")[0]
        unctag = syststr.split("dwn")[0]+"uncdwn"
        uncsign = -1.0

    sfname = sfname.replace(uncname,syststr)

    wcolname = "event_weight_"+uncname
    if wcolname not in colnames:
        print("You scale factor is not in the dataframe, aborting")
        exit()
    if not ("event_weight_"+unctag in colnames):
        print("You scale factor uncertainty is not in the dataframe, aborting")
        exit()
        
    unccol = df["event_weight_"+unctag]        
    wcol   = df[wcolname]
    #print(wcol)
    wnames.remove(wcolname)
    
    wdf = df[wnames].copy()
    wdf[wcolname] = wcol+uncsign*unccol
    ewdf = wdf.prod(axis=1)
    #print(ewdf)
    return ewdf,sfname
    
def boostUnc(values,weights,nbins,binstart,binstop):
    boosth = bh.Histogram(bh.axis.Regular(bins=nbins,start=binstart,stop=binstop),storage=bh.storage.Weight())
    boosth.fill(values,weight=weights)
    boostvar = boosth.view().variance
    boosterr = np.sqrt(boostvar)
    return boosterr

def boosthist(values,weights,nbins,binstart,binstop):
    boosth = bh.Histogram(bh.axis.Regular(bins=nbins,start=binstart,stop=binstop),storage=bh.storage.Weight())
    return boosth


def wrapPhi(phi):
    if phi < 0:
        wphi = -1*phi
    else:
        wphi = phi
    return wphi

def wrapDeltaPhi(dphi):
    if dphi > 3.14159:
        dp = 2*3.14159 - dphi
    else:
        dp = dphi
    return dp
 

def deltaPhi(v1phi,v2phi):
    dp = abs(v1phi-v2phi)
    wdp = dp.map(wrapDeltaPhi)
    #This returns a df with the same number of events as input
    #If a cut is introduced, need to carry weight column
    return wdp

def deltaR(v1phi,v2phi,v1eta,v2eta):
    dR = ((v2phi-v1phi)**2+(v2eta-v1eta)**2)**(1/2)
    return dR

def getContentsOfSample(listoffiles):
    fs = [x.split("/")[-1] for x in listoffiles]
    fronts = [x.split("_Zpt")[0] for x in fs]
    jectypes = [x.split("_")[-1] for x in fronts]
    names = [x.split("_topiary")[0] for x in fronts]

    names = list(set(names))
    jectypes = list(set(jectypes))
    return names,jectypes
    

if __name__=='__main__':
    parser.add_argument("-f","--sample",help = "sample file")
    parser.add_argument("-o","--output",help = "output file name")
    parser.add_argument("-b","--btagger",help = "btagger selection, need name part of key")
    parser.add_argument("-wp","--btagWP",type=float,default=0.6,help = "doulbeB tagger working point (default M1)")
    parser.add_argument("-zpt","--zPtCut",type=float,default = 100.0,help = "pT cut on Z")
    parser.add_argument("-hpt","--hPtCut",type=float,default = 250.0,help = "pT cut on h")
    parser.add_argument("-met","--metPtCut",type=float,default = 50.0,help = "pT cut on met")
    parser.add_argument("-sdm","--sdmCut",type=float,default = 10.0,help = "lowest soft drop mass cut")
    parser.add_argument("-d","--directory",type=str,help = "where are your topiary plots?")
    parser.add_argument("-sr","--signalregion",type=bool,help = "do you want a signal region plot?")
    parser.add_argument("-tot","--comboregion",type=bool,help = "do you want combined SR and SB?")
    parser.add_argument("-v","--validation",type=bool,help = "validation region bounds?")
    parser.add_argument("-syst","--systematics",type=str)
    parser.add_argument("-prefireup","--prefireup",type=bool,help = "up prefire uncs?")
    parser.add_argument("-prefiredwn","--prefiredwn",type=bool,help = "down prefire uncs?")
    parser.add_argument("-pdfup","--pdfup",type=bool,help = "up pdf uncs?")
    parser.add_argument("-pdfdwn","--pdfdwn",type=bool,help = "down pdf uncs?")
    parser.add_argument("-qcdup","--qcdup",type=bool,help = "up qcd uncs?")
    parser.add_argument("-qcddwn","--qcddwn",type=bool,help = "down qcd uncs?")
    parser.add_argument("-pu","--pileup",type=bool,help = "apply pu weights?")
    parser.add_argument("-puup","--pileupup",type=bool,help = "apply up pu weights?")
    parser.add_argument("-pudwn","--pileupdwn",type=bool,help = "apply down pu weights?")
    #special pileup weights
    parser.add_argument("-pumap","--pileupmap",type=bool,help = "apply pu weights from mapping?")
    parser.add_argument("-puupder","--pileupupder",type=bool,help = "apply up pu weights derived from mapped value?")
    parser.add_argument("-pudwnder","--pileupdwnder",type=bool,help = "apply down pu weights derived from mapped value?")
    parser.add_argument("-puupmap","--pileupupmap",type=bool,help = "apply up pu weights derived from up-mapped value?")
    parser.add_argument("-pudwnmap","--pileupdwnmap",type=bool,help = "apply down pu weights derived from up-mapped value?")
    #other
    parser.add_argument("-a","--alphar",type=bool,help = "alpha ratio regions?")
    parser.add_argument("-c","--chan",type=str)
    parser.add_argument("-unblind","--unblind",type=bool,help = "override unblinding protections and unblind")
    args = parser.parse_args()

    sampname   = args.sample
    sdmcut = args.sdmCut
    zptcut = args.zPtCut
    hptcut = args.hPtCut
    metcut = args.metPtCut
    btaggr = args.btagger
    btagwp = args.btagWP
    sr     = args.signalregion
    comb   = args.comboregion
    valid  = args.validation
    systl  = args.systematics
    alphatest = args.alphar
    channel = args.chan
    topdir = args.directory
    print(topdir)

    print(sampname)
    print(channel)
    inputfiles = [sampname]
    #inputfiles = ["root://cmseos.fnal.gov//store/group/lpcboostres/topiaries_systematics-dwnjer_2023-10-18/Autumn18.TTTo2L2Nu_TuneCP5_topiary_mumu_systjerdwn_Zptcut0.0_Hptcut250.0_metcut0.0_btagwp0.0.root"]
    #inputfiles = glob.glob(topdir+'/'+sampname+'*_topiary_'+channel+'*.root')
    #inputfiles = glob.glob('analysis_output_ZpAnomalon/2022-03-14/noLeadingReq/Run2016H-17Jul2018-v1.SingleElectron_topiary_emu_systnominal_elecTrigNoLeadingRequirement_Zptcut0.0_Hptcut250.0_metcut0.0_btagwp0.0*')
    #inputfiles = glob.glob('analysis_output_ZpAnomalon/2022-03-28/Autumn18.TTToSemiLeptonic_TuneCP5_13TV-powheg-pythia8_topiary_mumu_systnominal_Zptcut0.0_Hptcut250.0_metcut0.0_btagwp0.0*')
    print(inputfiles)
    innames,jecs = getContentsOfSample(inputfiles)
    samp = innames[0]
    print("The parameters you are invoking: ")
    print("    ",innames,jecs)

    for i,jectype in enumerate(jecs):
        stype,year = go.sampleType(samp)

        if valid:
            print("    Doing validation of alpha method cuts")
            if sr:
                print("    'signalr' labeled events are the stanard signalr in jet mass, but with an inverted met cut")
            else:
                print("    'sideband' labeled events are in soft drop mass bands [30,70],[150,5000], with inverted met cut")

        metstr = ''
        branches = [b'ZCandidate_*',
                    b'RunNum',
                    b'LumiBlockNum',
                    b'EvtNum',
                    b'hCandidate_*',
                    b'metsuable',
                    b'metphiusable',
                    b'ZPrime_mass_est',
                    b'ND_mass_est',
                    b'NS_mass_est',
                    b'event_weight*',
                    b'metxycorr',
                    b'metxyphicorr',
                    b'ZupCandidate_*',
                    b'ZdnCandidate_*',
                    #b'metmuup',
                    #b'metmudn',
        ]

        if (stype > 0):
            mcbranches = [b'ghCandidate_*',
                          b'gzCandidate_*',
                          b'pdfweight_*',
                          b'qcdweight_*',
                          #b'puWeight',
                          #b'puSysUp',
                          #b'puSysDown',
                          b'puweight*',
            ]
            branches.extend(mcbranches)

        if (channel == "emu"):
            mcbranches = [b'LMuCandidate_*',
                          b'sLMuCandidate_*',
                          b'LEleCandidate_*',
                          b'sLEleCandidate_*']
            branches.extend(mcbranches)

        if (channel == "mumu"):
            mcbranches = [b'LMuCandidate_*',
                          b'sLMuCandidate_*',]
            branches.extend(mcbranches)

        if (stype > 0) and ((year == 16) or (year == 17)):
            mcbranches = [b'NonPrefir*']
            branches.extend(mcbranches)


        if (("Run" in sampname) and ("emu" in channel) and not (("Single" in sampname) or ("EGamma" in sampname))):
            print(sampname)
            print("Going to combine datasets")
            muf = glob.glob(topdir+'/'+sampname+'*SingleMuon*_topiary_'+channel+'_*'+jectype+'*.root')
            if (year < 18):
                euf = glob.glob(topdir+'/'+sampname+'*SingleElectron*_topiary_'+channel+'_*'+jectype+'*.root')
            if (year == 18):
                euf = glob.glob(topdir+'/'+sampname+'*EGamma*_topiary_'+channel+'_*'+jectype+'*.root')
            print("Doing selections on:")
            print("    ",muf[0],euf[0])
            print("    ",sampname)
            
            upfmu = up3.open(muf[0])
            #upfmu = open(muf[0])
            dictupmu = dict(upf)
            mtree = dictupmu['PreSelection;1']
            mevents = mtree.pandas.df(branches=branches)
            upfe = up3.open(euf[0])
            #upfe = open(euf[0])
            dictupe = dict(upfe)
            etree = dictupe['PreSelection;1']
            eevents = etree.pandas.df(branches=branches)
            print("Number of events in muon dataset ",len(mevents))
            print("Number of events in elec dataset ",len(eevents))
            
            frames = [mevents,eevents]
            mixdf = pd.concat(frames)
            print("Number of events in straight mixed set ",len(mixdf))
            print("dropping duplicates")
            mixdf = mixdf.drop_duplicates(subset = ['RunNum','LumiBlockNum','EvtNum'])
            print("Number of events after dropping duplicates ",len(mixdf))
            events = mixdf
            goodname = sampname+"combined"

        else:
            print("No weird combos going on")
            #inputfiles = glob.glob(topdir+'/'+samp+'*_topiary_'+channel+'_*'+jectype+'*.root')
            print("Doing selections on:")
            print("    ",inputfiles[:1])
            print("    ",samp)
            #events = up3.pandas.iterate(inputfiles[:1],'PreSelection;1',branches=branches)
            upf = up3.open(inputfiles[0])
            dictup = dict(upf)
            tree = upf['PreSelection;1']
            events = tree.pandas.df(branches=branches)
            goodname = samp

        #print(events)
        print("    Concerning jet systematics:")
        print("    ",jectype)


            
        #events = tree.pandas.df()
        #print(events)
    
        # for b in events:
        #print(type(b))
        #print(events.keys())
        #print(b)
        #print("Doing SD mass lower cut of :" ,sdmcut)
        #do some cuts
        #print("Number of events in chunk ",len(events))

        if alphatest:
            metcut = 0.0



        #print("all prefire weights sum, topiary selecs    ",sum(events["NonPrefiringProb"]))
        #print("all prefireup weights sum, topiary selecs  ",sum(events["NonPrefiringProbUp"]))
        #print("all prefiredwn weights sum, topiary selecs ",sum(events["NonPrefiringProbDown"]))

        #Added for muon investigation
        #leadmudf = events[events['LMuCandidate_pt']-events['LMuCandidate_ptunc'] > 60.0]
        #subleadmudf = leadmudf[leadmudf['sLMuCandidate_pt']-leadmudf['sLMuCandidate_ptunc'] > 20.0]
        #leadmudf = events[events['LMuCandidate_pt']+events['LMuCandidate_ptunc'] > 60.0]
        #subleadmudf = leadmudf[leadmudf['sLMuCandidate_pt']+leadmudf['sLMuCandidate_ptunc'] > 20.0]
            
        #Usual analysis    
        zptdf  = events[events['ZCandidate_pt'] > zptcut]
        #zptdf  = subleadmudf[subleadmudf['ZupCandidate_pt'] > zptcut]
        #zptdf  = subleadmudf[subleadmudf['ZCandidate_pt'] > zptcut]
        #metdf   = zptdf[zptdf['metsuable'] > metcut]
        metdf   = zptdf[zptdf['metxycorr'] > metcut]
        #metdf   = zptdf[zptdf['metmudn'] > metcut]
        #metdf   = zptdf[zptdf['metmuup'] > metcut]
        hptdf  = metdf[metdf['hCandidate_pt'] > hptcut]
        btdf   = hptdf[hptdf['hCandidate_'+btaggr] > float(btagwp)]

        #Actual Analysis
        #if not (valid or alphatest):
        if not valid:
            srup   = btdf[btdf['hCandidate_sd'] > 70.]
            bldf   = srup[srup['hCandidate_sd'] < 150.]#full blinded region
            srdf   = bldf[bldf['hCandidate_sd'] > 110.]#Higgs Peak
            lowsbh = btdf[btdf['hCandidate_sd'] <= 70.]
            lowsb  = lowsbh[lowsbh['hCandidate_sd'] > 30.]
            highsb = btdf[btdf['hCandidate_sd'] >= 150.]
            sbdf   = pd.concat([lowsb,highsb])
            totdf = pd.concat([bldf,lowsb,highsb])
            
        #Validation region for alpha method for DY
        if valid and alphatest:
            #validdf = btdf[btdf['metsuable'] <= args.metPtCut]#metcut is set to zero before here
            validdf = btdf[btdf['metxycorr'] <= args.metPtCut]#metcut is set to zero before here
            srup   = validdf[validdf['hCandidate_sd'] > 70.]
            bldf   = srup[srup['hCandidate_sd'] < 150.]#full blinded region
            srdf   = bldf[bldf['hCandidate_sd'] > 110.]#Higgs Peak
            lowsbh = validdf[validdf['hCandidate_sd'] <= 70.]
            lowsb  = lowsbh[lowsbh['hCandidate_sd'] > 30.]
            highsb = validdf[validdf['hCandidate_sd'] >= 150.]
            sbdf   = pd.concat([lowsb,highsb])
            totdf = pd.concat([bldf,lowsb,highsb])
        
        if alphatest and not valid:
            normdf = totdf
            metcut = args.metPtCut#to return it for the signal region cut

        region = "sideband"
        if not (valid or alphatest):
            if stype > 0:
                if sr:
                    fdf = srdf
                    region = "signalr"
                    print("    using signal region selections")
                elif comb:
                    fdf = totdf
                    region = "totalr"
                    print("    using full region selections")
                else:
                    fdf = sbdf
                    region = "sideband"
                    print("    using sideband selections")
            elif stype <= 0 and ("emu" in channel) and sr:
                fdf = srdf
                region = "signalr"
                print("    using signal region selections")
            elif stype <= 0 and ("emu" in channel) and comb:
                fdf = totdf
                region = "totalr"
                print("    using full region selections")
            elif stype <= 0 and args.unblind:
                fdf = srdf
                region = "signalr"
                print("    using signal region selections")
                print("    YOU HAVE UNBLINDED YOU HAVE UNBLINDED -- YOU HAVE UNBLINDED IN CASE YOU STILL HAVE NOT NOTICED")
            else:
                fdf = sbdf
                print("    using sideband selections")
        if valid:
            if sr:
                fdf = srdf
                region = "validationr"
            elif comb:
                fdf = totdf
                region = "totalr"
            else:
                fdf = sbdf
                region = "validationsideband"
        if alphatest and not valid:
            if (sr and stype > 0):
                #fdf = srdf[srdf['metsuable'] > metcut]#added back in the met cut for the sr
                fdf = srdf[srdf['metxycorr'] > metcut]#added back in the met cut for the sr
                region = "signalr_alphat"
                print("    using signal region selections")
            elif (comb and stype > 0):
                fdf = normdf
                region = "totalr_alphat"
                print("    Using a total soft drop mass region without met cut")
            else:
                fdf = sbdf
                region = "sideband_alphat"
                print("    using a softdrop mass sideband without met cut")


        #print("number of btag passing events ",len(btdf))
        
        #calculated quantities
        deltaphizhdf   = deltaPhi(fdf['ZCandidate_phi'],fdf['hCandidate_phi'])
        #deltaphizmetdf = deltaPhi(fdf['ZCandidate_phi'],fdf['metphiusable'])
        #deltaphihmetdf = deltaPhi(fdf['hCandidate_phi'],fdf['metphiusable'])
        deltaphizmetdf = deltaPhi(fdf['ZCandidate_phi'],fdf['metxyphicorr'])
        deltaphihmetdf = deltaPhi(fdf['hCandidate_phi'],fdf['metxyphicorr'])
        deltaRzhdf     = deltaR(fdf['ZCandidate_phi'],fdf['hCandidate_phi'],fdf['ZCandidate_eta'],fdf['hCandidate_eta'])
        tscalmom       = fdf['ZCandidate_pt']+fdf['hCandidate_pt']+fdf['metxycorr']
        metscalmomfrac = fdf['metxycorr']/tscalmom
        zptscalmomfrac = fdf['ZCandidate_pt']/tscalmom
        hptscalmomfrac = fdf['hCandidate_pt']/tscalmom

        #the leading and subleading lepton stuff
        if (channel == "mumu"): 
            deltaRlmuhdf   = deltaR(fdf['LMuCandidate_phi'],fdf['hCandidate_phi'],fdf['LMuCandidate_eta'],fdf['hCandidate_eta'])
            deltaRslmuhdf   = deltaR(fdf['sLMuCandidate_phi'],fdf['hCandidate_phi'],fdf['sLMuCandidate_eta'],fdf['hCandidate_eta'])
            deltaRslmulmudf = deltaR(fdf['sLMuCandidate_phi'],fdf['LMuCandidate_phi'],fdf['sLMuCandidate_eta'],fdf['LMuCandidate_eta'])

        if (channel == "emu"):
            print("IN the emu chanel")
            #print(fdf['LEleCandidate_pt'])
            #print(fdf['LMuCandidate_pt'])
            print(" Length of final selection df ",len(fdf))

            leadmu = fdf[fdf['LMuCandidate_pt'] > fdf['LEleCandidate_pt']]
            leadel = fdf[fdf['LEleCandidate_pt'] > fdf['LMuCandidate_pt']]

            #print(leadmu)
            #print(leadel)

            print(" Length of leading muon df in final selection ",len(leadmu))
            print(" Length of leading electron df in",len(leadel))

            ###new way
            leadmuc = leadmu.copy()
            leadelc = leadel.copy()

            leadordered = pd.concat([leadmuc,leadelc])

            #gotta make the hCandidate stuff with the same variables
            #not really, could use the leadordered dataframe
            leadinglephcandphi = leadmu['hCandidate_phi'].append(leadel['hCandidate_phi'])
            leadinglephcandeta = leadmu['hCandidate_eta'].append(leadel['hCandidate_eta'])

            leadpt = leadmu['LMuCandidate_pt'].append(leadel['LEleCandidate_pt'])
            leadphi = leadmu['LMuCandidate_phi'].append(leadel['LEleCandidate_phi'])
            leadeta = leadmu['LMuCandidate_eta'].append(leadel['LEleCandidate_eta'])
            sleadpt = leadmu['sLEleCandidate_pt'].append(leadel['sLMuCandidate_pt'])
            sleadphi = leadmu['sLEleCandidate_phi'].append(leadel['sLMuCandidate_phi'])
            sleadeta = leadmu['sLEleCandidate_eta'].append(leadel['sLMuCandidate_eta'])
            allmupt =  leadmu['LMuCandidate_pt'].append(leadel['sLMuCandidate_pt'])
            allelpt = leadmu['sLEleCandidate_pt'].append(leadel['LEleCandidate_pt'])
            allmuphi =  leadmu['LMuCandidate_phi'].append(leadel['sLMuCandidate_phi'])
            allelphi = leadmu['sLEleCandidate_phi'].append(leadel['LEleCandidate_phi'])
            allmueta =  leadmu['LMuCandidate_eta'].append(leadel['sLMuCandidate_eta'])
            alleleta = leadmu['sLEleCandidate_eta'].append(leadel['LEleCandidate_eta'])
            deltaRllephdf   = deltaR(leadphi,leadinglephcandphi,leadeta,leadinglephcandeta)
            #deltaR(leadphi,leadordered['hCandidate_phi'],leadeta,leadordered['hCandidate_eta'])
            deltaRslephdf   = deltaR(sleadphi,leadinglephcandphi,sleadeta,leadinglephcandeta)
            deltaRsleplepdf = deltaR(sleadphi,leadphi,sleadeta,leadeta)
            deltaRemuLeadMuh = deltaR(leadmuc['LMuCandidate_phi'],leadmuc['hCandidate_phi'],leadmuc['LMuCandidate_eta'],leadmuc['hCandidate_eta'])
            deltaRemuSubLeadElh = deltaR(leadmuc['sLEleCandidate_phi'],leadmuc['hCandidate_phi'],leadmuc['sLEleCandidate_eta'],leadmuc['hCandidate_eta'])
            deltaRemuLeadElh = deltaR(leadelc['LEleCandidate_phi'],leadelc['hCandidate_phi'],leadelc['LEleCandidate_eta'],leadelc['hCandidate_eta'])
            deltaRemuSubLeadMuh = deltaR(leadelc['sLMuCandidate_phi'],leadelc['hCandidate_phi'],leadelc['sLMuCandidate_eta'],leadelc['hCandidate_eta'])

            if not systl:
                lepdrweights,notused  = makeEventWeightSeries(leadordered,systl)
                lmuweights,notused2   = makeEventWeightSeries(leadmuc,systl)
                lelweights,notused3   = makeEventWeightSeries(leadelc,systl)
            if systl:
                lepdrweights,notused = makeEventWeightSeriesWithUncertainties(leadordered,systl)
                lmuweights,notused2   = makeEventWeightSeriesWithUncertainties(leadmuc,systl)
                lelweights,notused3   = makeEventWeightSeriesWithUncertainties(leadelc,systl)

        if (stype > 0 and channel != "emu"):#reclustering comments
            deltaRlmughdf  = deltaR(fdf['LMuCandidate_phi'],fdf['ghCandidate_phi'],fdf['LMuCandidate_eta'],fdf['ghCandidate_eta'])
            deltaRslmughdf = deltaR(fdf['sLMuCandidate_phi'],fdf['ghCandidate_phi'],fdf['sLMuCandidate_eta'],fdf['ghCandidate_eta'])
            deltaRgzghdf   = deltaR(fdf['gzCandidate_phi'],fdf['ghCandidate_phi'],fdf['gzCandidate_eta'],fdf['ghCandidate_eta'])


        #To fill the table
        #Above actually returns a series, making a df to do what I want
        #deltaRzhdfrealdf = deltaRzhdf.to_frame()
        #deltaRzhdfrealdf.columns = ["drzh"]
        #print(deltaRzhdfrealdf)
        #largdrdf = deltaRzhdfrealdf[deltaRzhdfrealdf["drzh"] >=0.8]
        #smaldrdf = deltaRzhdfrealdf[deltaRzhdfrealdf["drzh"] < 0.8]
        #totindr = len(deltaRzhdfrealdf)
        #totbigdr = len(deltaRzhdfrealdf[deltaRzhdfrealdf["drzh"] >=0.8])
        #totsmaldr = len(deltaRzhdfrealdf[deltaRzhdfrealdf["drzh"] < 0.8])
        #print("Total Events         deltaR(Z,h) >= 0.8        deltaR(Z,h) < 0.8")
        #print("    {0}                       {1}                  {2}".format(totindr,totbigdr,totsmaldr))
        #print("    {0}                       {1}                  {2}".format(totindr,round(totbigdr/totindr,2),round(totsmaldr/totindr,2)))
        #print(type(deltaRzhdfrealdf))
        

        #calculate the event weight
        #initialize it. These are all ones if not DY
        eventweights = fdf['event_weight_kf']
        systname = 'btagsystdefaultname'
        if not systl and stype > 0:
            print("    Applying nominal sf")
            eventweights,systname = makeEventWeightSeries(fdf,systl)
        elif not systl and stype <= 0:
            print("    No systematic request and data, so no sf")
            systname = 'btagnom_muidnom'
        elif systl and stype <= 0:
            print("    Attempt at a systematic, but dat, so no sf")
            systname = 'btagnom_muidnom'
        else:
            print("Doing scale factor variations of type: ",systl)
            eventweights,systname = makeEventWeightSeriesWithUncertainties(fdf,systl)

        #print(fdf["NonPrefiringProb"])
        #print(fdf["NonPrefiringProbUp"])
        #print(fdf["NonPrefiringProbDown"])
        #print("all prefire weights sum, final selecs    ",sum(fdf["NonPrefiringProb"]))
        #print("all prefireup weights sum, final selecs  ",sum(fdf["NonPrefiringProbUp"]))
        #print("all prefiredwn weights sum, final selecs ",sum(fdf["NonPrefiringProbDown"]))
        #print("****event weight investiagtion***")
        #print("   basic weights, final selecs           ",sum(eventweights))
        #print("   basic weights w/ prefir, final selecs ",sum(eventweights*fdf["NonPrefiringProb"]))
        #print("   weights with prefireup, final selecs  ",sum(eventweights*fdf["NonPrefiringProbUp"]))
        #print("   weights with prefiredwn, final selecs ",sum(eventweights*fdf["NonPrefiringProbDown"]))



        #prefire weights
        if ((year == 16) or (year == 17)) and stype > 0 and args.prefiredwn:
            print("    Applying dwn prefiring weights!")
            eventweights = eventweights*fdf["NonPrefiringProbDown"]
            systname = systname+"_prfdwn"
        elif ((year == 16) or (year == 17)) and stype > 0 and args.prefireup:
            print("    Applying up prefiring weights!")
            eventweights = eventweights*fdf["NonPrefiringProbUp"]
            systname = systname+"_prfup"
        elif ((year == 16) or (year == 17)) and stype > 0:
            print("    Applying prefiring weights!")
            eventweights = eventweights*fdf["NonPrefiringProb"]
        else:
            print("This is not 2016 or 2017 MC, does not need prefire weights")

        #pile-up weights - real ones -- only in a few files
        #if stype > 0 and args.pileup:
        #    print("    Applying pileup weights!")
        #    eventweights = eventweights*fdf["puWeight"]
        #elif stype > 0 and args.pileupup:
        #    print("    Applying pileup weights - up!")
        #    eventweights = eventweights*fdf["puSysUp"]
        #elif stype > 0 and args.pileupdwn:
        #    print("    Applying pileup weights - down!")
        #    eventweights = eventweights*fdf["puSysDown"]
        #else:
        #    print("This is not MC or you do not want to apply pileup weights")

        #mapped pile-up weights
        f= dict(up3.open(inputfiles[0]))
        if stype > 0 and args.pileupmap and not args.pileupupder and not args.pileupdwnder and not args.pileupupmap and not args.pileupdwnmap:
            print("    Applying pileup weights from mapping!")
            skimtot     = f[b'hnskimed;1'].values[0]
            skimpu      = f[b'hnpu;1'].values[0]
            eventweights = eventweights*fdf["puweight"]*(skimtot/skimpu)
        elif stype > 0 and args.pileupupder:
            skimtot     = f[b'hnskimed;1'].values[0]
            skimpuup    = f[b'hnpuup;1'].values[0]
            print("    Applying pileup weights - up! from mapped value")
            eventweights = eventweights*fdf["puweight_up"]*(skimtot/skimpuup)
            systname = systname+"_puwup"
        elif stype > 0 and args.pileupdwnder:
            skimtot     = f[b'hnskimed;1'].values[0]
            skimpudn    = f[b'hnpudwn;1'].values[0]
            print("    Applying pileup weights - down! from mapped value")
            eventweights = eventweights*fdf["puweight_dwn"]*(skimtot/skimpudn)
            systname = systname+"_puwdwn"
        elif stype > 0 and args.pileupupmap:
            skimtot     = f[b'hnskimed;1'].values[0]
            skimpumapup    = f[b'hnpunumup;1'].values[0]
            print("    Applying pileup weights from up shifted map value")
            systname = systname+"_punumup"
            eventweights = eventweights*fdf["puweightvtx_up"]*(skimtot/skimpumapup)
        elif stype > 0 and args.pileupdwnmap:
            skimtot     = f[b'hnskimed;1'].values[0]
            skimpumapdn    = f[b'hnpunumdwn;1'].values[0]
            print("    Applying pileup weights from down shifted map value")
            eventweights = eventweights*fdf["puweightvtx_dwn"]*(skimtot/skimpumapdn)
            systname = systname+"_punumdwn"
        else:
            print("This is not MC or you do not want to apply pileup weights")




        f = dict(up3.open(inputfiles[0]))
        if stype > 0 and args.qcddwn:
            print("    Applying dwn qcd weights!")
            skimqcddown = f[b'hnskimeddwn;1'].values[0]
            skimtot     = f[b'hnskimed;1'].values[0]
            eventweights = eventweights*fdf["qcdweight_dwn"]*(skimtot/skimqcddown)
            systname = systname+"_qcddwn"
        if stype > 0 and args.qcdup:
            print("    Applying up qcd weights!")
            skimqcdup = f[b'hnskimedup;1'].values[0]
            skimtot     = f[b'hnskimed;1'].values[0]
            eventweights = eventweights*fdf["qcdweight_up"]*(skimtot/skimqcdup)
            
            systname = systname+"_qcdup"
        
        if stype > 0 and args.pdfdwn:
            print("    Applying dwn pdf weights!")
            eventweights = eventweights*fdf["pdfweight_dwn"]
            systname = systname+"_pdfdwn"
        if stype > 0 and args.pdfup:
            print("    Applying up pdf weights!")
            eventweights = eventweights*fdf["pdfweight_up"]
            systname = systname+"_pdfup"

                
        print("The event weights are wrong for general leading lepton dR plots")
        if stype > 0 and  not args.unblind:
            print("    number of passing events straight ",len(fdf))
            print("    number of passing events weighted ",eventweights.sum())
        #print("    max Zp mass estimator: ",fdf['ZPrime_mass_est'].max())

        #lets make some histograms.
        rootfilename  = go.makeOutFile(goodname,'upout_'+region+'_'+jectype+'_'+systname+'_'+btaggr,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))#need to update for btagger

        print("Saving file as ",rootfilename)
        npfilename    = go.makeOutFile(goodname,'totalevents_'+region+'_'+jectype+'_'+systname+'_'+btaggr,'.npy',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        pklfilename   = go.makeOutFile(goodname,'selected_errors_'+region+'_'+jectype+'_'+systname+'_'+btaggr,'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        rootOutFile   = up3.recreate(rootfilename,compression = None)
        npOutFile     = open(npfilename,'wb')

        rootOutFile["h_z_pt"]       = np.histogram(fdf['ZCandidate_pt'],bins=80,range=(0,800),weights=eventweights)
        rootOutFile["h_z_phi"]      = np.histogram(fdf['ZCandidate_phi'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
        rootOutFile["h_z_phiw"]     = np.histogram(fdf['ZCandidate_phi'].map(wrapPhi),bins=30,range=(0,3.14159),weights=eventweights)#wrapped version of phi
        rootOutFile["h_z_eta"]      = np.histogram(fdf['ZCandidate_eta'],bins=30,range=(-5,5),weights=eventweights)
        rootOutFile["h_z_m"]        = np.histogram(fdf['ZCandidate_m'],bins=100,range=(40,140),weights=eventweights)
        rootOutFile["h_h_pt"]       = np.histogram(fdf['hCandidate_pt'],bins=40,range=(200,1200),weights=eventweights)
        rootOutFile["h_h_phi"]      = np.histogram(fdf['hCandidate_phi'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
        rootOutFile["h_h_phiw"]     = np.histogram(fdf['hCandidate_phi'].map(wrapPhi),bins=30,range=(0,3.14159),weights=eventweights)#wrapped version of phi
        rootOutFile["h_h_eta"]      = np.histogram(fdf['hCandidate_eta'],bins=30,range=(-5,5),weights=eventweights)
        rootOutFile["h_h_m"]        = np.histogram(fdf['hCandidate_m'],bins=80,range=(0,400),weights=eventweights)
        rootOutFile["h_h_sd"]       = np.histogram(fdf['hCandidate_sd'],bins=80,range=(0,400),weights=eventweights)
        rootOutFile["h_metxy"]        = np.histogram(fdf['metxycorr'],bins=80,range=(0,2000),weights=eventweights)
        rootOutFile["h_metxy_phi"]    = np.histogram(fdf['metxyphicorr'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
        rootOutFile["h_metxy_phiw"]   = np.histogram(fdf['metxyphicorr'].map(wrapPhi),bins=30,range=(0,3.14159),weights=eventweights)#wrapped version of phi
        rootOutFile["h_met"]        = np.histogram(fdf['metsuable'],bins=80,range=(0,2000),weights=eventweights)
        rootOutFile["h_met_phi"]    = np.histogram(fdf['metphiusable'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
        rootOutFile["h_met_phiw"]   = np.histogram(fdf['metphiusable'].map(wrapPhi),bins=30,range=(0,3.14159),weights=eventweights)#wrapped version of phi
        #rootOutFile["h_zp_jigm"]    = np.histogram(fdf['ZPrime_mass_est'],bins=50,range=(500,5000),weights=eventweights)
        rootOutFile["h_zp_jigm"]    = np.histogram(fdf['ZPrime_mass_est'],bins=130,range=(0,13000),weights=eventweights)
        rootOutFile["h_nd_jigm"]    = np.histogram(fdf['ND_mass_est'],bins=35,range=(100,800),weights=eventweights)
        rootOutFile["h_ns_jigm"]    = np.histogram(fdf['NS_mass_est'],bins=25,range=(0,500),weights=eventweights)
        rootOutFile["h_btag"]       = np.histogram(fdf['hCandidate_'+btaggr],bins=110,range=(0,1.1),weights=eventweights)
        rootOutFile["h_weights"]    = np.histogram(eventweights,bins=40,range=(-1,7))
        rootOutFile["h_dphi_zh"]    = np.histogram(deltaphizhdf,bins=100,range=(0,3.14159),weights=eventweights)
        rootOutFile["h_dphi_zmet"]  = np.histogram(deltaphizmetdf,bins=100,range=(0,3.14159),weights=eventweights)
        rootOutFile["h_dphi_hmet"]  = np.histogram(deltaphihmetdf,bins=100,range=(0,3.14159),weights=eventweights)
        rootOutFile["h_dr_zh"]      = np.histogram(deltaRzhdf,bins=30,range=(0,6),weights=eventweights)
        rootOutFile["h_totmom"]     = np.histogram(tscalmom,bins=250,range=(0,2500),weights=eventweights)
        rootOutFile["h_metfrac"]    = np.histogram(metscalmomfrac,bins=100,range=(0,1),weights=eventweights)
        rootOutFile["h_zptfrac"]    = np.histogram(zptscalmomfrac,bins=100,range=(0,1),weights=eventweights)
        rootOutFile["h_hptfrac"]    = np.histogram(hptscalmomfrac,bins=100,range=(0,1),weights=eventweights)

        zpterrs      = boostUnc(fdf['ZCandidate_pt'],eventweights,80,0,800)
        zetaerrs     = boostUnc(fdf['ZCandidate_eta'],eventweights,30,-5,5)
        zphierrs     = boostUnc(fdf['ZCandidate_phi'],eventweights,30,-3.14159,3.14159)
        zphiwerrs    = boostUnc(fdf['ZCandidate_phi'].map(wrapPhi),eventweights,30,0,3.14159)
        zmerrs       = boostUnc(fdf['ZCandidate_m'],eventweights,100,40,140)
        hpterrs      = boostUnc(fdf['hCandidate_pt'],eventweights,40,200,1200)
        hetaerrs     = boostUnc(fdf['hCandidate_eta'],eventweights,30,-5,5)
        hphierrs     = boostUnc(fdf['hCandidate_phi'],eventweights,30,-3.14159,3.14159)
        hphiwerrs    = boostUnc(fdf['hCandidate_phi'].map(wrapPhi),eventweights,30,0,3.14159)
        hmerrs       = boostUnc(fdf['hCandidate_m'],eventweights,80,0,400)
        hsderrs      = boostUnc(fdf['hCandidate_sd'],eventweights,80,0,400)
        meterrs      = boostUnc(fdf['metsuable'],eventweights,80,0,2000)
        metphierrs   = boostUnc(fdf['metphiusable'],eventweights,30,-3.14159,3.14159)
        metphiwerrs  = boostUnc(fdf['metphiusable'].map(wrapPhi),eventweights,30,0,3.14159)
        meterrsxy      = boostUnc(fdf['metxycorr'],eventweights,80,0,2000)
        metphierrsxy   = boostUnc(fdf['metxyphicorr'],eventweights,30,-3.14159,3.14159)
        metphiwerrsxy  = boostUnc(fdf['metxyphicorr'].map(wrapPhi),eventweights,30,0,3.14159)
        zpjigerrs    = boostUnc(fdf['ZPrime_mass_est'],eventweights,130,0,13000)
        ndjigerrs    = boostUnc(fdf['ND_mass_est'],eventweights,35,100,800)
        nsjigerrs    = boostUnc(fdf['NS_mass_est'],eventweights,25,0,500)
        btagerrs     = boostUnc(fdf['hCandidate_'+btaggr],eventweights,110,0,1.1)
        dphizherrs   = boostUnc(deltaphizhdf,eventweights,100,0,3.14159)
        dphizmeterrs   = boostUnc(deltaphizmetdf,eventweights,100,0,3.14159)
        dphihmeterrs   = boostUnc(deltaphihmetdf,eventweights,100,0,3.14159)
        drzherrs       = boostUnc(deltaRzhdf,eventweights,30,0,6)
        tscalerrs      = boostUnc(tscalmom,eventweights,250,0,2500)
        metscalerrs    = boostUnc(metscalmomfrac,eventweights,100,0,1)
        zptscalerrs    = boostUnc(zptscalmomfrac,eventweights,100,0,1)
        hptscalerrs    = boostUnc(hptscalmomfrac,eventweights,100,0,1)
        
        #drlmuherrs     = boostUnc(deltaRlmuhdf,eventweights,30,0,6)
        #drslmuherrs    = boostUnc(deltaRslmuhdf,eventweights,30,0,6)
        #drslmulmuerrs  = boostUnc(deltaRslmulmudf,eventweights,30,0,6)
        #lmupterrs      = boostUnc(fdf['LMuCandidate_pt'],eventweights,50,0,500)
        #smupterrs      = boostUnc(fdf['sLMuCandidate_pt'],eventweights,50,0,500)
        #lmuphierrs     = boostUnc(fdf['LMuCandidate_phi'],eventweights,30,-3.14159,3.14159)
        #smuphierrs     = boostUnc(fdf['sLMuCandidate_phi'],eventweights,30,-3.14159,3.14159)
        #lmuetaerrs     = boostUnc(fdf['LMuCandidate_eta'],eventweights,30,-5,5)
        #smuetaerrs     = boostUnc(fdf['sLMuCandidate_eta'],eventweights,30,-5,5)
        
        unc_arrays = [zpterrs,
                      zetaerrs,
                      zphierrs,
                      zphiwerrs,
                      zmerrs,
                      hpterrs,
                      hetaerrs,
                      hphierrs,
                      hphiwerrs,
                      hmerrs,
                      hsderrs,
                      meterrs,
                      metphierrs,
                      metphiwerrs,
                      meterrsxy,
                      metphierrsxy,
                      metphiwerrsxy,
                      zpjigerrs,
                      ndjigerrs,
                      nsjigerrs,
                      btagerrs,
                      dphizherrs,
                      dphizmeterrs,
                      dphihmeterrs,
                      drzherrs,
                      tscalerrs,
                      metscalerrs,
                      zptscalerrs,
                      hptscalerrs,
        ]

        unc_names = ['h_z_pt',
                     'h_z_eta',
                     'h_z_phi',
                     'h_z_phiw',
                     'h_z_m',
                     'h_h_pt',
                     'h_h_eta',
                     'h_h_phi',
                     'h_h_phiw',
                     'h_h_m',
                     'h_h_sd',
                     'h_met',
                     'h_met_phi',
                     'h_met_phiw',
                     'h_metxy',
                     'h_metxy_phi',
                     'h_metxy_phiw',
                     'h_zp_jigm',
                     'h_nd_jigm',
                     'h_ns_jigm',
                     'h_btag',
                     'h_dphi_zh',
                     'h_dphi_zmet',
                     'h_dphi_hmet',
                     'h_dr_zh',
                     'h_totmom',
                     'h_metfrac',
                     'h_zptfrac',
                     'h_hptfrac',
        ]
        
        if (stype > 0 and channel != "emu"):
            rootOutFile["h_dr_lmu_gh"] = np.histogram(deltaRlmughdf,bins=30,range=(0,6),weights=eventweights)
            rootOutFile["h_dr_slmu_gh"] = np.histogram(deltaRslmughdf,bins=30,range=(0,6),weights=eventweights)
            rootOutFile["h_dr_gz_gh"] = np.histogram(deltaRgzghdf,bins=30,range=(0,6),weights=eventweights)
            drlmugherrs    = boostUnc(deltaRlmughdf,eventweights,30,0,6)
            drslmugherrs   = boostUnc(deltaRslmughdf,eventweights,30,0,6)
            drgzgherrs     = boostUnc(deltaRgzghdf,eventweights,30,0,6)

            addhists = [drlmugherrs,drslmugherrs,drgzgherrs]
            addnames = ["h_dr_lmu_gh","h_dr_slmu_gh","h_dr_gz_gh"]
            unc_arrays.extend(addhists)
            unc_names.extend(addnames)

        if (channel != "emu"):
            rootOutFile["h_dr_lmuh"]    = np.histogram(deltaRlmuhdf,bins=30,range=(0,6),weights=eventweights)
            rootOutFile["h_dr_slmuh"]   = np.histogram(deltaRslmuhdf,bins=30,range=(0,6),weights=eventweights)
            rootOutFile["h_dr_slmulmu"] = np.histogram(deltaRslmulmudf,bins=30,range=(0,6),weights=eventweights)
            rootOutFile["h_LMu_pt"]     = np.histogram(fdf['LMuCandidate_pt'],bins=200,range=(0,2000),weights=eventweights)
            rootOutFile["h_sLMu_pt"]     = np.histogram(fdf['sLMuCandidate_pt'],bins=200,range=(0,2000),weights=eventweights)
            rootOutFile["h_LMu_phi"]    = np.histogram(fdf['LMuCandidate_phi'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_sLMu_phi"]    = np.histogram(fdf['sLMuCandidate_phi'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_LMu_eta"]    = np.histogram(fdf['LMuCandidate_eta'],bins=30,range=(-5,5),weights=eventweights)
            rootOutFile["h_sLMu_eta"]    = np.histogram(fdf['sLMuCandidate_eta'],bins=30,range=(-5,5),weights=eventweights)
            drlmuherrs     = boostUnc(deltaRlmuhdf,eventweights,30,0,6)
            drslmuherrs    = boostUnc(deltaRslmuhdf,eventweights,30,0,6)
            drslmulmuerrs  = boostUnc(deltaRslmulmudf,eventweights,30,0,6)
            lmupterrs      = boostUnc(fdf['LMuCandidate_pt'],eventweights,200,0,2000)
            smupterrs      = boostUnc(fdf['sLMuCandidate_pt'],eventweights,200,0,2000)
            lmuphierrs     = boostUnc(fdf['LMuCandidate_phi'],eventweights,30,-3.14159,3.14159)
            smuphierrs     = boostUnc(fdf['sLMuCandidate_phi'],eventweights,30,-3.14159,3.14159)
            lmuetaerrs     = boostUnc(fdf['LMuCandidate_eta'],eventweights,30,-5,5)
            smuetaerrs     = boostUnc(fdf['sLMuCandidate_eta'],eventweights,30,-5,5)
            addhists  = [drlmuherrs,
                         drslmuherrs,
                         drslmulmuerrs,
                         lmupterrs,
                         smupterrs,
                         lmuphierrs,
                         smuphierrs,
                         lmuetaerrs,
                         smuetaerrs]
            addnames  = ['h_dr_lmuh',
                         'h_dr_slmuh',
                         'h_dr_slmulmu',
                         'h_LMu_pt',
                         'h_sLMu_pt',
                         'h_LMu_phi',
                         'h_sLMu_phi',
                         'h_LMu_eta',
                         'h_sLMu_eta']
            unc_arrays.extend(addhists)
            unc_names.extend(addnames)
            
        if (channel == "emu"):
            rootOutFile["h_dr_leadleph"]  = np.histogram(deltaRllephdf,bins=30,range=(0,6),weights=lepdrweights)
            rootOutFile["h_dr_sleadleph"] = np.histogram(deltaRslephdf,bins=30,range=(0,6),weights=lepdrweights)
            rootOutFile["h_dr_leps"]      = np.histogram(deltaRsleplepdf,bins=30,range=(0,6),weights=lepdrweights)
            rootOutFile["h_leadlep_pt"]   = np.histogram(leadpt,bins=200,range=(0,2000),weights=eventweights)
            rootOutFile["h_sleadlep_pt"]  = np.histogram(sleadpt,bins=200,range=(0,2000),weights=eventweights)
            rootOutFile["h_leadlep_phi"]  = np.histogram(leadphi,bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_sleadlep_phi"] = np.histogram(sleadphi,bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_leadlep_eta"]  = np.histogram(leadeta,bins=30,range=(-5,5),weights=eventweights)
            rootOutFile["h_sleadlep_eta"] = np.histogram(sleadeta,bins=30,range=(-5,5),weights=eventweights)
            rootOutFile["h_electron_pt"]  = np.histogram(allelpt,bins=50,range=(0,500),weights=eventweights)
            rootOutFile["h_electron_phi"]  = np.histogram(allelphi,bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_electron_eta"]  = np.histogram(alleleta,bins=30,range=(-5,5),weights=eventweights)
            rootOutFile["h_muon_pt"]  = np.histogram(allmupt,bins=200,range=(0,2000),weights=eventweights)
            rootOutFile["h_muon_phi"]  = np.histogram(allmuphi,bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_muon_eta"]  = np.histogram(allmueta,bins=30,range=(-5,5),weights=eventweights)
            #new kinematic plots
            rootOutFile["h_lmuon_pt"]  = np.histogram(leadmuc['LMuCandidate_pt'],bins=200,range=(0,2000),weights=lmuweights)
            rootOutFile["h_lmuon_phi"]  = np.histogram(leadmuc['LMuCandidate_phi'],bins=30,range=(-3.14159,3.14159),weights=lmuweights)
            rootOutFile["h_lmuon_eta"]  = np.histogram(leadmuc['LMuCandidate_eta'],bins=30,range=(-5,5),weights=lmuweights)
            rootOutFile["h_slelectron_pt"]  = np.histogram(leadmuc['sLEleCandidate_pt'],bins=50,range=(0,500),weights=lmuweights)
            rootOutFile["h_slelectron_phi"]  = np.histogram(leadmuc['sLEleCandidate_phi'],bins=30,range=(-3.14159,3.14159),weights=lmuweights)
            rootOutFile["h_slelectron_eta"]  = np.histogram(leadmuc['sLEleCandidate_eta'],bins=30,range=(-5,5),weights=lmuweights)
            rootOutFile["h_lelectron_pt"]  = np.histogram(leadelc['LEleCandidate_pt'],bins=50,range=(0,500),weights=lelweights)
            rootOutFile["h_lelectron_phi"]  = np.histogram(leadelc['LEleCandidate_phi'],bins=30,range=(-3.14159,3.14159),weights=lelweights)
            rootOutFile["h_lelectron_eta"]  = np.histogram(leadelc['LEleCandidate_eta'],bins=30,range=(-5,5),weights=lelweights)
            rootOutFile["h_slmuon_pt"]  = np.histogram(leadelc['sLMuCandidate_pt'],bins=200,range=(0,2000),weights=lelweights)
            rootOutFile["h_slmuon_phi"]  = np.histogram(leadelc['sLMuCandidate_phi'],bins=30,range=(-3.14159,3.14159),weights=lelweights)
            rootOutFile["h_slmuon_eta"]  = np.histogram(leadelc['sLMuCandidate_eta'],bins=30,range=(-5,5),weights=lelweights)
            rootOutFile["h_dr_leadmuonh"]  = np.histogram(deltaRemuLeadMuh,bins=30,range=(0,6),weights=lmuweights)
            rootOutFile["h_dr_subleadingeleh"]  = np.histogram(deltaRemuSubLeadElh,bins=30,range=(0,6),weights=lmuweights)
            rootOutFile["h_dr_leadeleh"]  = np.histogram(deltaRemuLeadElh,bins=30,range=(0,6),weights=lelweights)
            rootOutFile["h_dr_subleadingmuh"]  = np.histogram(deltaRemuSubLeadMuh,bins=30,range=(0,6),weights=lelweights)

            drllherrs   = boostUnc(deltaRllephdf,lepdrweights,30,0,6)
            drsllherrs  = boostUnc(deltaRslephdf,lepdrweights,30,0,6)
            drsllllerrs = boostUnc(deltaRsleplepdf,lepdrweights,30,0,6)
            llpterrs    = boostUnc(leadpt,eventweights,50,0,500)
            sllpterrs   = boostUnc(sleadpt,eventweights,50,0,500)
            llphierrs   = boostUnc(leadphi,eventweights,30,-3.14159,3.14159)
            sllphierrs  = boostUnc(sleadphi,eventweights,30,-3.14159,3.14159)
            lletaerrs   = boostUnc(leadeta,eventweights,30,-5,5)
            slletaerrs  = boostUnc(sleadeta,eventweights,30,-5,5)
            allelpterrs = boostUnc(allelpt,eventweights,50,0,500)
            allelphierrs = boostUnc(allelphi,eventweights,30,-3.14159,3.14159)
            alleletaerrs = boostUnc(alleleta,eventweights,30,-5,5)
            allmupterrs = boostUnc(allmupt,eventweights,50,0,500)
            allmuphierrs = boostUnc(allmuphi,eventweights,30,-3.14159,3.14159)
            allmuetaerrs = boostUnc(allmueta,eventweights,30,-5,5)
            emulmupterrs  = boostUnc(leadmuc['LMuCandidate_pt'],lmuweights,200,0,2000)
            emulmuphierrs = boostUnc(leadmuc['LMuCandidate_phi'],lmuweights,30,-3.14159,3.14159)
            emulmuetaerrs = boostUnc(leadmuc['LMuCandidate_eta'],lmuweights,30,-5,5)
            emuslelpterrs = boostUnc(leadmuc['sLEleCandidate_pt'],lmuweights,50,0,500)
            emuslelphierrs = boostUnc(leadmuc['sLEleCandidate_phi'],lmuweights,30,-3.14159,3.14159)
            emusleletaerrs = boostUnc(leadmuc['sLEleCandidate_eta'],lmuweights,30,-5,5)
            emulelpterrs = boostUnc(leadelc['LEleCandidate_pt'],lelweights,50,0,500)
            emulelphierrs = boostUnc(leadelc['LEleCandidate_phi'],lelweights,30,-3.14159,3.14159)
            emuleletaerrs = boostUnc(leadelc['LEleCandidate_eta'],lelweights,30,-5,5)
            emuslmupterrs = boostUnc(leadelc['sLMuCandidate_pt'],lelweights,200,0,2000)
            emuslmuphierrs = boostUnc(leadelc['sLMuCandidate_phi'],lelweights,30,-3.14159,3.14159)
            emuslmuetaerrs = boostUnc(leadelc['sLMuCandidate_eta'],lelweights,30,-5,5)
            emudrlmuherrs  = boostUnc(deltaRemuLeadMuh,lmuweights,30,0,6)
            emudrslelherrs = boostUnc(deltaRemuSubLeadElh,lmuweights,30,0,6)
            emudrleleherrs = boostUnc(deltaRemuLeadElh,lelweights,30,0,6)
            emudrslmuherrs = boostUnc(deltaRemuSubLeadMuh,lelweights,30,0,6)


            addhists  = [drllherrs,
                         drsllherrs,
                         drsllllerrs,
                         llpterrs,
                         sllpterrs,
                         llphierrs,
                         sllphierrs,
                         lletaerrs,
                         slletaerrs,
                         allelpterrs,
                         allelphierrs,
                         alleletaerrs,
                         allmupterrs,
                         allmuphierrs,
                         allmuetaerrs,
                         emulmupterrs,
                         emulmuphierrs,
                         emulmuetaerrs,
                         emuslelpterrs,
                         emuslelphierrs,
                         emusleletaerrs,
                         emulelpterrs,
                         emulelphierrs,
                         emuleletaerrs,
                         emuslmupterrs,
                         emuslmuphierrs,
                         emuslmuetaerrs,
                         emudrlmuherrs ,
                         emudrslelherrs,
                         emudrleleherrs,
                         emudrslmuherrs,
            ]
            
            addnames  = ["h_dr_leadleph",
                         "h_dr_sleadleph",
                         "h_dr_leps",
                         "h_leadlep_pt",
                         "h_sleadlep_pt",
                         "h_leadlep_phi",
                         "h_sleadlep_phi",
                         "h_leadlep_eta",
                         "h_sleadlep_eta",
                         "h_electron_pt",
                         "h_electron_phi",
                         "h_electron_eta",
                         "h_muon_pt",
                         "h_muon_phi",
                         "h_muon_eta",
                         "h_lmuon_pt",
                         "h_lmuon_phi",
                         "h_lmuon_eta",
                         "h_slelectron_pt",
                         "h_slelectron_phi",
                         "h_slelectron_eta",
                         "h_lelectron_pt",
                         "h_lelectron_phi",
                         "h_lelectron_eta",
                         "h_slmuon_pt",
                         "h_slmuon_phi",
                         "h_slmuon_eta",
                         "h_dr_leadmuonh",
                         "h_dr_subleadingeleh",
                         "h_dr_leadeleh",
                         "h_dr_subleadingmuh"
            ]
            unc_arrays.extend(addhists)
            unc_names.extend(addnames)
            

        max_length = len(max(unc_arrays,key = lambda ar : len(ar)))
        pad_arrays = [np.pad(arr,(0,max_length - len(arr)),'constant') for arr in unc_arrays]
        all_unc    = np.column_stack(pad_arrays)
        uncdf      = pd.DataFrame(all_unc,columns=unc_names)
        uncdf.to_pickle("./"+pklfilename)
    
        #Book Keeping
        #f = up3.open(inputfiles[0])
        np.save(npOutFile,np.array([f[b'hnorigevnts;1'].values[0]]))
        rootOutFile["hnevents"]      = str(f[b'hnorigevnts;1'].values[0])
        rootOutFile["hnevents_pMET"] = str((metdf['event_weight_kf']*metdf['event_weight_btag']*metdf['event_weight_muid']).sum())#str(len(metdf))
        rootOutFile["hnevents_pZ"]   = str((zptdf['event_weight_kf']*zptdf['event_weight_btag']*zptdf['event_weight_muid']).sum())#str(len(zptdf))
        rootOutFile["hnevents_ph"]   = str((hptdf['event_weight_kf']*hptdf['event_weight_btag']*hptdf['event_weight_muid']).sum())#str(len(hptdf))
        rootOutFile["hnevents_sb"]   = str((sbdf['event_weight_kf']*sbdf['event_weight_btag']*sbdf['event_weight_muid']).sum())#str(len(sbdf))
        rootOutFile["hnevents_btag"] = str((btdf['event_weight_kf']*btdf['event_weight_btag']*btdf['event_weight_muid']).sum())#str(len(btdf))

        #print("Passing events sb events, not sf              : ",len(sbdf))
        #print("Passing events sb events, sum of event weights: ",

        if stype > 0:
            if sr:
                rootOutFile["hnevents_sr"]   = str((srdf['event_weight_kf']*srdf['event_weight_btag']).sum())#str(len(srdf))
        print("Unweighted initial events :     ",len(events))
        print("Weighted Events, passing Z pT:  ",str((zptdf['event_weight_kf']*zptdf['event_weight_btag']*zptdf['event_weight_muid']).sum()))
        print("Weighted Events, passing MET:   ",str((metdf['event_weight_kf']*metdf['event_weight_btag']*metdf['event_weight_muid']).sum()))
        print("Weighted Events, passing H pT:  ",str((hptdf['event_weight_kf']*hptdf['event_weight_btag']*hptdf['event_weight_muid']).sum()))
        print("Weighted Events, passing btag:  ",str((btdf['event_weight_kf']*btdf['event_weight_btag']*btdf['event_weight_muid']).sum()))
        if "emu" in channel:
            print("Weighted Events, signal region: ",str((srdf['event_weight_kf']*srdf['event_weight_btag']*srdf['event_weight_muid']).sum()))
            print("Weighted Events, sideband:      ",str((sbdf['event_weight_kf']*sbdf['event_weight_btag']*sbdf['event_weight_muid']).sum()))

                
        #plots to understand stuff.
        #higgspt = fdf['LMuCandidate_pt']

        #print(fdf['LMuCandidate_pt'])
        #print(fdf['LMuCandidate_ptunc'])
        #print(fdf['LMuCandidate_ptunc']/fdf['LMuCandidate_pt'])
        #higgspt = fdf['sLMuCandidate_ptunc']/fdf['sLMuCandidate_pt']
        #higgspt = fdf['hCandidate_pt']
        #higgspt  = fdf['metxycorr']
        #higgspt  = fdf['ZCandidate_pt']
        #zpmest  = fdf['ZPrime_mass_est']
        #zpmest  = fdf['sLMuCandidate_pt']
        #scale = go.findScale(35000,137.6,1)
        #wscale = [scale for i in range(len(higgspt))]
        #wscale = np.array(wscale)
        #hboost = boosthist(hptscalmonfrac,eventweights,100,0,1000)
        #hboost = boosthist(fdf['hCandidate_pt'],eventweights,40,200,1200)
        #hboost = boosthist(fdf['metxycorr'],eventweights,80,0,2000)
        #hboost = boosthist(fdf['ZCandidate_pt'],eventweights,80,0,800)
        #hboost = boosthist(higgspt,eventweights,200,0,2000)#when muon pT alone y-axis
        #hboost  = boosthist(higgspt,eventweights,10,0,0.2)
        #hedg   = hboost.axes.edges
        #zpbst  = boosthist(fdf['ZPrime_mass_est'],eventweights,65,0,13000)#when RJR along x-axis
        #zpbst = boosthist(zpmest,eventweights,200,0,2000)
        #zpedg = zpbst.axes.edges
        #h2d, xedges, yedges = np.histogram2d(zpmest,higgspt,(zpedg[0],hedg[0]),weights=wscale)
        #h2d, xedges, yedges = np.histogram2d(zpmest,higgspt,(zpedg[0],hedg[0]))

        


        #h2d = h2d.T
        #xaxis,yaxis = np.meshgrid(xedges,yedges)
        #plt.pcolormesh(xaxis,yaxis,h2d)
        #plt.colorbar()
        #plt.xlabel("RJR Z Prime Mass Estimator")
        #plt.xlabel("Subleading Muon Momentum")
        #plt.ylabel("Leading Muon Momentum")
        #plt.ylabel("Relative Muon Momentum Unctertainty")
        #plt.title("Zp 5500 ND 1800 NS 200, signal region, xs = 1 fb")
        #plt.title("Zp 5500 ND 1800 NS 200, signal region, no xs scale")
        #plt.show()

        #highvals = [500.0,750.0,1000.0]

        #highzps = zpedg[0][zpedg[0] >= 5000]
        #highzps = highzps[highzps < 9800.0]

        #highval = 1000.0

        #for highval in highvals:
        #    highfracs = []
        #    numevents = []
        #    numhighs  = []
        #    for zp in highzps:
        #        zpdf = fdf[fdf['ZPrime_mass_est'] >= zp]
        #        trunzpdf = zpdf[zpdf['ZPrime_mass_est'] < 10000.0]
        #        himu = len(trunzpdf[trunzpdf['LMuCandidate_pt'] >= highval])
        #        numevents.append(len(trunzpdf))
        #        numhighs.append(himu)
        #        if len(trunzpdf) < 1:
        #            highfracs.append(0)
        #        else:
        #            highfrac = himu/len(trunzpdf)
        #            highfracs.append(highfrac)
                
        #    fig, axs = plt.subplots(2)

        #    scales = [x*scale for x in numhighs]
        #    axs[0].plot(highzps,scales)
        #    axs[0].set_title('Number of leading muons w/ pT > {0} GeV'.format(highval))
        #    axs[1].plot(highzps,highfracs)
        #    axs[1].set_title('Fraction of total events w/ leading muon pT > {0} GeV'.format(highval))
        #    axs[1].set(xlabel='Lower edge of RJR Estimator Overflow Bin')
        #    plt.ylim([0.0,1])

        #    for ax in axs.flat:
        #        ax.label_outer()

        #    figname = go.makeOutFile("leadingmuoninfo_vs_jigzp_scaled",'highptmucut'+str(int(highval))+'_'+region+'_'+jectype+'_'+systname+'_'+btaggr,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp))

        
        #    plt.savefig(figname)
        #    plt.close()
        
        
