import uproot as up4
import uproot3 as up3
import pandas as pd
import numpy as np
import boost_histogram as bh
import argparse
import glob
import gecorg_test as go


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
    parser.add_argument("-a","--alphar",type=bool,help = "alpha ratio regions?")
    parser.add_argument("-c","--chan",type=str)
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

    print(sampname)

    inputfiles = glob.glob(topdir+'/'+sampname+'*_topiary_'+channel+'*.root')
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
        ]

        if (stype > 0):
            mcbranches = [b'ghCandidate_*',
                          b'gzCandidate_*']
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
            
            mtree = up3.open(muf[0])['PreSelection;1']
            mevents = mtree.pandas.df(branches=branches)
            etree = up3.open(euf[0])['PreSelection;1']
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
            inputfiles = glob.glob(topdir+'/'+samp+'*_topiary_'+channel+'_*'+jectype+'*.root')
            print("Doing selections on:")
            print("    ",inputfiles[:1])
            print("    ",samp)
            #events = up3.pandas.iterate(inputfiles[:1],'PreSelection;1',branches=branches)
            tree = up3.open(inputfiles[0])['PreSelection;1']
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
            
        zptdf  = events[events['ZCandidate_pt'] > zptcut]
        #metdf   = zptdf[zptdf['metsuable'] > metcut]
        metdf   = zptdf[zptdf['metxycorr'] > metcut]
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

        if ((year == 16) or (year == 17)) and stype > 0:
            print("    Applying prefiring weights!")
            eventweights = eventweights*fdf["NonPrefiringProb"]
        if ((year == 16) or (year == 17)) and stype > 0 and args.prefireup:
            print("    Applying up prefiring weights!")
            eventweights = eventweights*fdf["NonPrefiringProbUp"]
            systname = systname+"_prfup"
        if ((year == 16) or (year == 17)) and stype > 0 and args.prefiredwn:
            print("    Applying dwn prefiring weights!")
            eventweights = eventweights*fdf["NonPrefiringProbDown"]
            systname = systname+"_prfdwn"

                
        print("    number of passing events straight ",len(fdf))
        print("    number of passing events weighted, before PDF scales ",eventweights.sum())
        #print("    max Zp mass estimator: ",fdf['ZPrime_mass_est'].max())

        #lets make some histograms.
        rootfilename  = go.makeOutFile(goodname,'upout_'+region+'_'+jectype+'_'+systname+'_'+btaggr,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))#need to update for btagger

        print("Saving file as ",rootfilename)
        npfilename    = go.makeOutFile(goodname,'totalevents_'+region+'_'+jectype+'_'+systname+'_'+btaggr,'.npy',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        pklfilename   = go.makeOutFile(goodname,'selected_errors_'+region+'_'+jectype+'_'+systname+'_'+btaggr,'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        rootOutFile   = up3.recreate(rootfilename,compression = None)
        npOutFile     = open(npfilename,'wb')

        #Do the QCD scale weights.
        unc_arrays = []
        unc_names = []

        for scale in range(9):
            print("checking scale ",scale)
            scaleseries = fdf.xs(scale,level=1,axis=0)#scales are 1 level into a multiindex
            

        
            rootOutFile["h_zp_jigm"+str(scale)]    = np.histogram(fdf['ZPrime_mass_est'],bins=130,range=(0,13000),weights=eventweights*scaleseries)
            zpjigerrs    = boostUnc(fdf['ZPrime_mass_est'],eventweights*scaleseries,130,0,13000)
            unc_arrays.append(zpjigerrs)
            unc_names.append("h_zp_jigm"+str(scale))
        



        max_length = len(max(unc_arrays,key = lambda ar : len(ar)))
        pad_arrays = [np.pad(arr,(0,max_length - len(arr)),'constant') for arr in unc_arrays]
        all_unc    = np.column_stack(pad_arrays)
        uncdf      = pd.DataFrame(all_unc,columns=unc_names)
        uncdf.to_pickle("./"+pklfilename)
    
        #Book Keeping
        f = up3.open(inputfiles[0])
        np.save(npOutFile,np.array([f['hnorigevnts'].values[0]]))
        rootOutFile["hnevents"]      = str(f['hnorigevnts'].values[0])
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

                
