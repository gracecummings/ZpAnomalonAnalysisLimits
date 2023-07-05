import uproot as up4
import uproot3 as up3
import pandas as pd
import numpy as np
import boost_histogram as bh
import argparse
import glob
import gecorg_test as go


parser = argparse.ArgumentParser()

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
    channel = args.chan
    topdir = args.directory

    inputfiles = glob.glob(topdir+'/'+sampname+'*_topiary_'+channel+'*.root')
    #inputfiles = glob.glob('analysis_output_ZpAnomalon/2022-03-14/noLeadingReq/Run2016H-17Jul2018-v1.SingleElectron_topiary_emu_systnominal_elecTrigNoLeadingRequirement_Zptcut0.0_Hptcut250.0_metcut0.0_btagwp0.0*')
    #inputfiles = glob.glob('analysis_output_ZpAnomalon/2022-03-28/Autumn18.TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_topiary_mumu_systnominal_Zptcut0.0_Hptcut250.0_metcut0.0_btagwp0.0*')
    innames,jecs = getContentsOfSample(inputfiles)
    samp = innames[0]
    print("The parameters you are invoking: ")
    print("    ",innames,jecs)

    for i,jectype in enumerate(jecs):
        stype,year = go.sampleType(samp)

        if valid:
            print("    Doing validation of alpha method cuts")
            if sr:
                print("    'signalr' labeled events are in soft drop mass bands (55,70]")
            else:
                print("    'sideband' labeled events are in soft drop mass bands [30,55],[150,5000]")

        metstr = ''
        branches = [b'ZCandidate_*',
                    #b'RunNum',
                    #b'LumiBlockNum',
                    #b'EvtNum',
                    b'hCandidate_*',
                    b'metsuable',
                    b'metphiusable',
                    b'ZPrime_mass_est',
                    b'ND_mass_est',
                    b'NS_mass_est',
                    b'event_weight*',
                    #b'LMuCandidate_*',
                    #b'sLMuCandidate_*',
                    #b'LEleCandidate_*',
                    #b'sLEleCandidate_*',
                    #b'event_weight_kf',
                    #b'event_weight_btag'
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
            #print(mevents)
            #print(eevents)
            print("Number of events in muon dataset ",len(mevents))
            print("Number of events in elec dataset ",len(eevents))
            
            frames = [mevents,eevents]
            mixdf = pd.concat(frames)
            print("Number of events in straight mixed set ",len(mixdf))
            print("dropping duplicates")
            mixdf.drop_duplicates(subset = ['RunNum','LumiBlockNum','EvtNum'])
            print("Number of events after dropping duplicates ",len(mixdf))
            events = mixdf
            #print(events)
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
        print("    Concerning systematics:")
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
        #sddf   = events[events['hCandidate_sd'] > sdmcut]
        #metdf  = sddf[sddf['metsuable'] > metcut]
        #zptdf  = metdf[metdf['ZCandidate_pt'] > zptcut]
        #hptdf  = zptdf[zptdf['hCandidate_pt'] > hptcut]
        #btdf   = hptdf[hptdf['hCandidate_'+btaggr] > float(btagwp)]

        #new cut order
        zptdf  = events[events['ZCandidate_pt'] > zptcut]
        metdf   = zptdf[zptdf['metsuable'] > metcut]
        hptdf  = metdf[metdf['hCandidate_pt'] > hptcut]
        btdf   = hptdf[hptdf['hCandidate_'+btaggr] > float(btagwp)]

        #Actual Analysis
        if not valid:
            srup   = btdf[btdf['hCandidate_sd'] > 70.]
            bldf   = srup[srup['hCandidate_sd'] < 150.]#full blinded region
            srdf   = bldf[bldf['hCandidate_sd'] > 110.]#Higgs Peak
            lowsbh = btdf[btdf['hCandidate_sd'] <= 70.]
            lowsb  = lowsbh[lowsbh['hCandidate_sd'] > 30.]
            #lowsb = btdf[btdf['hCandidate_sd'] <= 70.]#old way
            highsb = btdf[btdf['hCandidate_sd'] >= 150.]
            sbdf   = pd.concat([lowsb,highsb])
            totdf = pd.concat([sbdf,srdf])
            
        #Validation region for alpha method for DY
        if valid:
            srup   = btdf[btdf['hCandidate_sd'] > 55.]
            bldf   = srup[srup['hCandidate_sd'] < 150.]#full blinded region
            srdf   = bldf[bldf['hCandidate_sd'] <= 70.]#Validation region
            #lowsb  = btdf[btdf['hCandidate_sd'] <= 55.]
            lowsbh = btdf[btdf['hCandidate_sd'] <= 55.]
            lowsb  = lowsbh[lowsbh['hCandidate_sd'] > 30.]
            highsb = btdf[btdf['hCandidate_sd'] >= 150.]
            sbdf   = pd.concat([lowsb,highsb])
            totdf = pd.concat([sbdf,srdf])

        region = "sideband"
        if not valid:
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

        #print("number of btag passing events ",len(btdf))
        
        #calculated quantities
        deltaphizhdf   = deltaPhi(fdf['ZCandidate_phi'],fdf['hCandidate_phi'])
        deltaphizmetdf = deltaPhi(fdf['ZCandidate_phi'],fdf['metphiusable'])
        deltaphihmetdf = deltaPhi(fdf['hCandidate_phi'],fdf['metphiusable'])
        deltaRzhdf     = deltaR(fdf['ZCandidate_phi'],fdf['hCandidate_phi'],fdf['ZCandidate_eta'],fdf['hCandidate_eta'])

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

            print(" Length of leading muon df in final selection ",len(leadmu))
            print(" Length of leading electron df in",len(leadel))

            #get the event weights in order for these reordered series
            lepordrkfw         = leadmu['event_weight_kf'].append(leadel['event_weight_kf'])
            lepordrbtagw       = leadmu['event_weight_btag'].append(leadel['event_weight_btag'])
            lepordrbtaguncdwnw = leadmu['event_weight_btaguncdwn'].append(leadel['event_weight_btaguncdwn'])
            lepordrbtaguncupw  = leadmu['event_weight_btaguncup'].append(leadel['event_weight_btaguncup'])
            lepordrmuidw       = leadmu['event_weight_muid'].append(leadel['event_weight_muid'])
            lepordrmuidwuncdwn = leadmu['event_weight_muiduncdwn'].append(leadel['event_weight_muiduncdwn'])
            lepordrmuidwuncup  = leadmu['event_weight_muiduncup'].append(leadel['event_weight_muiduncup'])

            #gotta make the hCandidate stuff with the same variables
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
            deltaRslephdf   = deltaR(sleadphi,leadinglephcandphi,sleadeta,leadinglephcandeta)
            deltaRsleplepdf = deltaR(sleadphi,leadphi,sleadeta,leadeta)
            #print(" lenght of leading phi series ",len(leadphi))
            #print(" lenght of leading eta series ",len(leadeta))
            #print(" lenght of hcandidare phi ",len(fdf['hCandidate_phi']))
            #print(" lenght of hcandidare eta ",len(fdf['hCandidate_eta']))
            #print(" lenght of leading lep dr h series ",len(deltaRllephdf))
            #print("printing the weird series")
            #print(deltaRllephdf)
            #print("the leading series")
            #print(leadphi)
            #print(leadeta)
            #print(fdf['hCandidate_phi'])
            #print(fdf['hCandidate_eta'])
            #print(fdf)

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
        #print(fdf['event_weight_btag'])
        eventweights = fdf['event_weight_kf']
        
        #print("btagging systematics debug")
        #print(" systl: ",systl)
        #print(" stype: ",stype)
        systname = 'btagsystdefaultname'
        if not systl and stype > 0:
            #print("if testing is working")
            print("    Applying nominal sf")
            systname = 'btagnom_muidnom'
            eventweights = fdf['event_weight_kf']*fdf['event_weight_btag']*fdf['event_weight_muid']
            #print(eventweights)
        elif not systl and stype <= 0:
            #print("if testing is working")
            print("    No systematic request and data, so no sf")
            systname = 'btagnom_muidnom'
            #print(eventweights)
        elif systl and stype <= 0:
            #print("if testing is working")
            print("    Attempt at a systematic, but dat, so no sf")
            systname = 'btagnom_muidnom'

        else:
            if "btagup" == systl:
                print("    Applying uncUp btagging SF")
                eventweights = fdf['event_weight_kf']*(fdf['event_weight_btag']+fdf['event_weight_btaguncup'])*fdf['event_weight_muid']
                systname = systl+'_muidnom'
            elif "btagdwn" == systl:
                print("    Applying uncDwn btagging SF")
                eventweights = fdf['event_weight_kf']*(fdf['event_weight_btag']-fdf['event_weight_btaguncdwn'])*fdf['event_weight_muid']
                systname = systl+'_muidnom'
            elif "muidup" == systl:
                print("    Applying uncUp MuonID SF")
                eventweights = fdf['event_weight_kf']*fdf['event_weight_btag']*(fdf['event_weight_muid']+fdf['event_weight_muiduncup'])
                systname = 'btagnom_'+systl
            elif "muiddwn" == systl:
                print("    Applying uncDwn MuonID SF")
                eventweights = fdf['event_weight_kf']*fdf['event_weight_btag']*(fdf['event_weight_muid']-fdf['event_weight_muiduncdwn'])
                systname = 'btagnom_'+systl
            else:
                print("    SF unc not recognized, using the nominal values")
                systname = systl
                eventweights = fdf['event_weight_kf']*fdf['event_weight_btag']*fdf['event_weight_muid']

        print("The event weights are wrong for general leading lepton dR plots")
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
        rootOutFile["h_met"]        = np.histogram(fdf['metsuable'],bins=39,range=(50,2000),weights=eventweights)
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
        meterrs      = boostUnc(fdf['metsuable'],eventweights,39,50,2000)
        metphierrs   = boostUnc(fdf['metphiusable'],eventweights,30,-3.14159,3.14159)
        metphiwerrs  = boostUnc(fdf['metphiusable'].map(wrapPhi),eventweights,30,0,3.14159)
        zpjigerrs    = boostUnc(fdf['ZPrime_mass_est'],eventweights,130,0,13000)
        ndjigerrs    = boostUnc(fdf['ND_mass_est'],eventweights,35,100,800)
        nsjigerrs    = boostUnc(fdf['NS_mass_est'],eventweights,25,0,500)
        btagerrs     = boostUnc(fdf['hCandidate_'+btaggr],eventweights,110,0,1.1)
        dphizherrs   = boostUnc(deltaphizhdf,eventweights,100,0,3.14159)
        dphizmeterrs   = boostUnc(deltaphizmetdf,eventweights,100,0,3.14159)
        dphihmeterrs   = boostUnc(deltaphihmetdf,eventweights,100,0,3.14159)
        drzherrs       = boostUnc(deltaRzhdf,eventweights,30,0,6)
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
                      zpjigerrs,
                      ndjigerrs,
                      nsjigerrs,
                      btagerrs,
                      dphizherrs,
                      dphizmeterrs,
                      dphihmeterrs,
                      drzherrs,
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
                     'h_zp_jigm',
                     'h_nd_jigm',
                     'h_ns_jigm',
                     'h_btag',
                     'h_dphi_zh',
                     'h_dphi_zmet',
                     'h_dphi_hmet',
                     'h_dr_zh',
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
            rootOutFile["h_LMu_pt"]     = np.histogram(fdf['LMuCandidate_pt'],bins=50,range=(0,500),weights=eventweights)
            rootOutFile["h_sLMu_pt"]     = np.histogram(fdf['sLMuCandidate_pt'],bins=50,range=(0,500),weights=eventweights)
            rootOutFile["h_LMu_phi"]    = np.histogram(fdf['LMuCandidate_phi'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_sLMu_phi"]    = np.histogram(fdf['sLMuCandidate_phi'],bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_LMu_eta"]    = np.histogram(fdf['LMuCandidate_eta'],bins=30,range=(-5,5),weights=eventweights)
            rootOutFile["h_sLMu_eta"]    = np.histogram(fdf['sLMuCandidate_eta'],bins=30,range=(-5,5),weights=eventweights)
            drlmuherrs     = boostUnc(deltaRlmuhdf,eventweights,30,0,6)
            drslmuherrs    = boostUnc(deltaRslmuhdf,eventweights,30,0,6)
            drslmulmuerrs  = boostUnc(deltaRslmulmudf,eventweights,30,0,6)
            lmupterrs      = boostUnc(fdf['LMuCandidate_pt'],eventweights,50,0,500)
            smupterrs      = boostUnc(fdf['sLMuCandidate_pt'],eventweights,50,0,500)
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
            rootOutFile["h_dr_leadleph"]  = np.histogram(deltaRllephdf,bins=30,range=(0,6),weights=eventweights)
            rootOutFile["h_dr_sleadleph"] = np.histogram(deltaRslephdf,bins=30,range=(0,6),weights=eventweights)
            rootOutFile["h_dr_leps"]      = np.histogram(deltaRsleplepdf,bins=30,range=(0,6),weights=eventweights)
            rootOutFile["h_leadlep_pt"]   = np.histogram(leadpt,bins=50,range=(0,500),weights=eventweights)
            rootOutFile["h_sleadlep_pt"]  = np.histogram(sleadpt,bins=50,range=(0,500),weights=eventweights)
            rootOutFile["h_leadlep_phi"]  = np.histogram(leadphi,bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_sleadlep_phi"] = np.histogram(sleadphi,bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_leadlep_eta"]  = np.histogram(leadeta,bins=30,range=(-5,5),weights=eventweights)
            rootOutFile["h_sleadlep_eta"] = np.histogram(sleadeta,bins=30,range=(-5,5),weights=eventweights)
            rootOutFile["h_electron_pt"]  = np.histogram(allelpt,bins=50,range=(0,500),weights=eventweights)
            rootOutFile["h_electron_phi"]  = np.histogram(allelphi,bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_electron_eta"]  = np.histogram(alleleta,bins=30,range=(-5,5),weights=eventweights)
            rootOutFile["h_muon_pt"]  = np.histogram(allmupt,bins=50,range=(0,500),weights=eventweights)
            rootOutFile["h_muon_phi"]  = np.histogram(allmuphi,bins=30,range=(-3.14159,3.14159),weights=eventweights)
            rootOutFile["h_muon_eta"]  = np.histogram(allmueta,bins=30,range=(-5,5),weights=eventweights)


            drllherrs   = boostUnc(deltaRllephdf,eventweights,30,0,6)
            drsllherrs  = boostUnc(deltaRslephdf,eventweights,30,0,6)
            drsllllerrs = boostUnc(deltaRsleplepdf,eventweights,30,0,6)
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
                         allmuetaerrs
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
                         "h_muon_eta"
                ]
            unc_arrays.extend(addhists)
            unc_names.extend(addnames)
            

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

                
