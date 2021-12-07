import uproot as up4
import uproot3 as up3
import pandas as pd
import numpy as np
import boost_histogram as bh
import argparse
import glob
import gecorg as go


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

if __name__=='__main__':
    parser.add_argument("-f","--sample",help = "sample file")
    parser.add_argument("-o","--output",help = "output file name")
    parser.add_argument("-b","--btagger",help = "btagger selection, need name part of key")
    parser.add_argument("-wp","--btagWP",type=float,default=0.6,help = "doulbeB tagger working point (default M1)")
    parser.add_argument("-zpt","--zPtCut",type=float,default = 100.0,help = "pT cut on Z")
    parser.add_argument("-hpt","--hPtCut",type=float,default = 250.0,help = "pT cut on h")
    parser.add_argument("-met","--metPtCut",type=float,default = 50.0,help = "pT cut on met")
    parser.add_argument("-sdm","--sdmCut",type=float,default = 10.0,help = "lowest soft drop mass cut")
    parser.add_argument("-date","--date",type=str,help = "where are your topiary plots?")
    parser.add_argument("-sr","--signalregion",type=bool,help = "do you want a signal region plot?")
    parser.add_argument("-c","--comboregion",type=bool,help = "do you want combined SR and SB?")
    parser.add_argument("-v","--validation",type=bool,help = "validation region bounds?")
    parser.add_argument("-syst","--systematics",type=str)
    parser.add_argument("-channel","--chan",type=str)
    args = parser.parse_args()

    samp   = args.sample
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

    inputfiles = glob.glob(args.date+'/'+samp+'*_topiary_'+channel+'*.root')
    #inputfiles = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/'+samp+'*_topiary*systnominal*.root')
   # inputfiles = glob.glob(args.date+'/'+samp+'*_topiary*systnominal*.root')

    for fjec in inputfiles:
        #for jet systematics
        samp = fjec.split("/")[-1]
        front = samp.split("_Zpt")[0]
        jectype = front.split("_")[-1]
        samp = front.split("_topiary")[0]
        #inputfiles = glob.glob('analysis_output_ZpAnomalon/'+args.date+'/'+samp+'*_topiary*'+jectype+'*.root')
        inputfiles = glob.glob(args.date+'/'+samp+'*_topiary_'+channel+'_'+jectype+'*.root')
    
        print("Doing selections on:")
        print("    ",inputfiles[:1])
        print("    ",samp)
        print("    Concerning systematics:")
        print("    ",jectype)
        stype,year = go.sampleType(samp)

        if not valid:
            if sr and stype != 0:
                print("    using signal region selections")
            elif comb and stype != 0:
                print("    using full region selections")
            else:
                print("    using sideband selections")

        if valid:
            print("    Doing validation of alpha method cuts")
            if sr:
                print("    'signalr' labeled events are in soft drop mass bands (55,70]")
            else:
                print("    'sideband' labeled events are in soft drop mass bands [30,55],[150,5000]")

        metstr = ''
        branches = [b'ZCandidate_*',
                    b'hCandidate_*',
                    #b'MET',
                    #b'METPhi',
                    b'metsuable',
                    b'metphiusable',
                    #b'METclean',
                    #b'METPhiclean',
                    b'ZPrime_mass_est',
                    b'ND_mass_est',
                    b'NS_mass_est',
                    b'event_weight*',
                    b'LMuCandidate_*',
                    b'sLMuCandidate_*',
                    #b'event_weight_kf',
                    #b'event_weight_btag'
        ]

    
        #events = up3.pandas.iterate(inputfiles[:1],'PreSelection;1',branches=branches)
        tree = up3.open(inputfiles[0])['PreSelection;1']
        events = tree.pandas.df(branches=branches)
        #events = tree.pandas.df()
        #print(events)
    
        # for b in events:
        #print(type(b))
        #print(events.keys())
        #print(b)
        #print("Doing SD mass lower cut of :" ,sdmcut)
        #do some cuts
        #print("Number of events in chunk ",len(events))
        sddf   = events[events['hCandidate_sd'] > sdmcut]
        metdf  = sddf[sddf['metsuable'] > metcut]
        zptdf  = metdf[metdf['ZCandidate_pt'] > zptcut]
        hptdf  = zptdf[zptdf['hCandidate_pt'] > hptcut]
        btdf   = hptdf[hptdf['hCandidate_'+btaggr] > float(btagwp)]

        #Actual Analysis
        if not valid:
            srup   = btdf[btdf['hCandidate_sd'] > 70.]
            bldf   = srup[srup['hCandidate_sd'] < 150.]#full blinded region
            srdf   = bldf[bldf['hCandidate_sd'] > 110.]#Higgs Peak
            lowsb  = btdf[btdf['hCandidate_sd'] <= 70.]
            highsb = btdf[btdf['hCandidate_sd'] >= 150.]
            sbdf   = pd.concat([lowsb,highsb])
            
        #Validation region for alpha method for DY
        if valid:
            srup   = btdf[btdf['hCandidate_sd'] > 55.]
            bldf   = srup[srup['hCandidate_sd'] < 150.]#full blinded region
            srdf   = bldf[bldf['hCandidate_sd'] <= 70.]#Validation region
            lowsb  = btdf[btdf['hCandidate_sd'] <= 55.]
            highsb = btdf[btdf['hCandidate_sd'] >= 150.]
            sbdf   = pd.concat([lowsb,highsb])


        region = "sideband"
        if not valid:
            if stype != 0:
                if sr:
                    fdf = srdf
                    region = "signalr"
                elif comb:
                    fdf = btdf
                    region = "totalr"
                else:
                    fdf = sbdf
                    region = "sideband"
            else:
                fdf = sbdf

        if valid:
            if sr:
                fdf = srdf
                region = "validationr"
            elif comb:
                fdf = btdf
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
        deltaRlmuhdf   = deltaR(fdf['LMuCandidate_phi'],fdf['hCandidate_phi'],fdf['LMuCandidate_eta'],fdf['hCandidate_eta'])
        deltaRslmuhdf   = deltaR(fdf['sLMuCandidate_phi'],fdf['hCandidate_phi'],fdf['sLMuCandidate_eta'],fdf['hCandidate_eta'])
        deltaRslmulmudf   = deltaR(fdf['sLMuCandidate_phi'],fdf['LMuCandidate_phi'],fdf['sLMuCandidate_eta'],fdf['LMuCandidate_eta'])

        #calculate the event weight
        #print(fdf['event_weight_btag'])
        eventweights = fdf['event_weight_kf']
        #print(eventweights)

        #print("btagging systematics debug")
        #print(" systl: ",systl)
        #print(" stype: ",stype)
        systname = 'btagsystdefaultname'
        if not systl and stype != 0:
            #print("if testing is working")
            print("    Applying nominal btagging sf")
            systname = 'btagnom'
            eventweights = fdf['event_weight_kf']*fdf['event_weight_btag']
            #print(eventweights)
        elif not systl and stype == 0:
            #print("if testing is working")
            print("    No btag syst but data, so no btag sf")
            systname = 'btagnom'
            #print(eventweights)
        elif systl and stype == 0:
            #print("if testing is working")
            print("    No btag syst but data, so no btag sf")
            systname = 'btagnom'

        else:
            if "btagup" == systl:
                print("    Applying uncUp btagging SF")
                eventweights = fdf['event_weight_kf']*(fdf['event_weight_btag']+fdf['event_weight_btaguncup'])
                systname = systl
            elif "btagdwn" == systl:
                print("    Applying uncDwn btagging SF")
                eventweights = fdf['event_weight_kf']*(fdf['event_weight_btag']-fdf['event_weight_btaguncdwn'])
                systname = systl
            else:
                print("    Btag unc not recognized, using the nominal values")
                systname = systl
                eventweights = fdf['event_weight_kf']*fdf['event_weight_btag']

        print("    number of passing events straight ",len(fdf))
        print("    number of passing events weighted ",eventweights.sum())
        print("    max Zp mass estimator: ",fdf['ZPrime_mass_est'].max())

        #lets make some histograms.
        rootfilename  = go.makeOutFile(samp,'upout_'+region+'_'+jectype+'_'+systname+'_'+btaggr,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))#need to update for btagger
        npfilename    = go.makeOutFile(samp,'totalevents_'+region+'_'+jectype+'_'+systname+'_'+btaggr,'.npy',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        pklfilename   = go.makeOutFile(samp,'selected_errors_'+region+'_'+jectype+'_'+systname+'_'+btaggr,'.pkl',str(zptcut),str(hptcut),str(metcut),str(btagwp))
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
        rootOutFile["h_dr_lmuh"]      = np.histogram(deltaRlmuhdf,bins=30,range=(0,6),weights=eventweights)
        rootOutFile["h_dr_slmuh"]      = np.histogram(deltaRslmuhdf,bins=30,range=(0,6),weights=eventweights)
        rootOutFile["h_dr_slmulmu"]      = np.histogram(deltaRslmulmudf,bins=30,range=(0,6),weights=eventweights)

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
        drlmuherrs     = boostUnc(deltaRlmuhdf,eventweights,30,0,6)
        drslmuherrs    = boostUnc(deltaRslmuhdf,eventweights,30,0,6)
        drslmulmuerrs    = boostUnc(deltaRslmulmudf,eventweights,30,0,6)
        
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
                      drlmuherrs,
                      drslmuherrs,
                      drslmulmuerrs,
                      
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
                     'h_dr_lmuh',
                     'h_dr_slmuh',
                     'h_dr_slmulmu',
        ]

        max_length = len(max(unc_arrays,key = lambda ar : len(ar)))
        pad_arrays = [np.pad(arr,(0,max_length - len(arr)),'constant') for arr in unc_arrays]
        all_unc    = np.column_stack(pad_arrays)
        uncdf      = pd.DataFrame(all_unc,columns=unc_names)
        uncdf.to_pickle("./"+pklfilename)
    
        #Book Keeping
        f = up3.open(inputfiles[0])
        np.save(npOutFile,np.array([f['hnorigevnts'].values[0]]))
        rootOutFile["hnevents"]      = str(f['hnorigevnts'].values[0])
        rootOutFile["hnevents_pMET"] = str((metdf['event_weight_kf']*metdf['event_weight_btag']).sum())#str(len(metdf))
        rootOutFile["hnevents_pZ"]   = str((zptdf['event_weight_kf']*zptdf['event_weight_btag']).sum())#str(len(zptdf))
        rootOutFile["hnevents_ph"]   = str((hptdf['event_weight_kf']*hptdf['event_weight_btag']).sum())#str(len(hptdf))
        rootOutFile["hnevents_sb"]   = str((sbdf['event_weight_kf']*sbdf['event_weight_btag']).sum())#str(len(sbdf))
        rootOutFile["hnevents_btag"] = str((btdf['event_weight_kf']*btdf['event_weight_btag']).sum())#str(len(btdf))

        #print("Passing events sb events, not sf              : ",len(sbdf))
        #print("Passing events sb events, sum of event weights: ",
        
        if stype != 0:
            if sr:
                rootOutFile["hnevents_sr"]   = str((srdf['event_weight_kf']*srdf['event_weight_btag']).sum())#str(len(srdf))
