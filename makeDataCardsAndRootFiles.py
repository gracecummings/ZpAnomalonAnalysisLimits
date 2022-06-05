import ROOT
import glob
import os
import sys
import gecorg_test as go
import numpy as np
import pandas as pd
import configparser
import argparse

def doStatsUncertainty(hist):
    bins = hist.GetNbinsX()
    name = hist.GetName()
    histlist = []
    for b in range(bins+1):
        #print("The content of bin {0} is {1} with an error of {2}".format(b,hist.GetBinContent(b),hist.GetBinError(b)))
        if b == 0:
            continue
        hname = name+"_StatsUncBin"+str(b)
        bincont = hist.GetBinContent(b)
        binunc  = hist.GetBinError(b)
        hup = hist.Clone()
        hdwn = hist.Clone()
        hup.SetBinContent(b,bincont+binunc)
        hdwn.SetBinContent(b,bincont-binunc)
        hup.SetName(hname+"Up")
        hdwn.SetName(hname+"Down")
        histlist.append(hup)
        histlist.append(hdwn)
    return histlist
        
        

def makeSignalInfoDict(sigclass, region,sigxs):
    sigs = sigclass.getPreppedSig(region,sigxs)
    sigdict = {}
    for sig in sigs:
        sigdict[sig["name"]] = sig
    return sigdict

def checkBinNumbers(histlist):
    binnums = [hist.GetNbinsX() for hist in histlist]
    testl = [num == binnums[0] for num in binnums]
    aok = all(testl)
    return aok

def writeDataCard(processes,rootFileName,channel,yearstr,hdat):
    signame = processes["processnames"][0]
    nbkg = len(processes["processnames"])-1
    obsEvents = hdat.Integral()
    binname = signame+"_"+channel
    namestr = " ".join(processes["processnames"])
    rates   = [str(hist.Integral()) for hist in processes["hists"]]
    prepCardName = go.makeOutFile('Run2_'+yearstr+'_ZllHbbMET','datacard_'+chan+'_'+signame,'.txt',str(zptcut),str(hptcut),str(metcut),str(btagwp))
    card = open(prepCardName,"w")

    #Write the card
    card.write("imax 1\n")
    card.write("jmax {0}\n".format(nbkg))
    card.write("kmax *\n")              
    card.write("------------\n")
    card.write("shapes * * {0} $PROCESS $PROCESS_$SYSTEMATIC \n".format(rootFileName.split("/")[-1]))
    card.write("------------\n")
    card.write("bin {0} \n".format(binname))
    card.write("observation {0} \n".format(obsEvents))
    card.write("------------\n")
    card.write("bin {0} {1} {2} {3}\n".format(binname,binname,binname,binname))
    card.write("process "+namestr+"\n")
    card.write("process 0 1 2 3\n")#hardcode
    card.write("rate "+" ".join(rates)+"\n")
    card.write("------------\n")

    maxnamelength = max([len(x) for x in processes["systrates"].keys()])

    for syst in processes["systrates"].keys():
        vals = processes["systrates"][syst]["proc"]
        sval = [x if len(x) > 1 else "-" for x in vals]

        nameoffset = maxnamelength - len(syst)
        
        cardstr = "{0} {1} {2} {3} {4} {5}\n".format(str(syst)+" "*nameoffset,processes["systrates"][syst]["type"],sval[0],sval[1],sval[2],sval[3])
        card.write(cardstr)
    
    for syst in processes["systshapes"].keys():
        vals = processes["systshapes"][syst]["proc"]
        sval = [x if x != 0 else "-" for x in vals]
        nameoffset = maxnamelength - len(syst)
        cardstr = "{0} {1} {2} {3} {4} {5}\n".format(str(syst)+" "*nameoffset,processes["systshapes"][syst]["type"],sval[0],sval[1],sval[2],sval[3])
        card.write(cardstr)

    card.write("* autoMCStats 1\n")
    card.close()

#def gatherSystematics(confsecsyst):

def newNameAndStructure(hist,name,rebindiv,limrangelow,limrangehigh):
    hist.Rebin(rebindiv)
    nbins = hist.GetNbinsX()
    binw  = hist.GetBinWidth(1)
    newbins = [limrangelow+x*binw for x in range(int((limrangehigh-limrangelow)/binw))]
    nh = ROOT.TH1F(name,name,len(newbins),limrangelow,limrangehigh)
    for b,le in enumerate(newbins):
        bnum = hist.FindBin(le)
        bincontent = hist.GetBinContent(bnum)
        binerror   = hist.GetBinError(bnum)
        nh.SetBinContent(b+1,bincontent)
        nh.SetBinError(b+1,binerror)

    return nh
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m","--metcut", type=float,help = "met cut of samples")
    parser.add_argument("-z","--zptcut", type=float,help = "zpt cut of samples")
    parser.add_argument("-j","--hptcut", type=float,help = "hpt cut of samples")
    parser.add_argument("-wp","--btagwp", type=float,help = "btag working point")
    args = parser.parse_args()
    
    zptcut  = str(args.zptcut)#'150.0'
    hptcut  = str(args.hptcut)#'300.0'
    metcut  = str(args.metcut)#'200.0'
    btagwp  = str(args.btagwp)#'0.8'

    chan    = 'mumu'
    years   = [16,17,18]
    yearstr = go.yearFormatter(years)
    sigxs   = 1.0
    rebindiv = 2

    limrangelow = 1400
    limrangehigh = 3000

    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()

    ####load in the files with the nominal distributions
    bkgs = go.backgrounds(config.get('nominal','pathnom'),zptcut,hptcut,metcut,btagwp,config.get('nominal','strnom'))
    sig  = go.signal(config.get('nominal','pathsignom'),zptcut,hptcut,metcut,btagwp,sigxs,[16,17,18],config.get('nominal','strnom'))
    dyEst = ROOT.TFile(config.get('nominal','pathnom')+'/Run2_161718_dy_extraploation'+config.get('nominal','strnom')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')


    ####Prepping holders####
    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty1 = empty.Clone()
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    
    ####Getting the Estimations####
    hdy    = dyEst.Get("extrphistnoerrs").Clone()
    htt = bkgs.getAddedHist(empty1,"TT","sr","h_zp_jigm")
    hzz  = bkgs.getAddedHist(empty2,"ZZTo2L2Q","sr","h_zp_jigm")
    hwz  = bkgs.getAddedHist(empty3,"WZTo2L2Q","sr","h_zp_jigm")
    hvv  = hzz.Clone()
    hvv.Add(hwz)


    ####Rename and restucture
    htt = newNameAndStructure(htt,"TT",rebindiv,limrangelow,limrangehigh)
    hdy = newNameAndStructure(hdy,"DY",1,limrangelow,limrangehigh)
    hvv = newNameAndStructure(hvv,"VV",rebindiv,limrangelow,limrangehigh)
    hdat = hdy.Clone()
    hdat.SetName("data_obs")
    hdat.Add(htt)
    hdat.Add(hvv)

    ####Do bkg stats unc explicitly
    httstatsunc = doStatsUncertainty(htt)

    ####For Each signal, make a datacard, and a root file with all systematics
    siginfo = sig.getPreppedSig('sr',sigxs)
    signom = makeSignalInfoDict(sig,'sr',sigxs)
    sigcolors = go.colsFromPalette(siginfo,ROOT.kCMYK)
    nomsigs = [s["name"] for s in siginfo]
    for s in nomsigs:
        name = signom[s]["name"]
        signame = "holder"
        if "Tune" in name:
            strippedname = name.split("_Tune")[0]
            signame = strippedname.replace("-","")
        else:
            signame = name.replace("-","")
        print("------- Looking at signal sample ",signame)

        if "Zp4000ND800NS200" not in signame:
            continue

        ####Make Files
        prepRootName = go.makeOutFile('Run2_'+yearstr+'_ZllHbbMET',chan+'_'+signame,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        prepRootFile = ROOT.TFile(prepRootName,"recreate")

        ###Write bkg stats files
        for hist in httstatsunc:
            hist.Write()
            
        ####Nominal Signal
        hsigori = signom[s]["tfile"].Get("h_zp_jigm")
        hsigori.Sumw2(ROOT.kTRUE)#Throws a warning that it is already created
        hsig = hsigori.Clone()
        hsig.Scale(signom[s]["scale"])

        hsig = newNameAndStructure(hsig,signame,rebindiv,limrangelow,limrangehigh)
        prepRootFile.cd()

        #print("!!!!!!!!artificially scaling the nominal backgrounds!!!!")
        #bkghists = [htt,hvv,hdy,hdat]
        #for hist in bkghists:
        #    hist.Scale(10)
        htt.Write()
        hvv.Write()
        hdy.Write()
        hdat.Write()
        hsig.Write()

        #####Gather Systematics
        systdictrate = {"lumi_13TeV":{"type":"lnN","unc":1.018,"proc":["1.018","1.018","1.018","1.018"]}
                    }
        systdictshape = {}


        for syst in systs:
            #print("YOU HAVE CURRENTLY TURNED OFF SYSTEMATICS")
            #continue
            if syst == 'nominal':
                continue
            print("        Looking at systematic ",syst)

            
            systbkgsup  = go.backgrounds(config.get(syst,'pathup'),zptcut,hptcut,metcut,btagwp,config.get(syst,'strup'))
            systbkgsdwn = go.backgrounds(config.get(syst,'pathdwn'),zptcut,hptcut,metcut,btagwp,config.get(syst,'strdwn'))
            systsigup   = go.signal(config.get(syst,'pathsigup'),zptcut,hptcut,metcut,btagwp,sigxs,[16,17,18],config.get(syst,'strup'))
            systsigdwn  = go.signal(config.get(syst,'pathsigdwn'),zptcut,hptcut,metcut,btagwp,sigxs,[16,17,18],config.get(syst,'strdwn'))

            if len(systbkgsup.bkgs["DYJetsToLL"][18]["sb"][0]) < 1:
                print("        There are no DYJets entires for this systematic.")
                print("        moving on. This systematic will not included.")
                continue
            
            appcode = config.get(syst,'applist').split(',')
            ratenums = config.get(syst,'rate').split(',')
            ratelist = ratenums
            if len(appcode > 1):#If the shape systematic never gets applied
                applist = [float(x) for x in appcode]
                systdictshape[syst] = {"type":config.get(syst,'type'),"unc":1.0,"proc":applist}
                systdictrate[syst]  = {"type":config.get(syst,'typerate'),"proc":ratelist}
            else:
                systdictrate[syst]  = {"type":config.get(syst,'typerate'),"proc":ratelist}
            
            if rebindiv == 2:
                dyEstup     = ROOT.TFile(config.get(syst,'pathup')+'/Run2_161718_dy_extraploation'+config.get(syst,'strup')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
                dyEstdwn    = ROOT.TFile(config.get(syst,'pathdwn')+'/Run2_161718_dy_extraploation'+config.get(syst,'strdwn')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')

            ####Make it useful
            #sigup  = systsigup.getPreppedSig('sr',sigxs)#systsigup.prepsigsr
            #sigupdict = {}
            sigup = makeSignalInfoDict(systsigup,'sr',sigxs)
            sigdwn = makeSignalInfoDict(systsigdwn,'sr',sigxs)
            keyup = sigup[s]["tfile"].GetListOfKeys()
            keydwn = sigdwn[s]["tfile"].GetListOfKeys()
            keyup = [k.GetName() for k in keyup]
            keydwn = [k.GetName() for k in keydwn]

            #for k,key in enumerate(keyup):
            #    if "hnevents" in key:
            #        break
            #    print("Up key:   ",key)
            #    print("Down key: ",keydwn[k])
            #    print(" Up bins: ",sigup[s]["tfile"].Get(key).GetNbinsX())
            #    print("Dwn bins: ",sigdwn[s]["tfile"].Get(key).GetNbinsX())

            ####Check is the keys exist
            if (s in sigup) and (s in sigdwn): 
                ####Prepping holders####
                empty4 = empty.Clone()
                empty5 = empty.Clone()
                empty6 = empty.Clone()
                empty7 = empty.Clone()
                empty8 = empty.Clone()
                empty9 = empty.Clone()

                ####Gathering the Systematic
                #Background
                httup  = systbkgsup.getAddedHist(empty4,"TT","sr","h_zp_jigm")
                httdwn = systbkgsdwn.getAddedHist(empty5,"TT","sr","h_zp_jigm")
                hzzup  = systbkgsup.getAddedHist(empty6,"ZZTo2L2Q","sr","h_zp_jigm")
                hzzdwn = systbkgsdwn.getAddedHist(empty7,"ZZTo2L2Q","sr","h_zp_jigm")
                hwzup  = systbkgsup.getAddedHist(empty8,"WZTo2L2Q","sr","h_zp_jigm")
                hwzdwn = systbkgsdwn.getAddedHist(empty9,"WZTo2L2Q","sr","h_zp_jigm")
                hvvup  = hzzup.Clone()
                hvvdwn = hzzdwn.Clone()
                hvvup.Add(hwzup)
                hvvdwn.Add(hwzdwn)
                hdyup  = dyEstup.Get("extrphistnoerrs").Clone()
                hdydwn = dyEstdwn.Get("extrphistnoerrs").Clone()

                #Signal
                hsigupori = sigup[s]["tfile"].Get("h_zp_jigm")
                #print("              hup nbins  ",hsigupori.GetNbinsX())
                hsigup = hsigupori.Clone()
                hsigup.Scale(sigup[s]["scale"])
                hsigdwnori = sigdwn[s]["tfile"].Get("h_zp_jigm")
                #print("              hdwn nbins ",hsigdwnori.GetNbinsX())

                #if hsigupori.GetNbinsX() != hsigdwnori.GetNbinsX():
                #    print("___________________________________________________--bins weird")
                #    print("___________________________________________________-- ",sigdwn[s]["tfile"])

                hsigdwn = hsigdwnori.Clone()
                hsigdwn.Scale(sigdwn[s]["scale"])
                #print(sigdwn[s]["tfile"])
                #print(hsigdwnori)
                #print(hsigdwnori.Integral())
                #print(hsigdwn.Integral())

                #Rename and Restructure
                hsigdwn = newNameAndStructure(hsigdwn,signame+"_"+syst+"Down",rebindiv,limrangelow,limrangehigh)
                hsigup = newNameAndStructure(hsigup,signame+"_"+syst+"Up",rebindiv,limrangelow,limrangehigh)
                httup = newNameAndStructure(httup,"TT_"+syst+"Up",rebindiv,limrangelow,limrangehigh)
                httdwn = newNameAndStructure(httdwn,"TT_"+syst+"Down",rebindiv,limrangelow,limrangehigh)
                hvvup = newNameAndStructure(hvvup,"VV_"+syst+"Up",rebindiv,limrangelow,limrangehigh) 
                hvvdwn = newNameAndStructure(hvvdwn,"VV_"+syst+"Down",rebindiv,limrangelow,limrangehigh)
                hdyup  = newNameAndStructure(hdyup,"DY_"+syst+"Up",1,limrangelow,limrangehigh)
                hdydwn = newNameAndStructure(hdydwn,"DY_"+syst+"Down",1,limrangelow,limrangehigh)

                allsysthists = [hsigdwn,hsigup,httup,httdwn,hvvup,hvvdwn,hdyup,hdydwn]
                allthesame = checkBinNumbers(allsysthists)
                if not allthesame:
                    print("              NOT ALL OF THE BINS NUMBERS ARE THE SAME")
                    print("              NOT WRITING THE {0} UP/DWN HISTOGRAMS".format(syst))
                    continue


                #Write the histograms
                
                prepRootFile.cd()
                print("        Writing the up/dwn histograms")
                #print("!!!!!!!!artificially scaling the up/dwn backgrounds!!!!")
                #bkghists = [httup,hvvup,hdyup,httdwn,hvvdwn,hdydwn]
                #for hist in bkghists:
                #    hist.Scale(10)
                hsigup.Write()
                hsigdwn.Write()
                httup.Write()
                httdwn.Write()
                hvvup.Write()
                hvvdwn.Write()
                hdyup.Write()
                hdydwn.Write()
            else:
                print("****************{0} Does not have {1} systematic".format(s,syst))

        #For writing the datacard
        procdict = {"processnames":[signame,"DY","TT","VV"],
                    "hists":[hsig,hdy,htt,hvv],
                    "method":["mc","alpha","mc","mc"],
                    "systshapes":systdictshape,
                    "systrates":systdictrate,
        }
        print("        Defined the Datacard Dict")
        print("        Writing the Datacard")
        writeDataCard(procdict,prepRootName,chan,yearstr,hdat)
        print("        About to close the root file")
        prepRootFile.Close()

