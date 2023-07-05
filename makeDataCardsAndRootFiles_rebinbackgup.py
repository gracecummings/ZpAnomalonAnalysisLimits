import ROOT
import glob
import os
import sys
import gecorg_test as go
import numpy as np
import pandas as pd
import configparser
import argparse

def getDeviatedOverNominal(hist,histnom):
    interest = hist.Integral()
    intnom = histnom.Integral()
    devonom = interest/intnom
    return devonom

def applyStatsUncToSignal(hist,errseries,scale):
    alpha = 1.- 0.682689492
    for ibin in range(hist.GetNbinsX()+1):
        if ibin == 0:
            continue
        elif hist.GetBinContent(ibin) > 0:
            binerr = errseries[ibin-1]
            hist.SetBinError(ibin,binerr)
        else:
            binerrbase = ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin))+1,1)-hist.GetBinContent(ibin)
            binerr = binerrbase*scale
            hist.SetBinError(ibin,binerr)

            
    return hist

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

def gatherEmuStats(stattf,desiredbins,limrangelow,limrangehigh):
    #this function has to
    #1.) Get the stats up/downs from a file
    #2.) Only get the up/downs for shifted bins we care about
    #3.) rename and restructure those histograms to have the correct bin number
    #4.) make the dictionary for the datacard writing
    applist = [0,0,1.0,0]
    statsdict = {}
    huncs = []
    keys = stattf.GetListOfKeys()
    keys = [key.GetName() for key in keys]
    keys = sorted(keys,key=lambda k:int(k.split('StatsUncBin')[-1]))
    for key in keys:
        oribin = key.split('StatsUncBin')[-1]
        if int(oribin) in desiredbins:#only taking hists for bins we care about
            newbinnum = int(oribin)-7
            h = stattf.Get(key)#get histogram of interest
            direc = key.split('_TT_TT')[0]
            newname = "TT_TT_StatsUncBin"+str(newbinnum)+direc#new hist name
            trimmedh = newNameAndStructure(h,newname,1,limrangelow,limrangehigh)#only take bins we care about in hist
            #marking the hists to be saved
            huncs.append(trimmedh)
            uncdict  = {"TT_StatsUncBin"+str(newbinnum):{"type":"shape","unc":1.0,"proc":applist}}
            statsdict.update(uncdict)
        else:
            continue

    return huncs,statsdict

def doStatsUncertainty(hist):
    bins = hist.GetNbinsX()
    name = hist.GetName()
    appdict = {"TT":[0,0,1.0,0],"DY":[0,1.0,0,0],"VV":[0,0,0,1.0],"sig":[1.0,0,0,0]}
    napp = name
    if "Zp" in name:
        napp = "sig"
    applist = appdict[napp]
    histlist = []
    statsdict = {}
    for b in range(bins+1):
        #print("The content of bin {0} is {1} with an error of {2}".format(b,hist.GetBinContent(b),hist.GetBinError(b)))
        if b == 0:
            continue
        hname = name+"_"+name+"_StatsUncBin"+str(b)
        bincont = hist.GetBinContent(b)
        binunc  = hist.GetBinError(b)
        hup = hist.Clone()
        hdwn = hist.Clone()
        if bincont > 0:
            hup.SetBinContent(b,bincont+binunc)
            hdwn.SetBinContent(b,bincont-binunc)
            if bincont-binunc < 0:
                print("+++++++++++++++++++++++++++++++++++++++++++++++++++++negative input+++")
                print(hname)
                print(b)
                print(bincont)
                print(binunc)
                hdwn.SetBinContent(b,bincont-binunc)
        #else:
        #    c     = hist.GetBinContent(b)
        #    alpha = 1.- 0.682689492
        #    binuncup = ROOT.Math.gamma_quantile_c(alpha/2,int(c)+1,1)-c
        #    hup.SetBinContent(b,bincont+binuncup)
        else:
            hup.SetBinContent(b,bincont+binunc)#zero in bins, garwood interval scaled should be err

        hup.SetName(hname+"Up")
        hdwn.SetName(hname+"Down")
        histlist.append(hup)
        histlist.append(hdwn)
        uncdict  = {name+"_StatsUncBin"+str(b):{"type":"shape","unc":1.0,"proc":applist}}
        statsdict.update(uncdict)
        #print(updict)
        #print(dwndict)
        #print(statsdict)
    return histlist,statsdict
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

def writeDataCard(processes,rootFileName,channel,yearstr,hdat,noratel,noshapel):
    signame = processes["processnames"][0]
    nbkg = len(processes["processnames"])-1
    obsEvents = hdat.Integral()
    binname = signame+"_"+channel
    namestr = " ".join(processes["processnames"])
    rates   = [str(hist.Integral()) for hist in processes["hists"]]
    prepCardName = go.makeOutFile('Run2_'+yearstr+'_ZllHbbMET','datacard_extended_poisson_'+chan+'_'+signame,'.txt',str(zptcut),str(hptcut),str(metcut),str(btagwp))
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
        if syst in noratel:
            systname = "#"+str(syst)
        else:
            systname = str(syst)
        vals = processes["systrates"][syst]["proc"]
        sval = [x if len(x) > 1 else "-" for x in vals]

        nameoffset = maxnamelength - len(syst)
        
        cardstr = "{0} {1} {2} {3} {4} {5}\n".format(systname+" "*nameoffset,processes["systrates"][syst]["type"],sval[0],sval[1],sval[2],sval[3])
        card.write(cardstr)
    
    for syst in processes["systshapes"].keys():
        if syst in noshapel:
            systname = "#"+str(syst)
        else:
            systname = str(syst)

        vals = processes["systshapes"][syst]["proc"]
        sval = [x if x != 0 else "-" for x in vals]
        nameoffset = maxnamelength - len(syst)
        cardstr = "{0} {1} {2} {3} {4} {5}\n".format(systname+" "*nameoffset,processes["systshapes"][syst]["type"],sval[0],sval[1],sval[2],sval[3])
        card.write(cardstr)

    #card.write("* autoMCStats 1\n")
    card.close()

#def gatherSystematics(confsecsyst):

def getBinNumbersOfInterest(hist,rebindiv,limrangelow,limrangehigh):
    hist.Rebin(rebindiv)
    nbins = hist.GetNbinsX()
    binw  = hist.GetBinWidth(1)
    desiedges = [limrangelow+x*binw for x in range(int((limrangehigh-limrangelow)/binw))]
    desibins  = [edge/binw+1 for edge in desiedges]
    setdesi = set(desibins)
    #print("Bin numbers of interest ",desibins)
    #print("chec against hist finder:")
    #checked = []
    #for b,le in enumerate(desiedges):
    #    bnum = hist.FindBin(le)
    #    checked.append(bnum)
    #    print(bnum)
    #setcheck = set(checked)
    #inter = setcheck & setcheck
    #diff = inter-setcheck
    #print("Difference in the intersection of sets: ",diff)

    return setdesi

def gatherAlphaMethodUncs(dytf,nomdy,limrangelow,limrangehigh):
    rateuncdict = {}
    shapeuncdict = {}
    histlist = []

    keys = dytf.GetListOfKeys()
    keys = [key.GetName() for key in keys]
    unckeys = [key for key in keys if "extrap_" in key]
    unckeysup = [key for key in unckeys if "Up" in key]
    unckeysdn = [key for key in unckeys if "Down" in key]
    unckeysup = sorted(unckeysup)
    unckeysdn = sorted(unckeysdn)

    for i in range(len(unckeysup)):
        upkey = unckeysup[i]
        dnkey = unckeysdn[i]
        upr = dytf.Get(upkey)
        dnr = dytf.Get(dnkey)
        up = newNameAndStructure(upr,"DY_"+upkey,1,limrangelow,limrangehigh)
        dn = newNameAndStructure(dnr,"DY_"+dnkey,1,limrangelow,limrangehigh)

        #gather this histograms
        histlist.append(up)
        histlist.append(dn)

        #make the rate uncertainties
        uprate = getDeviatedOverNominal(up,nomdy)
        dwnrate = getDeviatedOverNominal(dn,nomdy)
        genname = upkey.split("Up")[0]
        ratestr = "-"
        if abs(round(1-uprate,2))-abs(round(1-dwnrate,2)) == 0:
            ratestr = str(max(round(uprate,2),round(dwnrate,2)))
        #elif "vv" not in genname:
        #     ratestr = str(round(dwnrate,2))+"/"+str(round(uprate,2))
        else:
            ratestr = str(round(dwnrate,3))+"/"+str(round(uprate,3))
        rateuncdict[genname] = {"type":"lnN","proc":["-",ratestr,"-","-"]}

        #make the shape uncertainties (can probably omit some later)
        shapeuncdict[genname] = {"type":"shape","unc":1.0,"proc":[0,1,0,0]}

    return histlist,shapeuncdict,rateuncdict

def gatherEmuUncs(procname,tf,hnom,limrangelow,limrangehigh):
    rateuncdict = {}
    shapeuncdict = {}
    histlist = []
    
    keys = tf.GetListOfKeys()
    keys = [key.GetName() for key in keys]
    unckeys = [key for key in keys if "nom" not in key]
    unckeysup = [key for key in unckeys if "Up" in key]
    unckeysdn = [key for key in unckeys if "Down" in key]
    unckeysup = sorted(unckeysup)
    unckeysdn = sorted(unckeysdn)

    for i in range(len(unckeysup)):
        upkey = unckeysup[i]
        dnkey = unckeysdn[i]
        upr = tf.Get(upkey)
        dnr = tf.Get(dnkey)
        up = newNameAndStructure(upr,procname+"_"+upkey,1,limrangelow,limrangehigh)
        dn = newNameAndStructure(dnr,procname+"_"+dnkey,1,limrangelow,limrangehigh)

        #gather this histograms
        histlist.append(up)
        histlist.append(dn)

        #make the rate uncertainties
        uprate = getDeviatedOverNominal(up,hnom)
        dwnrate = getDeviatedOverNominal(dn,hnom)
        genname = upkey.split("Up")[0]
        ratestr = "-"
        if abs(round(1-uprate,2))-abs(round(1-dwnrate,2)) == 0:
            ratestr = str(max(round(uprate,2),round(dwnrate,2)))
        #elif "vv" not in genname:
        #     ratestr = str(round(dwnrate,2))+"/"+str(round(uprate,2))
        else:
            ratestr = str(round(dwnrate,3))+"/"+str(round(uprate,3))
        rateuncdict[genname] = {"type":"lnN","proc":["-","-",ratestr,"-"]}

        #make the shape uncertainties (can probably omit some later)
        shapeuncdict[genname] = {"type":"shape","unc":1.0,"proc":[0,0,1,0]}

    return histlist,shapeuncdict,rateuncdict


def makeDeviatedTPave(histup,hist,histdwn,name):
    updiv  = getDeviatedOverNominal(histup,hist)
    dwndiv = getDeviatedOverNominal(histdwn,hist)
    updiv = "up deviated over nominal: "+str(round(updiv,3))
    dwndiv = "dwn deviated over nominal: "+str(round(dwndiv,3))
    lab = ROOT.TPaveText(.2,.5,.55,.65,"NBNDC")
    lab.AddText(updiv)
    lab.AddText(dwndiv)
    lab.SetFillColor(0)
    return lab


def plotSignalSystematics(hup,hdwn,hnom,name):
    tc = ROOT.TCanvas("tc","tc",700,600)
    p1 = ROOT.TPad("p1","plot",0,0.3,1.0,1.0)
    p2 = ROOT.TPad("p2","ratio",0,0,1.0,0.3)
    leg = ROOT.TLegend(0.15,0.7,0.55,0.85)
    lab = makeDeviatedTPave(hup,hnom,hdwn,name)

    #drawing style
    #ROOT.gStyle.SetErrorX(0)
    plotmax = max([hup.GetMaximum(),hdwn.GetMaximum(),hnom.GetMaximum()])*1.3
    hists = [hup,hdwn,hnom]
    for h in hists:
        h.GetYaxis().SetRangeUser(0,plotmax)
    hup.SetLineColor(ROOT.kRed)
    hdwn.SetLineColor(ROOT.kBlue)
    hnom.SetLineColor(ROOT.kBlack)
    leg.SetBorderSize(0)
    leg.AddEntry(hup,name+" up","l")
    leg.AddEntry(hnom,name+" nominal","l")
    leg.AddEntry(hdwn,name+" down","l")
    hdiv = hup.Clone()
    hdiv.Divide(hup,hdwn)
    hdiv.SetLineColor(ROOT.kBlack)
    hdiv.SetMarkerStyle(8)
    hdiv.SetMarkerSize(0.7)
    hdiv.SetMarkerColor(ROOT.kBlack)
    hdiv.SetStats(0)
    hdiv.GetXaxis().SetLabelSize(0.075)
    hdiv.GetYaxis().SetTitle("up/dwn")
    hdiv.GetYaxis().SetTitleSize(0.11)
    hdiv.GetYaxis().SetTitleOffset(.45)
    hdiv.GetYaxis().SetLabelSize(0.08)
    hdiv.GetYaxis().SetLabelOffset(0.02)
    hdiv.GetYaxis().SetNdivisions(503)
    

    #Draw
    tc.Draw()
    tc.cd()
    p1.Draw()
    p1.cd()
    hnom.Draw("hist")
    hup.Draw("histsame")
    hdwn.Draw("histsame")
    leg.Draw()
    lab.Draw()
    tc.cd()
    p2.Draw()
    p2.cd()
    hdiv.Draw("e1")
    tc.Update()
    tc.SaveAs(go.makeOutFile('Run2_161718_ZllHbbMET_sigalSystematics',name,'.png',str(zptcut),str(hptcut),str(metcut),str(btagwp)))
    

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
    limrangehigh = 10000

    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()

    ####load in the files with the nominal distributions
    bkgs = go.backgrounds(config.get('nominal','pathnom'),zptcut,hptcut,metcut,btagwp,config.get('nominal','strnom'))
    sig  = go.signal(config.get('nominal','pathsignom'),zptcut,hptcut,metcut,btagwp,sigxs,[16,17,18],config.get('nominal','strnom'))
    dyEst = ROOT.TFile(config.get('nominal','pathnom')+'/Run2_161718_dy_extraploationalphat_'+config.get('nominal','strnom')+'__Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')

    #print("hardcoded the DY files!!!!!!")
    #dyEst = ROOT.TFile(config.get('nominal','pathnom')+'/Run2_161718_dy_extraploationsystnominal_kfnom_btagnom_muidnom_elidnom_elreconom__Zptcut100.0_Hptcut300.0_metcut0.0_btagwp0.8.root')
    #dyEst = ROOT.TFile(config.get('nominal','pathnom')+'/backup_before_flip_Run2_161718_dy_extraploationsystnominal_kfnom_btagnom_muidnom_elidnom_elreconom_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root')
    
    ttEst = ROOT.TFile(config.get('nominal','pathemu')+'/Run2_161718_emu_extrapolation_'+config.get('nominal','stremu')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
    ttStats = ROOT.TFile(config.get('nominal','pathemu')+'/Run2_161718_emu_extrapStats_'+config.get('nominal','stremu')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')

    ####Prepping holders####
    tf1 = ROOT.TFile(bkgs.f17dyjetsb[0])
    empty = tf1.Get('h_zp_jigm')
    empty.Reset("ICESM")#creates an empty hist with same structure
    empty1 = empty.Clone()
    empty2 = empty.Clone()
    empty3 = empty.Clone()
    
    ####Getting the Estimations####
    hdy    = dyEst.Get("extrphistnoerrs").Clone()
    htt    = ttEst.Get("TT_scalenominal").Clone()
    #htt = bkgs.getAddedHist(empty1,"TT","sr","h_zp_jigm")#mc way
    hzz  = bkgs.getAddedHist(empty2,"ZZTo2L2Q","sr","h_zp_jigm")
    hwz  = bkgs.getAddedHist(empty3,"WZTo2L2Q","sr","h_zp_jigm")
    hvv  = hzz.Clone()
    hvv.Add(hwz)

    ####Rename and restucture
    desibins = getBinNumbersOfInterest(empty.Clone(),2,limrangelow,limrangehigh)
    #htt = newNameAndStructure(htt,"TT",rebindiv,limrangelow,limrangehigh)#mc way
    htt = newNameAndStructure(htt,"TT",1,limrangelow,limrangehigh)
    hdy = newNameAndStructure(hdy,"DY",1,limrangelow,limrangehigh)
    hvv = newNameAndStructure(hvv,"VV",rebindiv,limrangelow,limrangehigh)

    ####Make data obs hist
    #tfd = ROOT.TFile(config.get('nominal','pathnom')+"/Run2_161718_gofplots_"+config.get('nominal','strnom')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
    #hdat = tfd.Get("h_zp_jigm").Clone()
    #hdat = newNameAndStructure(hdat,"data_obs",1,limrangelow,limrangehigh)
    
    #Old way of making data hist
    hdat = hdy.Clone()
    hdat.SetName("data_obs")
    hdat.Add(htt)
    hdat.Add(hvv)

    ####Do bkg stats unc explicitly
    #httstatsunc,httuncdict = doStatsUncertainty(htt)
    httstatsunc,httuncdict = gatherEmuStats(ttStats,desibins,limrangelow,limrangehigh)
    hvvstatsunc,hvvuncdict = doStatsUncertainty(hvv)

    ###Do alpha method unc
    dymethunchists,dymethshapedict,dymethratedict = gatherAlphaMethodUncs(dyEst,hdy,limrangelow,limrangehigh)

    ###Do emu tt extrp unc
    emuunchists,emushapedict,emuratedict = gatherEmuUncs("TT",ttEst,htt,limrangelow,limrangehigh)

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

        if "Zp5500ND1800NS200" not in signame:
            continue

        ####Make Files
        prepRootName = go.makeOutFile('Run2_'+yearstr+'_ZllHbbMET','extended_'+chan+'_'+signame,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        prepRootFile = ROOT.TFile(prepRootName,"recreate")

        ####Nominal Signal
        hsigori = signom[s]["tfile"].Get("h_zp_jigm")
        hsigori.Sumw2(ROOT.kTRUE)#Throws a warning that it is already created
        hsig = hsigori.Clone()
        hsig.Scale(signom[s]["scale"])
        hsig = applyStatsUncToSignal(hsig,signom[s]["errdf"]["h_zp_jigm"]*signom[s]["scale"],signom[s]["scale"])

        hsig = newNameAndStructure(hsig,signame,rebindiv,limrangelow,limrangehigh)
        prepRootFile.cd()


        #bkghists = [htt,hvv,hdy,hdat]
        bkghists = [htt,hvv,hdy]
        #print("!!!!!!!!artificially scaling the nominal backgrounds!!!!")
        #for hist in bkghists:
        #    hist.Scale(10)

        ###Write the Standard Stuff
        htt.Write()
        hvv.Write()
        hdy.Write()
        hdat.Write()
        hsig.Write()

        ###Do signal statistical uncs
        hsigstatsunc,hsiguncdict = doStatsUncertainty(hsig)

        ###Write bkg stats files
        print("about to do the stats uncs and this is the range  tt: ",len(httstatsunc))
        print("about to do the stats uncs and this is the range  vv: ",len(hvvstatsunc))
        print("about to do the stats uncs and this is the range sig: ",len(hsigstatsunc))
        
        for h in range(len(hvvstatsunc)):#need to use vv or sig, since tt has all bins
            #print("!!!!!!!!artificially scaling the nominal background stats unc!!!!")
            #httstatsunc[h].Scale(10)
            #hvvstatsunc[h].Scale(10)
            
            httstatsunc[h].Write()
            hvvstatsunc[h].Write()
            hsigstatsunc[h].Write()
            
        #####Gather Systematics
        systdictrate = {"lumi_13TeV":{"type":"lnN","unc":1.016,"proc":["1.016","1.016","1.016","1.016"]}
                    }
        systdictshape = {}

        ####Add the statistical uncertainites
        systdictshape.update(httuncdict)
        systdictshape.update(hvvuncdict)
        systdictshape.update(hsiguncdict)

        ####Add the alphamethod uncertainties
        systdictshape.update(dymethshapedict)
        systdictrate.update(dymethratedict)
        for hist in dymethunchists:
            #print("ARTIFICIALLY SCALING THE ALPHA METHOD HISTS")
            #hist.Scale(10)
            hist.Write()

        ###Add the emu uncertainties
        systdictshape.update(emushapedict)
        systdictrate.update(emuratedict)
        for hist in emuunchists:
            hist.Write()
        
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

            if len(systbkgsup.bkgs["WZTo2L2Q"][18]["sr"][0]) < 1:
                print("        There are no WZ entires for this systematic.")
                print("        moving on. This systematic will not included.")
                continue
            
            appcode = config.get(syst,'applist').split(',')
            ratenums = config.get(syst,'rate').split(',')
            ratelist = ratenums
            if (len(appcode) > 1):#If the shape systematic never gets applied
                applist = [float(x) for x in appcode]
                systdictshape[syst] = {"type":config.get(syst,'type'),"unc":1.0,"proc":applist}
                systdictrate[syst]  = {"type":config.get(syst,'typerate'),"proc":ratelist}
            else:
                systdictrate[syst]  = {"type":config.get(syst,'typerate'),"proc":ratelist}
            
            if rebindiv == 2:
                dyEstup     = ROOT.TFile(config.get(syst,'pathup')+'/Run2_161718_dy_extraploationalphat_'+config.get(syst,'strup')+'__Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
                dyEstdwn    = ROOT.TFile(config.get(syst,'pathdwn')+'/Run2_161718_dy_extraploationalphat_'+config.get(syst,'strdwn')+'__Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')

            ####Make it useful
            #sigup  = systsigup.getPreppedSig('sr',sigxs)#systsigup.prepsigsr
            #sigupdict = {}
            sigup = makeSignalInfoDict(systsigup,'sr',sigxs)
            sigdwn = makeSignalInfoDict(systsigdwn,'sr',sigxs)

            #Some sort of debug
            #keyup = sigup[s]["tfile"].GetListOfKeys()
            #keydwn = sigdwn[s]["tfile"].GetListOfKeys()
            #keyup = [k.GetName() for k in keyup]
            #keydwn = [k.GetName() for k in keydwn]

            #for k,key in enumerate(keyup):
            #    if "hnevents" in key:
            #        break
            #    print("Up key:   ",key)
            #    print("Down key: ",keydwn[k])
            #    print(" Up bins: ",sigup[s]["tfile"].Get(key).GetNbinsX())
            #    print("Dwn bins: ",sigdwn[s]["tfile"].Get(key).GetNbinsX())

            ####Check is the keys exist
            if (s in sigup) and (s in sigdwn): #Checks that the syst key exist for SIGNAl
                #Signal
                hsigupori = sigup[s]["tfile"].Get("h_zp_jigm")
                hsigup = hsigupori.Clone()
                hsigup.Scale(sigup[s]["scale"])
                hsigup = applyStatsUncToSignal(hsigup,sigup[s]["errdf"]["h_zp_jigm"]*sigup[s]["scale"],sigup[s]["scale"])
                hsigdwnori = sigdwn[s]["tfile"].Get("h_zp_jigm")
                hsigdwn = hsigdwnori.Clone()
                hsigdwn.Scale(sigdwn[s]["scale"])
                hsigdwn = applyStatsUncToSignal(hsigdwn,sigdwn[s]["errdf"]["h_zp_jigm"]*sigdwn[s]["scale"],sigdwn[s]["scale"])
                #Rename and Restructure
                hsigdwn = newNameAndStructure(hsigdwn,signame+"_"+syst+"Down",rebindiv,limrangelow,limrangehigh)
                hsigup = newNameAndStructure(hsigup,signame+"_"+syst+"Up",rebindiv,limrangelow,limrangehigh)
                #Get the proper signal rate systematics
                siguprate = getDeviatedOverNominal(hsigup,hsig)
                sigdwnrate = getDeviatedOverNominal(hsigdwn,hsig)
                sigratestr = str(round(sigdwnrate,4))+"/"+str(round(siguprate,4))
                if (round(sigdwnrate,4) == 1.0) or( round(siguprate,4) == 1.0):
                    sigratestr = "-"
                ratelist[0] = sigratestr

                #Write the Signal Hists
                hsigup.Write()
                hsigdwn.Write()
            else:
                print("****************{0} Does not have {1} systematic for signal".format(s,syst))

                
            ####Prepping holders####
            empty4 = empty.Clone()
            empty5 = empty.Clone()
            empty6 = empty.Clone()
            empty7 = empty.Clone()
            empty8 = empty.Clone()
            empty9 = empty.Clone()

            ####Gathering the Systematic
            #Background
            #httup  = systbkgsup.getAddedHist(empty4,"TT","sr","h_zp_jigm")
            #httdwn = systbkgsdwn.getAddedHist(empty5,"TT","sr","h_zp_jigm")
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

            
            ####Rename and Restructure
            #httup = newNameAndStructure(httup,"TT_"+syst+"Up",rebindiv,limrangelow,limrangehigh)
            #httdwn = newNameAndStructure(httdwn,"TT_"+syst+"Down",rebindiv,limrangelow,limrangehigh)
            hvvup = newNameAndStructure(hvvup,"VV_"+syst+"Up",rebindiv,limrangelow,limrangehigh) 
            hvvdwn = newNameAndStructure(hvvdwn,"VV_"+syst+"Down",rebindiv,limrangelow,limrangehigh)
            hdyup  = newNameAndStructure(hdyup,"DY_"+syst+"Up",1,limrangelow,limrangehigh)
            hdydwn = newNameAndStructure(hdydwn,"DY_"+syst+"Down",1,limrangelow,limrangehigh)
            
            #allsysthists = [httup,httdwn,hvvup,hvvdwn,hdyup,hdydwn]
            allsysthists = [hvvup,hvvdwn,hdyup,hdydwn]
            allthesame = checkBinNumbers(allsysthists)
            if not allthesame:
                print("              NOT ALL OF THE BINS NUMBERS ARE THE SAME")
                print("              NOT WRITING THE {0} UP/DWN HISTOGRAMS".format(syst))
                continue

                
            #Write the histograms
            prepRootFile.cd()
            print("        Writing the up/dwn histograms")
                
            #bkghists = [httup,hvvup,hdyup,httdwn,hvvdwn,hdydwn]
            bkghists = [hvvup,hdyup,hvvdwn,hdydwn]
            print("!!!!!!!!!!!!put back in ttbar up/downs!!!!!!!!!!!!!!!!!!!!!!!")
            #print("!!!!!!!!artificially scaling the up/dwn backgrounds!!!!")
            #for hist in bkghists:
            #    hist.Scale(10)
            #httup.Write()
            #httdwn.Write()
            hvvup.Write()
            hvvdwn.Write()
            hdyup.Write()
            hdydwn.Write()

            #plotSignalSystematics(hsigup,hsigdwn,hsig,signame+"_"+syst)

        #For writing the datacard
        procdict = {"processnames":[signame,"DY","TT","VV"],
                    "hists":[hsig,hdy,htt,hvv],
                    "method":["mc","alpha","mc","mc"],
                    "systshapes":systdictshape,
                    "systrates":systdictrate,
        }
        print("        Defined the Datacard Dict")
        print("        Writing the Datacard")

        #non-extended ignoring
        #ignorerates = ["extrap_subdatasb_vvfit_par1","extrap_subdatasb_vvfit_par0","extrap_subdatasb_ttfit_par0"]
        #ignoreshapes = ["extrap_subdatasb_vvfit_par1","extrap_subdatasb_vvfit_par0","extrap_subdatasb_ttfit_par1","extrap_subdatasb_ttfit_par0","extrap_alpha_DY_sr_par0","extrap_alpha_DY_sb_par1","extrap_alpha_DY_sb_par0","muonid","TT_scale"]
        ignorerates = []
        ignoreshapes = ["muonid","TT_scale"]
        writeDataCard(procdict,prepRootName,chan,yearstr,hdat,ignorerates,ignoreshapes)
        print("        About to close the root file")
        prepRootFile.Close()

