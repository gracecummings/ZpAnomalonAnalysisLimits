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
    #alpha = 1.- 0.682689492
    for ibin in range(hist.GetNbinsX()+1):
        if ibin == 0:
            continue
        #elif hist.GetBinContent(ibin) > 0:
        binerr = errseries[ibin-1]
        hist.SetBinError(ibin,binerr)
        #else:
        #    binerrbase = ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin))+1,1)-hist.GetBinContent(ibin)
        #    binerr = binerrbase*scale
        #    hist.SetBinError(ibin,binerr)

            
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

def gatherEmuStats(stattf,desiredbins,limrangelow,limrangehigh,newedges):
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
    #print(len(newedges))
    for key in keys:
        oribin = key.split('StatsUncBin')[-1]
        if int(oribin) in desiredbins:#only taking hists for bins we care about
            newbinnum = int(oribin)-(int(min(desibins)-1))#move lowest bin # to 1
            if newbinnum > (len(newedges)-1):#One less bin than the bin edges
                continue
            h = stattf.Get(key)#get histogram of interest
            #print("Number of bins in the ttstats hist: ",h.GetNbinsX())
            direc = key.split('_TT_TT')[0]
            newname = "TT_TT_StatsUncBin"+str(newbinnum)+direc#new hist name
            trimmedh = newNameAndStructure(h,newname,1,limrangelow,limrangehigh)#only take bins we care about in hist
            trimmedh = trimmedh.Rebin(len(newbinedges)-1,newname,newbinedges)
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
                hdwn.SetBinContent(b,0)
        else:
            hup.SetBinContent(b,bincont+binunc)#zero in bins, garwood interval scaled should be err

        hup.SetName(hname+"Up")
        hdwn.SetName(hname+"Down")
        histlist.append(hup)
        histlist.append(hdwn)
        uncdict  = {name+"_StatsUncBin"+str(b):{"type":"shape","unc":1.0,"proc":applist}}
        statsdict.update(uncdict)
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
    prepCardName = go.makeOutFile('Run2_'+yearstr+'_ZllHbbMET','datacard_interpolation_'+chan+'_'+signame,'.txt',str(zptcut),str(hptcut),str(metcut),str(btagwp))
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


def makeBinLowEdges(hist,lastnormbin):
    #this takes the histogram ou want to rebin, and you just give it the last normal bin
    nbins = hist.GetNbinsX()
    binw  = hist.GetBinWidth(1)
    iniedges = [hist.GetBinLowEdge(i) for i in range(nbins+2)]#+2 for the actual highedge
    if lastnormbin not in iniedges:
        print("desired bin edge not possible, try again")
    newedges = [edg for edg in iniedges if (edg <= lastnormbin)]
    newedges.append(iniedges[-1])
    newedges = newedges[1:]
    return np.array(newedges)


def getBinNumbersOfInterest(hist,rebindiv,limrangelow,limrangehigh):
    #for the ttbar stats gathering
    hist.Rebin(rebindiv)
    nbins = hist.GetNbinsX()
    binw  = hist.GetBinWidth(1)
    desiedges = [limrangelow+x*binw for x in range(int((limrangehigh-limrangelow)/binw))]
    desibins  = [edge/binw+1 for edge in desiedges]
    setdesi = set(desibins)
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
        up = up.Rebin(len(newbinedges)-1,"DY_"+upkey,newbinedges)
        dn = dn.Rebin(len(newbinedges)-1,"DY_"+dnkey,newbinedges)
        
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

def gatherEmuUncs(procname,tf,hnom,limrangelow,limrangehigh,newbinedges):
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
        up = up.Rebin(len(newbinedges)-1,procname+"_"+upkey,newbinedges)
        dn = dn.Rebin(len(newbinedges)-1,procname+"_"+dnkey,newbinedges)

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
    limrangehigh = 5000

    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()

    ####load in the files with the nominal distributions
    bkgs = go.backgrounds(config.get('nominal','pathnom'),zptcut,hptcut,metcut,btagwp,config.get('nominal','strnom'))
    #sig  = go.signal(config.get('nominal','pathsignom'),zptcut,hptcut,metcut,btagwp,sigxs,[16,17,18],config.get('nominal','strnom'))
    sig   = ROOT.TFile('analysis_output_ZpAnomalon/2022-10-05/interpolation_'+config.get('nominal','strnom')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
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

    ####Rebin with the new edges
    newbinedges = makeBinLowEdges(hvv,2800)#3000 is last bin we want to be normal
    hvv = hvv.Rebin(len(newbinedges)-1,"VV",newbinedges)
    htt = htt.Rebin(len(newbinedges)-1,"TT",newbinedges)
    hdy = hdy.Rebin(len(newbinedges)-1,"DY",newbinedges)

    ####Make data obs hist - GoF test
    #tfd = ROOT.TFile(config.get('nominal','pathnom')+"/Run2_161718_gofplots_"+config.get('nominal','strnom')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
    #hdat = tfd.Get("h_zp_jigm").Clone()
    #hdat = newNameAndStructure(hdat,"data_obs",1,limrangelow,limrangehigh)
    #hdat = hdat.Rebin(len(newbinedges)-1,"data_obs",newbinedges)
    
    #Old way of making data hist
    hdat = hdy.Clone()
    hdat.SetName("data_obs")
    hdat.Add(htt)
    hdat.Add(hvv)

    ####Do bkg stats unc explicitly
    #httstatsunc,httuncdict = doStatsUncertainty(htt)
    httstatsunc,httuncdict = gatherEmuStats(ttStats,desibins,limrangelow,limrangehigh,newbinedges)
    hvvstatsunc,hvvuncdict = doStatsUncertainty(hvv)#appropriately handling VV MS errs
    
    ###Do alpha method unc
    dymethunchists,dymethshapedict,dymethratedict = gatherAlphaMethodUncs(dyEst,hdy,limrangelow,limrangehigh)

    ###Do emu tt extrp unc
    emuunchists,emushapedict,emuratedict = gatherEmuUncs("TT",ttEst,htt,limrangelow,limrangehigh,newbinedges)

    ####For Each signal, make a datacard, and a root file with all systematics
    #siginfo = sig.getPreppedSig('sr',sigxs)
    #signom = makeSignalInfoDict(sig,'sr',sigxs)
    #sigcolors = go.colsFromPalette(siginfo,ROOT.kCMYK)
    #nomsigs = [s["name"] for s in siginfo]

    sigkeys = sig.GetListOfKeys()
    nomsigs = [key.GetName() for key in sigkeys]

    for name in nomsigs:
        #name = signom[s]["name"]
        #signame = "holder"
        #if "Tune" in name:
            #strippedname = name.split("_Tune")[0]
            #signame = strippedname.replace("-","")
       # else:
            #signame = name.replace("-","")
        print("------- Looking at signal sample ",name)

        #if "Zp4000ND600NS200" not in name:
        #    continue
        if "NS1" in name:
            continue

        ####Make Files
        prepRootName = go.makeOutFile('Run2_'+yearstr+'_ZllHbbMET','interpolation_'+chan+'_'+name,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        prepRootFile = ROOT.TFile(prepRootName,"recreate")

        ####Nominal Signal
        hsigori = sig.Get(name)#signom[s]["tfile"].Get("h_zp_jigm")
        #hsigori.Sumw2(ROOT.kTRUE)#Throws a warning that it is already created
        hsig = hsigori.Clone()
        scale = go.findScale(35000,137.6,sigxs)#interpolation based on 35000...
        hsig.Scale(scale)
        #hsig = applyStatsUncToSignal(hsig,signom[s]["errdf"]["h_zp_jigm"]*signom[s]["scale"],signom[s]["scale"])#no stats unc on interpolation

        hsig = newNameAndStructure(hsig,name,1,limrangelow,limrangehigh)
        hsig = hsig.Rebin(len(newbinedges)-1,name,newbinedges)
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

        ###Write bkg stats files
        print("about to do the stats uncs and this is the range  tt: ",len(httstatsunc))
        print("about to do the stats uncs and this is the range  vv: ",len(hvvstatsunc))
        
        for h in range(len(hvvstatsunc)):#need to use vv or sig, since tt has all bins
            httstatsunc[h].Write()
            hvvstatsunc[h].Write()
             
        #####Gather Systematics
        systdictrate = {"lumi_13TeV":{"type":"lnN","unc":1.016,"proc":["1.016","1.016","1.016","1.016"]},
                        "prefire":{"type":"lnN","unc":1.05,"proc":["1.05","-","-","-"]}
                    }
        systdictshape = {}

        ####Add the statistical uncertainites
        systdictshape.update(httuncdict)
        systdictshape.update(hvvuncdict)
 
        ####Add the alphamethod uncertainties
        systdictshape.update(dymethshapedict)
        systdictrate.update(dymethratedict)
        for hist in dymethunchists:
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
            systsigup   = ROOT.TFile('analysis_output_ZpAnomalon/2022-10-05/interpolation_'+config.get(syst,'strup')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
            systsigdwn  = ROOT.TFile('analysis_output_ZpAnomalon/2022-10-05/interpolation_'+config.get(syst,'strdwn')+'_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')

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
            #sigup = systsigup.Get(name)
           # sigdwn = systsigdwn.Get(name)

            sigsystshists = []

            #Signal
            hsigupori = systsigup.Get(name)
            hsigup = hsigupori.Clone()
            hsigup.SetDirectory(0)
            systsigup.Close()
            hsigup.Scale(scale)
            hsigdwnori = systsigdwn.Get(name)
            hsigdwn = hsigdwnori.Clone()
            hsigdwn.SetDirectory(0)
            systsigdwn.Close()
            hsigdwn.Scale(scale)

            #Rename and Restructure
            hsigdwn = newNameAndStructure(hsigdwn,name+"_"+syst+"Down",1,limrangelow,limrangehigh)
            #print('Down signal integral after rebin1 ',hsigdwn.Integral())
            hsigup = newNameAndStructure(hsigup,name+"_"+syst+"Up",1,limrangelow,limrangehigh)
            hsigup  = hsigup.Rebin(len(newbinedges)-1,name+"_"+syst+"Up",newbinedges)
            hsigdwn = hsigdwn.Rebin(len(newbinedges)-1,name+"_"+syst+"Down",newbinedges)

            #Get the proper signal rate systematics
            siguprate = getDeviatedOverNominal(hsigup,hsig)
            sigdwnrate = getDeviatedOverNominal(hsigdwn,hsig)
            sigratestr = str(round(sigdwnrate,4))+"/"+str(round(siguprate,4))
            if (round(sigdwnrate,4) == 1.0) and (round(siguprate,4) == 1.0):
                sigratestr = "-"
            if (round(sigdwnrate,4) == 1.0):
                sigratestr = str(round(siguprate,4))
            if (round(siguprate,4) == 1.0):
                sigratestr = str(round(sigdwnrate,4))
            ratelist[0] = sigratestr

            #CollectSignal Hists
            sigsystshists.append(hsigup)
            sigsystshists.append(hsigdwn)
                
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
            hvvup = hvvup.Rebin(len(newbinedges)-1,"VV_"+syst+"Up",newbinedges)
            hvvdwn = hvvdwn.Rebin(len(newbinedges)-1,"VV_"+syst+"Down",newbinedges)
            hdyup = hdyup.Rebin(len(newbinedges)-1,"DY_"+syst+"Up",newbinedges)
            hdydwn = hdydwn.Rebin(len(newbinedges)-1,"DY_"+syst+"Down",newbinedges)
            
            #allsysthists = [httup,httdwn,hvvup,hvvdwn,hdyup,hdydwn]
            allsysthists = [hdy,hvv,hvvup,hvvdwn,hdyup,hdydwn]
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
            histstosave = sigsystshists+bkghists
            print("!!!!!!!!!!!!put back in ttbar up/downs!!!!!!!!!!!!!!!!!!!!!!!")
            #print("!!!!!!!!artificially scaling the up/dwn backgrounds!!!!")
            #for hist in bkghists:
            #    hist.Scale(10)

            for hist in histstosave:
                hist.Write()
            #httup.Write()
            #httdwn.Write()
            #hvvup.Write()
            #hvvdwn.Write()
            #hdyup.Write()
            #hdydwn.Write()

            #plotSignalSystematics(hsigup,hsigdwn,hsig,name+"_"+syst)

        #For writing the datacard
        procdict = {"processnames":[name,"DY","TT","VV"],
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

