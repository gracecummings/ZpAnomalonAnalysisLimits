import os
import sys
import glob
import ROOT
import configparser
import pandas as pd
import numpy as np
from datetime import date

def getDeviatedOverNominal(hist,histnom):
    interest = hist.Integral()
    intnom = histnom.Integral()
    devonom = interest/intnom
    return devonom

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
        nh.SetBinError(b+1,binerror)#

    return nh

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

def makeSignalErrorDictionary(flist):
    sigsrerrsdict = {}
    for sig in flist:
        name1 = nameSignal(sig)
        name = name1.split("_selected")[0]#This is the name to match the errs to files
        sigsrerrsdict[name] = sig
    return sigsrerrsdict
    

def lumiFormatter(yearlist):
    lumidict = {16:36.31,17:41.53,18:59.74}
    lumi = 0
    for year in yearlist:
        lumi += lumidict[year]

    lumi = round(lumi,2)
    lumistr = str(lumi)+" fb^{-1}"
    return lumistr

def yearFormatter(yearlist):
    yearstr =''
    for year in yearlist:
        yearstr = yearstr+str(year)
    return yearstr


def sampleType(sampstring,givejecs=False):
    #Make numerical code for type of sample
    if "Run" in sampstring:
        if "SingleMuon" in sampstring:
            samptype = -1
        elif ("SingleElectron" in sampstring) or ("EGamma" in sampstring):
            samptype = -2
        else:
            samptype = 0
    elif "ZpAnomalon" in sampstring:
        samptype = 1
    elif "DYJetsToLL" in sampstring:
         samptype = 2
    elif "TTTo" in sampstring:
        samptype = 3
    elif "WZTo" in sampstring:
        samptype = 4
    elif "ZZTo" in sampstring:
        samptype = 5
    elif "WWTo" in sampstring:
        samptype = 6

    else:
        samptype = -1000

    if "2018" in sampstring:
        year = 18
        if givejecs:
            if "Run2018A" in sampstring:
                year = 180
            if "Run2018B" in sampstring:
                year = 181
            if "Run2018C" in sampstring:
                year = 182
            if "Run2018D" in sampstring:
                year = 183
    if "Autumn18" in sampstring:
        year = 18
    if "2017" in sampstring:
        year = 17
        if givejecs:
            if "Run2017B" in sampstring:
                year = 170
            if "Run2017C" in sampstring:
                year = 171
            if "Run2017D" in sampstring:
                year = 172
            if "Run2017E" in sampstring:
                year = 172
            if "Run2017F" in sampstring:
                year = 173
    if "Fall17" in sampstring:
        year = 17
    if "2016" in sampstring:
        year = 16
        if givejecs:
            if "Run2016B" in sampstring:
                year = 160
            if "Run2016C" in sampstring:
                year = 160
            if "Run2016D" in sampstring:
                year = 160
            if "Run2016E" in sampstring: 
                year = 161
            if "Run2016F" in sampstring:
                year = 161
            if "Run2016G" in sampstring:
                year = 162
            if "Run2016H" in sampstring:
                year = 162
    if "Summer16" in sampstring:
        year = 16
    return samptype,year

def makeOutFile(sampstring,descrip,ftype,zptcut,hptcut,metcut,btagwp):
    if not os.path.exists("analysis_output_ZpAnomalon/"+str(date.today())+"/"):
        os.makedirs("analysis_output_ZpAnomalon/"+str(date.today())+"/")
    outFile = "analysis_output_ZpAnomalon/"+str(date.today())+"/"+sampstring+"_"+descrip+"_Zptcut"+zptcut+"_Hptcut"+hptcut+"_metcut"+metcut+"_btagwp"+btagwp+ftype
    return outFile

def orderDY(histFile):
    s0 = histFile.split("/")[-1]
    s1 = s0.split("to")[0]
    s2 = s1.split("HT-")[1]
    return int(s2)       

def orderFall17TT(histFile):#NROKEN AT THE MOMENT
    s1 = histFile.split("Events")[0]
    s2 = s1.split("_")[-1]
    return int(s2)

def massPoints(nameSig):#do this for full name
    s1  = nameSig.split("-")
    mzp = int(s1[0].split("Zp")[1])
    mnd = int(s1[1].split("ND")[1])
    if "Tune" in nameSig:
        ms1 = s1[2].split("_Tune")[0]
        mns = int(ms1.split("NS")[1])
    else:
        ms1 = s1[2].split("_upout")[0]
        mns = int(ms1.split("NS")[1])
    return mzp,mnd,mns

def nameSignal(histFile):
    s1 = histFile.split("ZpAnomalonHZ_UFO-")[1]
    s2 = s1.split("_upout")[0]
    return s2

def findScale(prodnum,lumi,xsec):
    expecnum = xsec*lumi
    scalefac = expecnum/prodnum
    return  scalefac

def colsFromPalette(samplist,palname):
    collist = []
    ROOT.gStyle.SetPalette(palname)
    cols = ROOT.TColor.GetPalette()
    colsnum = cols.GetSize()
    for i in range(len(samplist)):
        collist.append(cols.At(0+i*int(colsnum/len(samplist))))
    collist.reverse()
    return collist

def gatherBkg(bkg_dir,descrip,zptcut,hptcut,metcut,btagwp,year):
    if year == 18:
        mcprefix = 'Autumn18'
    if year == 17:
        mcprefix = 'Fall17'
        
    DYJetsToLL = glob.glob(str(bkg_dir)+'/'+mcprefix+'.DYJetsToLL_M-50_HT*'+descrip+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
    TT         = glob.glob(str(bkg_dir)+'/'+mcprefix+'.TTT*_'+descrip+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')                                         
    WZTo2L2Q   = glob.glob(str(bkg_dir)+'/'+mcprefix+'.WZTo2L2Q*'+descrip+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')                                    
    ZZTo2L2Q   = glob.glob(str(bkg_dir)+'/'+mcprefix+'.ZZTo2L2Q*'+descrip+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')                                    
    bkgfiles   = [DYJetsToLL,TT,WZTo2L2Q,ZZTo2L2Q]
    return bkgfiles

def prepSig(sigfiles,sig_colors,sig_xsec,lumi,sigerrs):
    sig_info = []
    print("          Fix your hack for ntuple signal counting, setting sample size to 35000")
    for s,sig in enumerate(sigfiles):
        sig_dict = {}
        sig_dict["tfile"] = ROOT.TFile.Open(sig,"read")
        sig_samplesize    = 0.0#str(sig_dict["tfile"].Get('hnevents').GetString())
        if float(sig_samplesize) == 0.0:
            sig_samplesize = "35000.0"###HACK
        sig_dict["scale"] = findScale(float(sig_samplesize),sig_xsec,lumi)
        sig_dict["name"]  = nameSignal(sig)
        mzp,mnd,mns       = massPoints(sig_dict["name"])
        sig_dict["mzp"]   = mzp#
        sig_dict["mnd"]   = mnd
        sig_dict["mns"]   = mns

        #Do the errors
        sigerrfilename = sigerrs[sig_dict["name"]]
        sigerrdf = pd.read_pickle(sigerrfilename)
        ssigerrdf = sigerrdf#*sig_dict["scale"]
        sig_dict["errdf"] = ssigerrdf
        sig_info.append(sig_dict)                                                                    
    #Sort Signals by ND mass, then by Zp mass                                                           
    sig_info = sorted(sig_info,key = lambda sig: (sig["mnd"],sig["mzp"],sig["mns"]))                    
    for s,sig in enumerate(sig_info):                                                                   
        sig["color"] = sig_colors[s]

    return sig_info

def prepBkg(bkgfiles,bkgnames,bkg_colors,ini_file,lumi,flag="yes"):
    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open(ini_file)
    config.read_file(fp)
    bkg_info = []
    for b,bkg in enumerate(bkgfiles):
        bkg_binsum   = {}
        bkg_binlist  = []
        bkg_channel  = bkgnames[b]
        bkg_expyield = 0
        # bkg xs from .ini file
        bkgbin_xs_pairs = config.items(bkg_channel)
        normscale = 1.
        if bkg_channel == "DYJetsToLL":
            #orders smallest HT to largest
            bkg.sort(key = orderDY)
            #normscale = np.load('BkgInputsNominalJECBtag/Run2_2017_2018_dynormalization_systnominal_btagnom_signalblind_Zptcut150.0_Hptcut300.0_metcut200.0_btagwp0.8.npy')[0]
            #normscale = np.load('pfMETNominal_CorrCounting_Optimization/Run2_2017_2018_dynormalization_systnominal_btagnom_signalblind_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.npy')[0]
            #print("applying a normalization,watch out")
            normscale = 1
            print("Put back in derived norm")
        elif bkg_channel == "TT":
            #sorts in alphabetical order 
            bkg.sort()                                                     
        else:
            if flag == "no":
                break
        #loop through each process bin or categrory
        for s,bkgbin in enumerate(bkg):
            bkgbin_dict = {}
            #print(bkgbin)
            bkgbin_dict["binname"] = bkgbin_xs_pairs[s][0]
            #print(bkgbin_dict["binname"])
            bkgbin_dict["tfile"]   = ROOT.TFile.Open(bkgbin)
            bkgbin_sampsize        = str(bkgbin_dict["tfile"].Get('hnevents').GetString())
            bkgbin_xs              = float(bkgbin_xs_pairs[s][1].split()[0])*1000#Into Femtobarn
            #print(bkgbin_sampsize)
            #print(bkgbin_xs)
            bkgbin_dict["scale"]   = findScale(float(bkgbin_sampsize),bkgbin_xs,lumi)*normscale
            bkgbin_dict["color"]   = bkg_colors[b]
            #get the number of passing events
            bkgbin_yield           = float(str(bkgbin_dict["tfile"].Get('hnevents_pZ').GetString()))#last cut hist ##goes back to btag
            bkg_expyield          += bkgbin_yield*bkgbin_dict["scale"]
            bkg_binlist.append(bkgbin_dict)
            bkg_binsum["expyield"] = bkg_expyield

        bkg_binsum["binlist"] = bkg_binlist
        bkg_binsum["name"]    = bkg_channel
        bkg_info.append(bkg_binsum)

    #Sort the backgrounds from the smallest yields to largest        
    bkg_info = sorted(bkg_info, key = lambda bkg:bkg["expyield"])
    return bkg_info
#
def stackBkg(bkg_info,hist_to_stack,hsbkg,legend,stack_max,stack_min):
    for bkg in bkg_info:
        for b,bkgbin in enumerate(bkg["binlist"]):
            hbkg = bkgbin["tfile"].Get(hist_to_stack)
            hbkg.SetStats(0)
            hbkg.Scale(bkgbin["scale"])
            hbkg.SetFillColor(bkgbin["color"])
            hbkg.SetLineColor(bkgbin["color"])
            hbkg.SetMaximum(stack_max)
            hbkg.SetMinimum(stack_min)
            hsbkg.Add(hbkg)
            hsbkg.Draw("HIST")
            hsbkg.SetMaximum(stack_max)
            hsbkg.SetMinimum(stack_min)
            if b == len(bkg["binlist"])-1:
                legend.AddEntry(hbkg,bkg["name"],"f")
                
def stackBkgMultiYear(bkginfo0,bkginfo1,hist_to_stack,hsbkg,legend,stack_max,stack_min):
    for i,bkg in enumerate(bkginfo0):
        bkgdictmatch = bkginfo1[i]
        if bkg["name"] != bkginfo1[i]["name"]:
            #print("background orders do not match!!!!!")
            bkgdictmatch = [bkgdict for bkgdict in bkginfo1 if bkgdict["name"] == bkg["name"]][0]
        for b,bkgbin in enumerate(bkg["binlist"]):
            #print("new bin!")
            #print(bkgbin["binname"])
            #print(bkginfo1[i]["binlist"][b]["binname"])
            #print(bkgdictmatch["binlist"][b]["binname"])
            hbkg0 = bkgbin["tfile"].Get(hist_to_stack)
            hbkg0.SetStats(0)
            hbkg0.Scale(bkgbin["scale"])
            hbkg0.SetFillColor(bkgbin["color"])
            hbkg0.SetLineColor(bkgbin["color"])
            hbkg0.SetMaximum(stack_max)
            hbkg0.SetMinimum(stack_min)
            hbkg1 = bkgdictmatch["binlist"][b]["tfile"].Get(hist_to_stack)
            hbkg1.SetStats(0)
            hbkg1.Scale(bkgdictmatch["binlist"][b]["scale"])
            hbkg1.SetFillColor(bkgdictmatch["binlist"][b]["color"])
            hbkg1.SetLineColor(bkgdictmatch["binlist"][b]["color"])
            hbkg1.SetMaximum(stack_max)
            hbkg1.SetMinimum(stack_min)

            hsbkg.Add(hbkg0)
            hsbkg.Add(hbkg1)
            hsbkg.Draw("HIST")
            hsbkg.SetMaximum(stack_max)
            hsbkg.SetMinimum(stack_min)
            
            if b == len(bkg["binlist"])-1:
                legend.AddEntry(hbkg1,bkg["name"],"f")

def saveNpUncertainties(uncdf,filename):
    npF = open(filename,'wb')
    np.savez(npF,
             h_z_pt  = uncdf['h_z_pt'].values,
             h_z_eta = uncdf['h_z_eta'].values,
             h_z_m   = uncdf['h_z_m'].values,
             h_h_pt  = uncdf['h_h_pt'].values,
             h_h_eta = uncdf['h_h_eta'].values,
             h_h_m   = uncdf['h_h_m'].values,
             h_h_sd  = uncdf['h_h_sd'].values,
             h_met   = uncdf['h_met'].values,
             h_zp_jigm = uncdf['h_zp_jigm'].values,
             h_nd_jigm = uncdf['h_nd_jigm'].values,
             h_ns_jigm = uncdf['h_ns_jigm'].values,
             h_btag    = uncdf['h_btag'].values
             )

class backgrounds_local:
    def __init__(self,path,zptcut,hptcut,metcut,btagwp,systr=""):
        self.path = path
        self.zptcut = zptcut
        self.hptcut = hptcut
        self.metcut = metcut
        self.btagwp = btagwp
        self.systr = systr

        if systr != "":
            systr = "_"+systr

        #gather background MC files
        #2017
        self.f17dyjetsb = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17dyjetsr = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17dyjettr = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsb = glob.glob(str(path)+'/Fall17.TTT*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsr = glob.glob(str(path)+'/Fall17.TTT*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17tttr = glob.glob(str(path)+'/Fall17.TTT*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsb = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsr = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wztr = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsb = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsr = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zztr = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')

        #2018
        self.a18dyjetsb = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18dyjetsr = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18dyjettr = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsb = glob.glob(str(path)+'/Autumn18.TTT*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsr = glob.glob(str(path)+'/Autumn18.TTT*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18tttr = glob.glob(str(path)+'/Autumn18.TTT*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsb = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsr = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wztr = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsb = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsr = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zztr = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')

        #2016
        self.s16dyjetsb = glob.glob(str(path)+'/Summer16v3.DYJetsToLL_M-50_HT*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16dyjetsr = glob.glob(str(path)+'/Summer16v3.DYJetsToLL_M-50_HT*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16dyjettr = glob.glob(str(path)+'/Summer16v3.DYJetsToLL_M-50_HT*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16ttsb = glob.glob(str(path)+'/Summer16v3.TTT*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16ttsr = glob.glob(str(path)+'/Summer16v3.TTT*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16tttr = glob.glob(str(path)+'/Summer16v3.TTT*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16wzsb = glob.glob(str(path)+'/Summer16v3.WZTo2L2Q*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16wzsr = glob.glob(str(path)+'/Summer16v3.WZTo2L2Q*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16wztr = glob.glob(str(path)+'/Summer16v3.WZTo2L2Q*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16zzsb = glob.glob(str(path)+'/Summer16v3.ZZTo2L2Q*_upout_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16zzsr = glob.glob(str(path)+'/Summer16v3.ZZTo2L2Q*_upout_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16zztr = glob.glob(str(path)+'/Summer16v3.ZZTo2L2Q*_upout_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')



        #gather errors
        self.f17dyjetsberrs = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17dyjetsrerrs = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17dyjettrerrs = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsberrs = glob.glob(str(path)+'/Fall17.TTT*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsrerrs = glob.glob(str(path)+'/Fall17.TTT*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17tttrerrs = glob.glob(str(path)+'/Fall17.TTT*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsberrs = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsrerrs = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wztrerrs = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsberrs = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsrerrs = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zztrerrs = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')

        self.a18dyjetsberrs = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18dyjetsrerrs = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18dyjettrerrs = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsberrs = glob.glob(str(path)+'/Autumn18.TTT*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsrerrs = glob.glob(str(path)+'/Autumn18.TTT*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18tttrerrs = glob.glob(str(path)+'/Autumn18.TTT*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsberrs = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsrerrs = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wztrerrs = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsberrs = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsrerrs = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zztrerrs = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        
        self.s16dyjetsberrs = glob.glob(str(path)+'/Summer16v3.DYJetsToLL_M-50_HT*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16dyjetsrerrs = glob.glob(str(path)+'/Summer16v3.DYJetsToLL_M-50_HT*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16dyjettrerrs = glob.glob(str(path)+'/Summer16v3.DYJetsToLL_M-50_HT*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16ttsberrs = glob.glob(str(path)+'/Summer16v3.TTT*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16ttsrerrs = glob.glob(str(path)+'/Summer16v3.TTT*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16tttrerrs = glob.glob(str(path)+'/Summer16v3.TTT*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16wzsberrs = glob.glob(str(path)+'/Summer16v3.WZTo2L2Q*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16wzsrerrs = glob.glob(str(path)+'/Summer16v3.WZTo2L2Q*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16wztrerrs = glob.glob(str(path)+'/Summer16v3.WZTo2L2Q*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16zzsberrs = glob.glob(str(path)+'/Summer16v3.ZZTo2L2Q*_selected_errors_sideband'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16zzsrerrs = glob.glob(str(path)+'/Summer16v3.ZZTo2L2Q*_selected_errors_signalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16zztrerrs = glob.glob(str(path)+'/Summer16v3.ZZTo2L2Q*_selected_errors_totalr'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')


        self.bkgs = {"DYJetsToLL":
                {18:
                 {"sb":[self.a18dyjetsb,self.a18dyjetsberrs],
                  "sr":[self.a18dyjetsr,self.a18dyjetsrerrs],
                  "tr":[self.a18dyjettr,self.a18dyjettrerrs]},
                 17:
                 {"sb":[self.f17dyjetsb,self.f17dyjetsberrs],
                  "sr":[self.f17dyjetsr,self.f17dyjetsrerrs],
                  "tr":[self.f17dyjettr,self.f17dyjettrerrs]},
                 16:
                 {"sb":[self.s16dyjetsb,self.s16dyjetsberrs],
                  "sr":[self.s16dyjetsr,self.s16dyjetsrerrs],
                  "tr":[self.s16dyjettr,self.s16dyjettrerrs]},
                },
                "TT":
                {18:
                 {"sb":[self.a18ttsb,self.a18ttsberrs],
                  "sr":[self.a18ttsr,self.a18ttsrerrs],
                  "tr":[self.a18tttr,self.a18tttrerrs]},
                 17:
                 {"sb":[self.f17ttsb,self.f17ttsberrs],
                  "sr":[self.f17ttsr,self.f17ttsrerrs],
                  "tr":[self.f17tttr,self.f17tttrerrs]},
                 16:
                 {"sb":[self.s16ttsb,self.s16ttsberrs],
                  "sr":[self.s16ttsr,self.s16ttsrerrs],
                  "tr":[self.s16tttr,self.s16tttrerrs]},
                },
                "WZTo2L2Q":
                {18:
                 {"sb":[self.a18wzsb,self.a18wzsberrs],
                  "sr":[self.a18wzsr,self.a18wzsrerrs],
                  "tr":[self.a18wztr,self.a18wztrerrs]},
                 17:
                 {"sb":[self.f17wzsb,self.f17wzsberrs],
                  "sr":[self.f17wzsr,self.f17wzsrerrs],
                  "tr":[self.f17wztr,self.f17wztrerrs]},
                 16:
                 {"sb":[self.s16wzsb,self.s16wzsberrs],
                  "sr":[self.s16wzsr,self.s16wzsrerrs],
                  "tr":[self.s16wztr,self.s16wztrerrs]},
                },
                "ZZTo2L2Q":
                {18:
                 {"sb":[self.a18zzsb,self.a18zzsberrs],
                  "sr":[self.a18zzsr,self.a18zzsrerrs],
                  "tr":[self.a18zztr,self.a18zztrerrs]},
                 17:
                 {"sb":[self.f17zzsb,self.f17zzsberrs],
                  "sr":[self.f17zzsr,self.f17zzsrerrs],
                  "tr":[self.f17zztr,self.f17zztrerrs]},
                 16:
                 {"sb":[self.s16zzsb,self.s16zzsberrs],
                  "sr":[self.s16zzsr,self.s16zzsrerrs],
                  "tr":[self.s16zztr,self.s16zztrerrs]},
                }
        }

        self.config = configparser.RawConfigParser()
        self.config.optionxform = str
        fp = open('xsects_2017.ini')#2017 and 2018dy+jets xs same
        self.config.read_file(fp)

        
    def getAddedHist(self,hist,samp,region,hname,years = [16,17,18]):
        bkg = self.bkgs[samp]
        xspairs = self.config.items(samp)
        #print(xspairs)
        bkgdfs  = []
        scales = []
        for year in years:
            if year == 16:
                lumi = 36.31
            if year == 17:
                lumi = 41.53
            if year == 18:
                lumi = 59.74
                
            files = bkg[year][region][0]
            errs  = bkg[year][region][1]
            if "DYJetsToLL" in samp:
                files.sort(key = orderDY)
                errs.sort(key = orderDY)
            if "TT" in samp:
                files.sort()
                errs.sort()

            #print(files)
            for i,f in enumerate(files):
                fparts = f.split("/")
                name = fparts[-1]
                #print("The MC sample in question: ",name)
                tf = ROOT.TFile.Open(f)
                numevents = float(str(tf.Get('hnevents').GetString()))
                #print("The number of miniAOD events processed: ",numevents)
                xs = float(xspairs[i][1].split()[0])*1000#Into Femtobarn
                #print("The cross section read from ini file, in femtobarns: ",xs)
                scale = findScale(numevents,xs,lumi)
                scales.append(scale)
                #print("The scaling applied: ",scale)
                h = tf.Get(hname)
                hscaled = h.Clone()
                hscaled.Scale(scale)
                #print("The unscaled contribution: ",h.Integral())
                #print("The   scaled contribution: ",hscaled.Integral())
                hist.Add(hscaled)
                
                #calc hist errors
                df = pd.read_pickle(errs[i])
                sdf = df*scale
                sqrddf = sdf**2
                bkgdfs.append(sqrddf)

                debugstring = name+" "+str(xs)+" "+str(scale)
                #print(debugstring)
                #print("   integral of unscaled histogram: ",h.Integral())
                #print("                            scale: ",scale)
                #print("      intgral of scaled histogram: ",hscaled.Integral())

        uncsqdDYJetsdf = sum(bkgdfs)
        uncDYJetsdf    = uncsqdDYJetsdf**(1/2)

        #alpha = 1.- 0.682689492
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
            #elif hist.GetBinContent(ibin) > 0:
                binerr = uncDYJetsdf[hname][ibin-1]
                hist.SetBinError(ibin,binerr)
            #for garwood interavls    
            #else:
                #binerrbase = ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin))+1,1)-hist.GetBinContent(ibin)
                #errs = [binerrbase*scale for scale in scales]
                #erssq = [err*err for err in errs]
                #sumer = sum(erssq)
                #binerr = (sumer)**(1/2)
                #hist.SetBinError(ibin,binerr)
                

        #print("    integral added hist: ",hist.Integral())
        return hist

    def getAddedHistXSErr(self,hist,samp,region,hname,xserradddir,years = [16,17,18]):
        bkg = self.bkgs[samp]
        xspairs = self.config.items(samp)
        xsuncpairs = self.config.items(samp+'_unc')
        bkgdfs  = []
        scales = []
        for year in years:
            if year == 16:
                lumi = 36.31
            if year == 17:
                lumi = 41.53
            if year == 18:
                lumi = 59.74
                
            files = bkg[year][region][0]
            errs  = bkg[year][region][1]
            if "DYJetsToLL" in samp:
                files.sort(key = orderDY)
                errs.sort(key = orderDY)
            if "TT" in samp:
                files.sort()
                errs.sort()

            #print(files)
            for i,f in enumerate(files):
                fparts = f.split("/")
                name = fparts[-1]
                #print("The MC sample in question: ",name)
                tf = ROOT.TFile.Open(f)
                numevents = float(str(tf.Get('hnevents').GetString()))
                #print("The number of miniAOD events processed: ",numevents)
                #xsunc = float(xsuncpairs[i][1].split()[0])*1000#Into Femtobarn#the uncertainty in the database
                xsunc  = float(xspairs[i][1].split()[0])*1000*0.06#6% uncertainty
                xs = float(xspairs[i][1].split()[0])*1000+xserradddir*xsunc#Into Femtobarn
                #print("The cross section read from ini file, in femtobarns: ",xs)
                scale = findScale(numevents,xs,lumi)
                scales.append(scale)
                #print("The scaling applied: ",scale)
                h = tf.Get(hname)
                hscaled = h.Clone()
                hscaled.Scale(scale)
                #print("The unscaled contribution: ",h.Integral())
                #print("The   scaled contribution: ",hscaled.Integral())
                hist.Add(hscaled)
                
                #calc hist errors
                df = pd.read_pickle(errs[i])
                sdf = df*scale
                sqrddf = sdf**2
                bkgdfs.append(sqrddf)

                debugstring = name+" "+str(xs)+" "+str(scale)
                #print(debugstring)
                #print("   integral of unscaled histogram: ",h.Integral())
                #print("                            scale: ",scale)
                #print("      intgral of scaled histogram: ",hscaled.Integral())

        uncsqdDYJetsdf = sum(bkgdfs)
        uncDYJetsdf    = uncsqdDYJetsdf**(1/2)

        #alpha = 1.- 0.682689492
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
            #elif hist.GetBinContent(ibin) > 0:
                binerr = uncDYJetsdf[hname][ibin-1]
                hist.SetBinError(ibin,binerr)
            #for garwood interavls    
            #else:
                #binerrbase = ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin))+1,1)-hist.GetBinContent(ibin)
                #errs = [binerrbase*scale for scale in scales]
                #erssq = [err*err for err in errs]
                #sumer = sum(erssq)
                #binerr = (sumer)**(1/2)
                #hist.SetBinError(ibin,binerr)
                

        #print("    integral added hist: ",hist.Integral())
        return hist

class run2_local:
    def __init__(self,path,zptcut,hptcut,metcut,btagwp,systr=""):
        self.path = path
        self.zptcut = zptcut
        self.hptcut = hptcut
        self.metcut = metcut
        self.btagwp = btagwp
        self.systr = systr

        #gather data files
        self.run16sb = glob.glob(str(path)+'/Run2016*upout_sideband_'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run17sb = glob.glob(str(path)+'/Run2017*upout_sideband_'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run18sb = glob.glob(str(path)+'/Run2018*upout_sideband_'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run16sr = glob.glob(str(path)+'/Run2016*upout_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run17sr = glob.glob(str(path)+'/Run2017*upout_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run18sr = glob.glob(str(path)+'/Run2018*upout_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run16tr = glob.glob(str(path)+'/Run2016*upout_totalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run17tr = glob.glob(str(path)+'/Run2017*upout_totalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run18tr = glob.glob(str(path)+'/Run2018*upout_totalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')

        #gather errors
        self.run16sberrs = glob.glob(str(path)+'/Run2016*selected_errors_sideband_'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run17sberrs = glob.glob(str(path)+'/Run2017*selected_errors_sideband_'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run18sberrs = glob.glob(str(path)+'/Run2018*selected_errors_sideband_'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run16srerrs = glob.glob(str(path)+'/Run2016*selected_errors_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run17srerrs = glob.glob(str(path)+'/Run2017*selected_errors_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run18srerrs = glob.glob(str(path)+'/Run2018*selected_errors_signalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run16trerrs = glob.glob(str(path)+'/Run2016*selected_errors_totalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run17trerrs = glob.glob(str(path)+'/Run2017*selected_errors_totalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run18trerrs = glob.glob(str(path)+'/Run2018*selected_errors_totalr*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')

        self.data = {18:
                      {"sb":[self.run18sb,self.run18sberrs],
                       "sr":[self.run18sr,self.run18srerrs],
                       "tr":[self.run18tr,self.run18trerrs],
                       },
                      17:
                      {"sb":[self.run17sb,self.run17sberrs],
                       "sr":[self.run17sr,self.run17srerrs],
                       "tr":[self.run17tr,self.run17trerrs],
                       },
                     16:
                      {"sb":[self.run16sb,self.run16sberrs],
                       "sr":[self.run16sr,self.run16srerrs],
                       "tr":[self.run16tr,self.run16trerrs],
                       }
                      }

    def getAddedHist(self,hist,region,hname,years = [16,17,18]):
        data = self.data
        datadfs = []
        for year in years:
            files = data[year][region][0]
            errs  = data[year][region][1]
            files.sort()
            errs.sort()
            for i,f in enumerate(files):
                tf = ROOT.TFile.Open(f)
                h = tf.Get(hname)
                #print(h.Integral())
                hist.Add(h)

                #calc errs
                df = pd.read_pickle(errs[i])
                sqrddf = df**2
                datadfs.append(sqrddf)

        uncsqddf = sum(datadfs)
        uncdf   = uncsqddf**(1/2)
        
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
                binerr = uncdf[hname][ibin-1]
                hist.SetBinError(ibin,binerr)

        return hist

    def getAddedHistPoissonErrors(self,hist,region,hname,years = [16,17,18]):
        data = self.data
        datadfs = []
        alpha = 1.- 0.682689492
        for year in years:
            files = data[year][region][0]
            errs  = data[year][region][1]
            files.sort()
            errs.sort()
            for i,f in enumerate(files):
                tf = ROOT.TFile.Open(f)
                h = tf.Get(hname)
                #print(h.Integral())
                hist.Add(h)
        
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
                binerr = ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin))+1,1)-hist.GetBinContent(ibin)
                #binerrdn = hist.GetBinContent(ibin) - ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin)),1)
                #if int(hist.GetBinContent(ibin)) == 0:
                #    binerrdn = 0
                
                hist.SetBinError(ibin,binerr)

        return hist


class signal_local:
    def __init__(self,path,zptcut,hptcut,metcut,btagwp,xs,years = [16,17,18],systr=""):
        self.path = path
        self.zptcut = zptcut
        self.hptcut = hptcut
        self.metcut = metcut
        self.btagwp = btagwp
        self.systr = systr

        if systr != "":
            systr = "_"+systr

        lumidict = {16:36.31,17:41.53,18:59.74}
        lumi = 0
        for year in years:
            lumi += lumidict[year]

        
        
        #gather signas
        self.sigsr = glob.glob(str(path)+'/Autumn18*Zp*_upout_signalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.sigsb = glob.glob(str(path)+'/Autumn18*Zp*_upout_sideband'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.sigtr = glob.glob(str(path)+'/Autumn18*Zp*_upout_totalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        #gather signas
        self.sigsrerrs = glob.glob(str(path)+'/Autumn18*Zp*_selected_errors_signalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.sigsberrs = glob.glob(str(path)+'/Autumn18*Zp*_selected_errors_sideband'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.sigtrerrs = glob.glob(str(path)+'/Autumn18*Zp*_selected_errors_totalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')

        self.sigsrerrsdict = makeSignalErrorDictionary(self.sigsrerrs)
        self.sigsberrsdict = makeSignalErrorDictionary(self.sigsberrs)
        self.sigtrerrsdict = makeSignalErrorDictionary(self.sigtrerrs)

            

        self.sigs ={"sb":[self.sigsb,self.sigsberrsdict],
                    "sr":[self.sigsr,self.sigsrerrsdict],
                    "tr":[self.sigtr,self.sigtrerrsdict]
                    }
        
        #sig_colors = colsFromPalette(self.sigsr,ROOT.kCMYK)
        
        #prep signals
        #self.prepsigsr = prepSig(self.sigsr,sig_colors,xs,lumi)
        
    def getPreppedSig(self,region,xs,years=[16,17,18]):
        signal = self.sigs[region][0]
        #print(signal)
        signalerrs = self.sigs[region][1]
        sig_colors = colsFromPalette(signal,ROOT.kCMYK)
        lumidict = {16:36.31,17:41.53,18:59.74}
        l = 0
        for year in years:
            l += lumidict[year]

        #print("          for signal you are using a lumi of ",l)
        #print("          for signal you are using the region ",region)
    
        preppedsig = prepSig(signal,sig_colors,xs,l,signalerrs)
        return preppedsig

class signal_run2_local:
    def __init__(self,path,zptcut,hptcut,metcut,btagwp,xs,years = [16,17,18],systr=""):
        self.path = path
        self.zptcut = zptcut
        self.hptcut = hptcut
        self.metcut = metcut
        self.btagwp = btagwp
        self.systr = systr
        
        if systr != "":
            systr = "_"+systr
            
        lumidict = {16:36.31,17:41.53,18:59.74}
        lumi = 0
        for year in years:
            lumi += lumidict[year]
            
        #gather signals
        #2016
        self.sig16sr = sorted(glob.glob(str(path)+'/Summer16.Zp*_upout_signalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root'))
        self.sig16sb = sorted(glob.glob(str(path)+'/Summer16.Zp*_upout_sideband'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root'))
        self.sig16tr = sorted(glob.glob(str(path)+'/Summer16.Zp*_upout_totalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root'))
        #2017
        self.sig17sr = sorted(glob.glob(str(path)+'/Fall17.Zp*_upout_signalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root'))
        self.sig17sb = sorted(glob.glob(str(path)+'/Fa;;17.Zp*_upout_sideband'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root'))
        self.sig17tr = sorted(glob.glob(str(path)+'/Fall17.Zp*_upout_totalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root'))
        #2018
        self.sig18sr = sorted(glob.glob(str(path)+'/Autumn18.Zp*_upout_signalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root'))
        self.sig18sb = sorted(glob.glob(str(path)+'/Autumn18.Zp*_upout_sideband'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root'))
        self.sig18tr = sorted(glob.glob(str(path)+'/Autunmn18.Zp*_upout_totalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root'))
        
        #gather signal errors
        #2016
        self.sig16srerrs = sorted(glob.glob(str(path)+'/Summer16.Zp*_selected_errors_signalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl'))
        self.sig16sberrs = sorted(glob.glob(str(path)+'/Summer16.Zp*_selected_errors_sideband'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl'))
        self.sig16trerrs = sorted(glob.glob(str(path)+'/Summer16.Zp*_selected_errors_totalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl'))
        #2017
        self.sig17srerrs = sorted(glob.glob(str(path)+'/Fall17.Zp*_selected_errors_signalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl'))
        self.sig17sberrs = sorted(glob.glob(str(path)+'/Fall17.Zp*_selected_errors_sideband'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl'))
        self.sig17trerrs = sorted(glob.glob(str(path)+'/Fall17.Zp*_selected_errors_totalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl'))
        #2018
        self.sig18srerrs = sorted(glob.glob(str(path)+'/Autumn18.Zp*_selected_errors_signalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl'))
        self.sig18sberrs = sorted(glob.glob(str(path)+'/Autumn18.Zp*_selected_errors_sideband'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl'))
        self.sig18trerrs = sorted(glob.glob(str(path)+'/Autumn18.Zp*_selected_errors_totalr'+systr+'*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl'))


        self.signames18 = [x.split('ZpAnomalonHZ_UFO-')[-1].split("_Tune")[0] for x in self.sig18sr]
        self.signames17 = [x.split('ZpAnomalonHZ_UFO-')[-1].split("_Tune")[0] for x in self.sig17sr]
        self.signames16 = [x.split('ZpAnomalonHZ_UFO-')[-1].split("_Tune")[0] for x in self.sig16sr]

        int1817 = set(self.signames18).intersection(set(self.signames17))
        int1816 = set(self.signames18).intersection(set(self.signames16))
        not1718 = set(self.signames18) - int1817
        not1618 = set(self.signames18) - int1816
        neither1617 = not1718.union(not1618)
        allthreeyearssig = list(set(self.signames18) - neither1617)
        #print(self.signames18)
        #print(self.signames17)
        #print(self.signames16)
        #print(allthreeyearssig)

        self.sigs ={18:
                    {"sb":[self.sig18sb,self.sig18sberrs],
                     "sr":[self.sig18sr,self.sig18srerrs],
                     "tr":[self.sig18tr,self.sig18trerrs]},
                    17:
                    {"sb":[self.sig17sb,self.sig17sberrs],
                     "sr":[self.sig17sr,self.sig17srerrs],
                     "tr":[self.sig17tr,self.sig17trerrs]},
                    16:
                    {"sb":[self.sig16sb,self.sig16sberrs],
                     "sr":[self.sig16sr,self.sig16srerrs],
                     "tr":[self.sig16tr,self.sig16trerrs]},
        }
        
        #removing names with no match in 2018
        for name in neither1617:
            badidx = self.signames18.index(name)
            files2018 = [self.sigs[18]["sb"],self.sigs[18]["sr"],self.sigs[18]["tr"]]
            for region in files2018:
                for group in region:
                    if len(group)  > badidx:
                        badval = group[badidx]
                        if name in badval:
                            del group[badidx]

        self.sigsbyname = {}
        for i,name in enumerate(sorted(allthreeyearssig)):
            #print(i)
            #print(name)
            self.sigsbyname[name] = {'sb':
                                     {18:[[],[]],
                                      17:[[],[]],
                                      16:[[],[]]
                                      },
                                     'sr':
                                     {18:[[],[]],
                                      17:[[],[]],
                                      16:[[],[]]
                                      },
                                     'tr':
                                     {18:[[],[]],
                                      17:[[],[]],
                                      16:[[],[]]
                                      },
            }
            for year in lumidict.keys():
                #print(year)
                for region in self.sigs[year].keys():
                    #print(region)
                    for group in self.sigs[year][region]:
                        #print(group)
                        groupidx = self.sigs[year][region].index(group)
                        if len(group) > 0:
                            for fstr in group:
                                if name in fstr:
                                    self.sigsbyname[name][region][year][groupidx].append(fstr)
                        #    print(group[i])
                        #    if name in group[i]:
                        #        self.sigsbyname[name][region][year][groupidx].append(group[i])
    #DIeally that whole shenaniga would hae just been a dataframe...
                        

        
    def getAddedHist(self,signalname,xs,hist,region,hname,years = [16,17,18]):
        sig = self.sigsbyname[signalname]
        #print(sig)
        sigdfs  = []
        scales = []
        for year in years:
            if year == 16:
                lumi = 36.31
            if year == 17:
                lumi = 41.53
            if year == 18:
                lumi = 59.74
                
            f = sig[region][year][0][0]
            errs  = sig[region][year][1][0]
            fparts = f.split("/")
            name = fparts[-1]
            tf = ROOT.TFile.Open(f)
            numevents = float(str(tf.Get('hnevents').GetString()))
            scale = findScale(numevents,xs,lumi)
            scales.append(scale)
            h = tf.Get(hname)
            hscaled = h.Clone()
            hscaled.Scale(scale)
            hist.Add(hscaled)
                
            #calc hist errors
            df = pd.read_pickle(errs)
            sdf = df*scale
            sqrddf = sdf**2
            sigdfs.append(sqrddf)

            debugstring = name+" "+str(xs)+" "+str(scale)
            #print(debugstring)
            #print("   integral of unscaled histogram: ",h.Integral())
            #print("                            scale: ",scale)
            #print("      intgral of scaled histogram: ",hscaled.Integral())

        uncsqdDYJetsdf = sum(sigdfs)
        uncDYJetsdf    = uncsqdDYJetsdf**(1/2)

        #alpha = 1.- 0.682689492
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
            #elif hist.GetBinContent(ibin) > 0:
                binerr = uncDYJetsdf[hname][ibin-1]
                hist.SetBinError(ibin,binerr)
            #for garwood interavls    
            #else:
                #binerrbase = ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin))+1,1)-hist.GetBinContent(ibin)
                #errs = [binerrbase*scale for scale in scales]
                #erssq = [err*err for err in errs]
                #sumer = sum(erssq)
                #binerr = (sumer)**(1/2)
                #hist.SetBinError(ibin,binerr)
                

        #print("    integral added hist: ",hist.Integral())
        return hist




        

    

class validation_local:
    def __init__(self,path,zptcut,hptcut,metcut,btagwp,systr):
        self.path = path
        self.zptcut = zptcut
        self.hptcut = hptcut
        self.metcut = metcut
        self.btagwp = btagwp
        self.systr = systr

        # background MC files
        #2017
        self.f17dyjetsb = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17dyjetsr = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsb = glob.glob(str(path)+'/Fall17.TTT*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsr = glob.glob(str(path)+'/Fall17.TTT*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsb = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsr = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsb = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsr = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')

        #2018
        self.a18dyjetsb = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18dyjetsr = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsb = glob.glob(str(path)+'/Autumn18.TTT*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsr = glob.glob(str(path)+'/Autumn18.TTT*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsb = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsr = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsb = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsr = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')

        #2016
        self.s16dyjetsb = glob.glob(str(path)+'/Summer16v3.DYJetsToLL_M-50_HT*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16dyjetsr = glob.glob(str(path)+'/Summer16v3.DYJetsToLL_M-50_HT*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16ttsb = glob.glob(str(path)+'/Summer16v3.TTT*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16ttsr = glob.glob(str(path)+'/Summer16v3.TTT*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16wzsb = glob.glob(str(path)+'/Summer16v3.WZTo2L2Q*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16wzsr = glob.glob(str(path)+'/Summer16v3.WZTo2L2Q*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16zzsb = glob.glob(str(path)+'/Summer16v3.ZZTo2L2Q*_upout_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16zzsr = glob.glob(str(path)+'/Summer16v3.ZZTo2L2Q*_upout_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')


        #gather errors
        self.f17dyjetsberrs = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17dyjetsrerrs = glob.glob(str(path)+'/Fall17.DYJetsToLL_M-50_HT*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsberrs = glob.glob(str(path)+'/Fall17.TTT*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17ttsrerrs = glob.glob(str(path)+'/Fall17.TTT*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsberrs = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17wzsrerrs = glob.glob(str(path)+'/Fall17.WZTo2L2Q*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsberrs = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.f17zzsrerrs = glob.glob(str(path)+'/Fall17.ZZTo2L2Q*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')

        self.a18dyjetsberrs = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18dyjetsrerrs = glob.glob(str(path)+'/Autumn18.DYJetsToLL_M-50_HT*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsberrs = glob.glob(str(path)+'/Autumn18.TTT*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18ttsrerrs = glob.glob(str(path)+'/Autumn18.TTT*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsberrs = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18wzsrerrs = glob.glob(str(path)+'/Autumn18.WZTo2L2Q*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsberrs = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.a18zzsrerrs = glob.glob(str(path)+'/Autumn18.ZZTo2L2Q*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')

        self.s16dyjetsberrs = glob.glob(str(path)+'/Summer16v3.DYJetsToLL_M-50_HT*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16dyjetsrerrs = glob.glob(str(path)+'/Summer16v3.DYJetsToLL_M-50_HT*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16ttsberrs = glob.glob(str(path)+'/Summer16v3.TTT*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16ttsrerrs = glob.glob(str(path)+'/Summer16v3.TTT*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16wzsberrs = glob.glob(str(path)+'/Summer16v3.WZTo2L2Q*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16wzsrerrs = glob.glob(str(path)+'/Summer16v3.WZTo2L2Q*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16zzsberrs = glob.glob(str(path)+'/Summer16v3.ZZTo2L2Q*_selected_errors_sideband_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')
        self.s16zzsrerrs = glob.glob(str(path)+'/Summer16v3.ZZTo2L2Q*_selected_errors_validationr_'+systr+'_DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'*')


        self.bkgs = {"DYJetsToLL":
                {18:
                 {"sb":[self.a18dyjetsb,self.a18dyjetsberrs],
                  "sr":[self.a18dyjetsr,self.a18dyjetsrerrs]},
                 17:
                 {"sb":[self.f17dyjetsb,self.f17dyjetsberrs],
                  "sr":[self.f17dyjetsr,self.f17dyjetsrerrs]},
                 16:
                 {"sb":[self.s16dyjetsb,self.s16dyjetsberrs],
                  "sr":[self.s16dyjetsr,self.s16dyjetsrerrs]},
                },
                "TT":
                {18:
                 {"sb":[self.a18ttsb,self.a18ttsberrs],
                  "sr":[self.a18ttsr,self.a18ttsrerrs]},
                 17:
                 {"sb":[self.f17ttsb,self.f17ttsberrs],
                  "sr":[self.f17ttsr,self.f17ttsrerrs]},
                 16:
                 {"sb":[self.s16ttsb,self.s16ttsberrs],
                  "sr":[self.s16ttsr,self.s16ttsrerrs]},
                },
                "WZTo2L2Q":
                {18:
                 {"sb":[self.a18wzsb,self.a18wzsberrs],
                  "sr":[self.a18wzsr,self.a18wzsrerrs]},
                 17:
                 {"sb":[self.f17wzsb,self.f17wzsberrs],
                  "sr":[self.f17wzsr,self.f17wzsrerrs]},
                 16:
                 {"sb":[self.s16wzsb,self.s16wzsberrs],
                  "sr":[self.s16wzsr,self.s16wzsrerrs]},
                },
                "ZZTo2L2Q":
                {18:
                 {"sb":[self.a18zzsb,self.a18zzsberrs],
                  "sr":[self.a18zzsr,self.a18zzsrerrs]},
                 17:
                 {"sb":[self.f17zzsb,self.f17zzsberrs],
                  "sr":[self.f17zzsr,self.f17zzsrerrs]},
                 16:
                 {"sb":[self.s16zzsb,self.s16zzsberrs],
                  "sr":[self.s16zzsr,self.s16zzsrerrs]},
                }
        }

        #gather data files
        self.run16sb = glob.glob(str(path)+'/Run2016*upout_sideband_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run17sb = glob.glob(str(path)+'/Run2017*upout_sideband_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run18sb = glob.glob(str(path)+'/Run2018*upout_sideband_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run16vr = glob.glob(str(path)+'/Run2016*upout_validationr_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run17vr = glob.glob(str(path)+'/Run2017*upout_validationr_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        self.run18vr = glob.glob(str(path)+'/Run2018*upout_validationr_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.root')
        
        #gather errors
        self.run16sberrs = glob.glob(str(path)+'/Run2016*selected_errors_sideband_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run17sberrs = glob.glob(str(path)+'/Run2017*selected_errors_sideband_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run18sberrs = glob.glob(str(path)+'/Run2018*selected_errors_sideband_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run16vrerrs = glob.glob(str(path)+'/Run2016*selected_errors_validationr_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run17vrerrs = glob.glob(str(path)+'/Run2017*selected_errors_validationr_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        self.run18vrerrs = glob.glob(str(path)+'/Run2018*selected_errors_validationr_*_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)+'.pkl')
        
        self.data = {18:
                     {"sb":[self.run18sb,self.run18sberrs],
                      "vr":[self.run18vr,self.run18vrerrs],
                     },
                     17:
                     {"sb":[self.run17sb,self.run17sberrs],
                      "vr":[self.run17vr,self.run17vrerrs],
                     },
                     16:
                     {"sb":[self.run16sb,self.run16sberrs],
                      "vr":[self.run16vr,self.run16vrerrs],
                     }
        }


        self.config = configparser.RawConfigParser()
        self.config.optionxform = str
        fp = open('xsects_2017.ini')#2017 and 2018dy+jets xs same
        self.config.read_file(fp)

        
    def getAddedHist(self,hist,samp,region,hname,years = [16,17,18]):
        bkg = self.bkgs[samp]
        xspairs = self.config.items(samp)
        bkgdfs  = []
        
        for year in years:
            if year == 16:
                lumi = 36.31
            if year == 17:
                lumi = 41.53
            if year == 18:
                lumi = 59.74
                
            files = bkg[year][region][0]
            errs  = bkg[year][region][1]
            if "DYJetsToLL" in samp:
                files.sort(key = orderDY)
                errs.sort(key = orderDY)
            if "TT" in samp:
                files.sort()
                errs.sort()
            for i,f in enumerate(files):
                fparts = f.split("/")
                name = fparts[-1]
                tf = ROOT.TFile.Open(f)
                numevents = float(str(tf.Get('hnevents').GetString()))
                xs = float(xspairs[i][1].split()[0])*1000#Into Femtobarn
                scale = findScale(numevents,xs,lumi)
                h = tf.Get(hname)
                h.Scale(scale)
                hist.Add(h)
                
                #calc hist errors
                df = pd.read_pickle(errs[i])
                sdf = df*scale
                sqrddf = sdf**2
                bkgdfs.append(sqrddf)
                
            debugstring = name+" "+str(xs)+" "+str(scale)
                
            uncsqdDYJetsdf = sum(bkgdfs)
            uncDYJetsdf    = uncsqdDYJetsdf**(1/2)
        
            for ibin in range(hist.GetNbinsX()+1):
                if ibin == 0:
                    continue
                else:
                    binerr = uncDYJetsdf[hname][ibin-1]
                    hist.SetBinError(ibin,binerr)

        return hist

    def getAddedHistData(self,hist,region,hname,years = [16,17,18]):
        data = self.data
        datadfs = []
        
        for year in years:
            files = data[year][region][0]
            errs  = data[year][region][1]
            files.sort()
            errs.sort()
            for i,f in enumerate(files):
                tf = ROOT.TFile.Open(f)
                h = tf.Get(hname)
                hist.Add(h)
                
                #calc errs
                df = pd.read_pickle(errs[i])
                sqrddf = df**2
                datadfs.append(sqrddf)
                
        uncsqddf = sum(datadfs)
        uncdf   = uncsqddf**(1/2)
            
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
                binerr = uncdf[hname][ibin-1]
                hist.SetBinError(ibin,binerr)
                
        return hist

sampnames_dy2018 = ['Autumn18.DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights',
                    'Autumn18.DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights',
                    'Autumn18.DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights',
                    'Autumn18.DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights',
                    'Autumn18.DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights',
                    'Autumn18.DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights',
                    'Autumn18.DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights'
                ]
sampnames_tt2018 = ['Autumn18.TTTo2L2Nu_TuneCP5',
                   'Autumn18.TTToHadronic_TuneCP5',
                   'Autumn18.TTToSemiLeptonic_TuneCP5'
               ]

sampnames_wz2018 = ['Autumn18.WZTo2L2Q']

sampnames_zz2018 = ['Autumn18.ZZTo2L2Q']

sampnames_dy2017 = ['Fall17.DYJetsToLL_M-50_HT-100to200_TuneCP5',
                    'Fall17.DYJetsToLL_M-50_HT-1200to2500_TuneCP5',
                    'Fall17.DYJetsToLL_M-50_HT-200to400_TuneCP5',
                    'Fall17.DYJetsToLL_M-50_HT-2500toInf_TuneCP5',
                    'Fall17.DYJetsToLL_M-50_HT-400to600_TuneCP5',
                    'Fall17.DYJetsToLL_M-50_HT-600to800_TuneCP5',
                    'Fall17.DYJetsToLL_M-50_HT-800to1200_TuneCP5'
                ]

sampnames_tt2017 = ['Fall17.TTTo2L2Nu_TuneCP5_PSweights',
                    'Fall17.TTToHadronic_TuneCP5',
                    'Fall17.TTToSemiLeptonic_TuneCP5',
                ]

sampnames_wz2017 = ['Fall17.WZTo2L2Q']

sampnames_zz2017 = ['Fall17.ZZTo2L2Q']

sampnames_dy2016 = ['Summer16v3.DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1',
                    'Summer16v3.DYJetsToLL_M-50_HT-1200to2500_TuneCUETP8M1',
                    'Summer16v3.DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1',
                    'Summer16v3.DYJetsToLL_M-50_HT-2500toInf_TuneCUETP8M1',
                    'Summer16v3.DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1',
                    'Summer16v3.DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1',
                    'Summer16v3.DYJetsToLL_M-50_HT-800to1200_TuneCUETP8M1'
                ]

sampnames_tt2016 = ['Summer16v3.TTTo2L2Nu_TuneCP5_PSweights',
                    'Summer16v3.TTToHadronic_TuneCP5_PSweights',
                    'Summer16v3.TTToSemiLeptonic_TuneCP5_PSweights'
                ]

sampnames_wz2016 = ['Summer16v3.WZTo2L2Q']

sampnames_zz2016 = ['Summer16v3.ZZTo2L2Q']


class backgrounds:
    def __init__(self,path,zptcut,hptcut,metcut,btagwp,systr=""):
        self.path = path
        self.zptcut = zptcut
        self.hptcut = hptcut
        self.metcut = metcut
        self.btagwp = btagwp
        self.systr = systr
        self.ftype   = 'upout'
        self.cutstring = 'DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)

        #gather background MC files
        #2017
        self.f17dyjetsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_dy2017]
        self.f17dyjetsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_dy2017]
        self.f17dyjettr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_dy2017]
        self.f17ttsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_tt2017]
        self.f17ttsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_tt2017]
        self.f17tttr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_tt2017]
        self.f17wzsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_wz2017]
        self.f17wzsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_wz2017]
        self.f17wztr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_wz2017]
        self.f17zzsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_zz2017]
        self.f17zzsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_zz2017]
        self.f17zztr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_zz2017]

        #2018
        self.a18dyjetsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_dy2018]
        self.a18dyjetsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_dy2018]
        self.a18dyjettr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_dy2018]
        self.a18ttsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_tt2018]
        self.a18ttsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_tt2018]
        self.a18tttr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_tt2018]
        self.a18wzsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_wz2018]
        self.a18wzsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_wz2018]
        self.a18wztr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_wz2018]
        self.a18zzsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_zz2018]
        self.a18zzsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_zz2018]
        self.a18zztr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_zz2018]

        #2016
        self.s16dyjetsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_dy2016]
        self.s16dyjetsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_dy2016]
        self.s16dyjettr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_dy2016]
        self.s16ttsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_tt2016]
        self.s16ttsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_tt2016]
        self.s16tttr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_tt2016]
        self.s16wzsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_wz2016]
        self.s16wzsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_wz2016]
        self.s16wztr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_wz2016]
        self.s16zzsb = [path+s+"_"+self.ftype+"_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_zz2016]
        self.s16zzsr = [path+s+"_"+self.ftype+"_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_zz2016]
        self.s16zztr = [path+s+"_"+self.ftype+"_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_zz2016]

        #gather errors
        self.f17dyjetsberrs = [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_dy2017]
        self.f17dyjetsrerrs = [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_dy2017]
        self.f17dyjettrerrs = [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_dy2017]
        self.f17ttsberrs =    [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_tt2017]
        self.f17ttsrerrs =    [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_tt2017]
        self.f17tttrerrs =    [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_tt2017]
        self.f17wzsberrs =    [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_wz2017]
        self.f17wzsrerrs =    [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_wz2017]
        self.f17wztrerrs =    [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_wz2017]
        self.f17zzsberrs =    [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_zz2017]
        self.f17zzsrerrs =    [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_zz2017]
        self.f17zztrerrs =    [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_zz2017]

        self.a18dyjetsberrs = [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_dy2018]
        self.a18dyjetsrerrs = [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_dy2018]
        self.a18dyjettrerrs = [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_dy2018]
        self.a18ttsberrs =    [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_tt2018]
        self.a18ttsrerrs =    [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_tt2018]
        self.a18tttrerrs =    [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_tt2018]
        self.a18wzsberrs =    [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_wz2018]
        self.a18wzsrerrs =    [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_wz2018]
        self.a18wztrerrs =    [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_wz2018]
        self.a18zzsberrs =    [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_zz2018]
        self.a18zzsrerrs =    [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_zz2018]
        self.a18zztrerrs =    [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_zz2018]  
        
        self.s16dyjetsberrs = [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_dy2016]
        self.s16dyjetsrerrs = [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_dy2016] 
        self.s16dyjettrerrs = [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_dy2016]  
        self.s16ttsberrs =    [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_tt2016]
        self.s16ttsrerrs =    [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_tt2016] 
        self.s16tttrerrs =    [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_tt2016]  
        self.s16wzsberrs =    [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_wz2016]
        self.s16wzsrerrs =    [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_wz2016] 
        self.s16wztrerrs =    [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_wz2016]  
        self.s16zzsberrs =    [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_zz2016]
        self.s16zzsrerrs =    [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_zz2016] 
        self.s16zztrerrs =    [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_zz2016]  


        self.bkgs = {"DYJetsToLL":
                {18:
                 {"sb":[self.a18dyjetsb,self.a18dyjetsberrs],
                  "sr":[self.a18dyjetsr,self.a18dyjetsrerrs],
                  "tr":[self.a18dyjettr,self.a18dyjettrerrs]},
                 17:
                 {"sb":[self.f17dyjetsb,self.f17dyjetsberrs],
                  "sr":[self.f17dyjetsr,self.f17dyjetsrerrs],
                  "tr":[self.f17dyjettr,self.f17dyjettrerrs]},
                 16:
                 {"sb":[self.s16dyjetsb,self.s16dyjetsberrs],
                  "sr":[self.s16dyjetsr,self.s16dyjetsrerrs],
                  "tr":[self.s16dyjettr,self.s16dyjettrerrs]},
                },
                "TT":
                {18:
                 {"sb":[self.a18ttsb,self.a18ttsberrs],
                  "sr":[self.a18ttsr,self.a18ttsrerrs],
                  "tr":[self.a18tttr,self.a18tttrerrs]},
                 17:
                 {"sb":[self.f17ttsb,self.f17ttsberrs],
                  "sr":[self.f17ttsr,self.f17ttsrerrs],
                  "tr":[self.f17tttr,self.f17tttrerrs]},
                 16:
                 {"sb":[self.s16ttsb,self.s16ttsberrs],
                  "sr":[self.s16ttsr,self.s16ttsrerrs],
                  "tr":[self.s16tttr,self.s16tttrerrs]},
                },
                "WZTo2L2Q":
                {18:
                 {"sb":[self.a18wzsb,self.a18wzsberrs],
                  "sr":[self.a18wzsr,self.a18wzsrerrs],
                  "tr":[self.a18wztr,self.a18wztrerrs]},
                 17:
                 {"sb":[self.f17wzsb,self.f17wzsberrs],
                  "sr":[self.f17wzsr,self.f17wzsrerrs],
                  "tr":[self.f17wztr,self.f17wztrerrs]},
                 16:
                 {"sb":[self.s16wzsb,self.s16wzsberrs],
                  "sr":[self.s16wzsr,self.s16wzsrerrs],
                  "tr":[self.s16wztr,self.s16wztrerrs]},
                },
                "ZZTo2L2Q":
                {18:
                 {"sb":[self.a18zzsb,self.a18zzsberrs],
                  "sr":[self.a18zzsr,self.a18zzsrerrs],
                  "tr":[self.a18zztr,self.a18zztrerrs]},
                 17:
                 {"sb":[self.f17zzsb,self.f17zzsberrs],
                  "sr":[self.f17zzsr,self.f17zzsrerrs],
                  "tr":[self.f17zztr,self.f17zztrerrs]},
                 16:
                 {"sb":[self.s16zzsb,self.s16zzsberrs],
                  "sr":[self.s16zzsr,self.s16zzsrerrs],
                  "tr":[self.s16zztr,self.s16zztrerrs]},
                }
        }

        self.config = configparser.RawConfigParser()
        self.config.optionxform = str
        fp = open('xsects_2017.ini')#2017 and 2018dy+jets xs same
        self.config.read_file(fp)

        
    def getAddedHist(self,hist,samp,region,hname,years = [16,17,18]):
        bkg = self.bkgs[samp]
        xspairs = self.config.items(samp)
        #print(xspairs)
        bkgdfs  = []
        scales = []
        for year in years:
            if year == 16:
                lumi = 36.31
            if year == 17:
                lumi = 41.53
            if year == 18:
                lumi = 59.74
                
            files = bkg[year][region][0]
            errs  = bkg[year][region][1]
            if "DYJetsToLL" in samp:
                files.sort(key = orderDY)
                errs.sort(key = orderDY)
            if "TT" in samp:
                files.sort()
                errs.sort()

            #print(files)
            for i,f in enumerate(files):
                fparts = f.split("/")
                name = fparts[-1]
                #print("The MC sample in question: ",name)
                tf = ROOT.TFile.Open(f)
                numevents = float(str(tf.Get('hnevents').GetString()))
                #print("The number of miniAOD events processed: ",numevents)
                xs = float(xspairs[i][1].split()[0])*1000#Into Femtobarn
                #print("The cross section read from ini file, in femtobarns: ",xs)
                scale = findScale(numevents,xs,lumi)
                scales.append(scale)
                #print("The scaling applied: ",scale)
                h = tf.Get(hname)
                h.SetDirectory(0)
                hscaled = h.Clone()
                hscaled.Scale(scale)
                #print("The unscaled contribution: ",h.Integral())
                #print("The   scaled contribution: ",hscaled.Integral())
                hist.Add(hscaled)
                tf.Close()
                
                #calc hist errors
                df = pd.read_pickle(errs[i])
                sdf = df*scale
                sqrddf = sdf**2
                bkgdfs.append(sqrddf)

                debugstring = name+" "+str(xs)+" "+str(scale)
                #print(debugstring)
                #print("   integral of unscaled histogram: ",h.Integral())
                #print("                            scale: ",scale)
                #print("      intgral of scaled histogram: ",hscaled.Integral())

        uncsqdDYJetsdf = sum(bkgdfs)
        uncDYJetsdf    = uncsqdDYJetsdf**(1/2)

        #alpha = 1.- 0.682689492
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
            #elif hist.GetBinContent(ibin) > 0:
                binerr = uncDYJetsdf[hname][ibin-1]
                hist.SetBinError(ibin,binerr)
            #for garwood interavls    
            #else:
                #binerrbase = ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin))+1,1)-hist.GetBinContent(ibin)
                #errs = [binerrbase*scale for scale in scales]
                #erssq = [err*err for err in errs]
                #sumer = sum(erssq)
                #binerr = (sumer)**(1/2)
                #hist.SetBinError(ibin,binerr)
                

        #print("    integral added hist: ",hist.Integral())
        return hist

    def getAddedHistXSErr(self,hist,samp,region,hname,xserradddir,years = [16,17,18]):
        bkg = self.bkgs[samp]
        xspairs = self.config.items(samp)
        xsuncpairs = self.config.items(samp+'_unc')
        bkgdfs  = []
        scales = []
        for year in years:
            if year == 16:
                lumi = 36.31
            if year == 17:
                lumi = 41.53
            if year == 18:
                lumi = 59.74
                
            files = bkg[year][region][0]
            errs  = bkg[year][region][1]
            if "DYJetsToLL" in samp:
                files.sort(key = orderDY)
                errs.sort(key = orderDY)
            if "TT" in samp:
                files.sort()
                errs.sort()

            #print(files)
            for i,f in enumerate(files):
                fparts = f.split("/")
                name = fparts[-1]
                #print("The MC sample in question: ",name)
                tf = ROOT.TFile.Open(f)
                numevents = float(str(tf.Get('hnevents').GetString()))
                #print("The number of miniAOD events processed: ",numevents)
                #xsunc = float(xsuncpairs[i][1].split()[0])*1000#Into Femtobarn#the uncertainty in the database
                xsunc  = float(xspairs[i][1].split()[0])*1000*0.06#6% uncertainty
                xs = float(xspairs[i][1].split()[0])*1000+xserradddir*xsunc#Into Femtobarn
                #print("The cross section read from ini file, in femtobarns: ",xs)
                scale = findScale(numevents,xs,lumi)
                scales.append(scale)
                #print("The scaling applied: ",scale)
                h = tf.Get(hname)
                hscaled = h.Clone()
                hscaled.Scale(scale)
                #print("The unscaled contribution: ",h.Integral())
                #print("The   scaled contribution: ",hscaled.Integral())
                hist.Add(hscaled)
                
                #calc hist errors
                df = pd.read_pickle(errs[i])
                sdf = df*scale
                sqrddf = sdf**2
                bkgdfs.append(sqrddf)

                debugstring = name+" "+str(xs)+" "+str(scale)
                #print(debugstring)
                #print("   integral of unscaled histogram: ",h.Integral())
                #print("                            scale: ",scale)
                #print("      intgral of scaled histogram: ",hscaled.Integral())

        uncsqdDYJetsdf = sum(bkgdfs)
        uncDYJetsdf    = uncsqdDYJetsdf**(1/2)

        #alpha = 1.- 0.682689492
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
            #elif hist.GetBinContent(ibin) > 0:
                binerr = uncDYJetsdf[hname][ibin-1]
                hist.SetBinError(ibin,binerr)
            #for garwood interavls    
            #else:
                #binerrbase = ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin))+1,1)-hist.GetBinContent(ibin)
                #errs = [binerrbase*scale for scale in scales]
                #erssq = [err*err for err in errs]
                #sumer = sum(erssq)
                #binerr = (sumer)**(1/2)
                #hist.SetBinError(ibin,binerr)
                

        #print("    integral added hist: ",hist.Integral())
        return hist

sampnames_run2018 =['Run2018A-17Sep2018-v1.SingleMuon',
                    'Run2018B-17Sep2018-v1.SingleMuon',
                    'Run2018C-17Sep2018-v1.SingleMuon',
                    'Run2018D-22Jan2019-v2.SingleMuon'
                ]
sampnames_run2017 =['Run2017B-31Mar2018-v1.SingleMuon',
                    'Run2017C-31Mar2018-v1.SingleMuon',
                    'Run2017D-31Mar2018-v1.SingleMuon',
                    'Run2017E-31Mar2018-v1.SingleMuon',
                    'Run2017F-31Mar2018-v1.SingleMuon'
                ]
sampnames_run2016 =['Run2016B-17Jul2018-ver2-v1.SingleMuons',
                    'Run2016C-17Jul2018-v1.SingleMuon',
                    'Run2016D-17Jul2018-v1.SingleMuon',
                    'Run2016E-17Jul2018-v1.SingleMuon',
                    'Run2016F-17Jul2018-v1.SingleMuon',
                    'Run2016G-17Jul2018-v1.SingleMuon'
                ]


class run2:
    def __init__(self,path,zptcut,hptcut,metcut,btagwp,systr=""):
        self.path = path
        self.zptcut = zptcut
        self.hptcut = hptcut
        self.metcut = metcut
        self.btagwp = btagwp
        self.systr = systr
        self.cutstring = 'DeepMassDecorrelTagZHbbvsQCD_Zptcut'+str(zptcut)+'_Hptcut'+str(hptcut)+'_metcut'+str(metcut)+'_btagwp'+str(btagwp)

        #gather data files
        self.run16sb = [path+s+"_upout_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_run2016]
        self.run17sb = [path+s+"_upout_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_run2017]
        self.run18sb = [path+s+"_upout_sideband_"+systr+"_"+self.cutstring+".root" for s in sampnames_run2018]
        self.run16sr = [path+s+"_upout_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_run2016]
        self.run17sr = [path+s+"_upout_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_run2017]
        self.run18sr = [path+s+"_upout_signalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_run2018]
        self.run16tr = [path+s+"_upout_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_run2016]
        self.run17tr = [path+s+"_upout_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_run2017]
        self.run18tr = [path+s+"_upout_totalr_"+systr+"_"+self.cutstring+".root" for s in sampnames_run2018]

        #gather errors
        self.run16sberrs = [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_run2016]
        self.run17sberrs = [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_run2017]
        self.run18sberrs = [path+s+"_selected_errors_sideband_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_run2018]
        self.run16srerrs = [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_run2016]
        self.run17srerrs = [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_run2017]
        self.run18srerrs = [path+s+"_selected_errors_signalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_run2018]
        self.run16trerrs = [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_run2016]
        self.run17trerrs = [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_run2017]
        self.run18trerrs = [path+s+"_selected_errors_totalr_"+systr+"_"+self.cutstring+".pkl" for s in sampnames_run2018]

        self.data = {18:
                      {"sb":[self.run18sb,self.run18sberrs],
                       "sr":[self.run18sr,self.run18srerrs],
                       "tr":[self.run18tr,self.run18trerrs],
                       },
                      17:
                      {"sb":[self.run17sb,self.run17sberrs],
                       "sr":[self.run17sr,self.run17srerrs],
                       "tr":[self.run17tr,self.run17trerrs],
                       },
                     16:
                      {"sb":[self.run16sb,self.run16sberrs],
                       "sr":[self.run16sr,self.run16srerrs],
                       "tr":[self.run16tr,self.run16trerrs],
                       }
                      }

    def getAddedHist(self,hist,region,hname,years = [16,17,18]):
        data = self.data
        datadfs = []
        for year in years:
            files = data[year][region][0]
            errs  = data[year][region][1]
            files.sort()
            errs.sort()
            for i,f in enumerate(files):
                tf = ROOT.TFile.Open(f)
                h = tf.Get(hname)
                #print(h.Integral())
                hist.Add(h)

                #calc errs
                df = pd.read_pickle(errs[i])
                sqrddf = df**2
                datadfs.append(sqrddf)

        uncsqddf = sum(datadfs)
        uncdf   = uncsqddf**(1/2)
        
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
                binerr = uncdf[hname][ibin-1]
                hist.SetBinError(ibin,binerr)

        return hist

    def getAddedHistPoissonErrors(self,hist,region,hname,years = [16,17,18]):
        data = self.data
        datadfs = []
        alpha = 1.- 0.682689492
        for year in years:
            files = data[year][region][0]
            errs  = data[year][region][1]
            files.sort()
            errs.sort()
            for i,f in enumerate(files):
                tf = ROOT.TFile.Open(f)
                h = tf.Get(hname)
                #print(h.Integral())
                hist.Add(h)
        
        for ibin in range(hist.GetNbinsX()+1):
            if ibin == 0:
                continue
            else:
                binerr = ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin))+1,1)-hist.GetBinContent(ibin)
                #binerrdn = hist.GetBinContent(ibin) - ROOT.Math.gamma_quantile_c(alpha/2,int(hist.GetBinContent(ibin)),1)
                #if int(hist.GetBinContent(ibin)) == 0:
                #    binerrdn = 0
                
                hist.SetBinError(ibin,binerr)

        return hist
