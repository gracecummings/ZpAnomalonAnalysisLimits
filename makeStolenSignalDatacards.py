import ROOT
import gecorg_test as go
import glob
import numpy as np

if __name__=='__main__':
    #Years
    years   = [16,17,18]
    yearstr = go.yearFormatter(years)

    #Starting parameters
    zptcut  = '100.0'
    hptcut  = '300.0'
    metcut  = '75.0'
    btagwp  = '0.8'
    strdescrip = 'unblindttshenanigans2stolen'
    chan = 'mumu'

    template_card = 'test_datacards_2024-10-03/Run2_161718_ZllHbbMET_datacard_unblindttshenanigans2_mumu_Zp2000ND400NS200_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.txt'
    template_root = 'test_datacards_2024-10-03/Run2_161718_ZllHbbMET_unblindttshenanigans_withsysts_mumu_Zp2000ND400NS200_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.root'

    tosteal_cards = glob.glob('test_datacards_2024-10-03/Run2_161718_ZllHbbMET_datacard_vvuncbugfix_mumu_Zp*.txt')
    tosteal_roots = glob.glob('test_datacards_2024-10-03/Run2_161718_ZllHbbMET_vvuncbugfix_mumu_Zp*.root')
    tosteal_cards.sort()
    tosteal_roots.sort()
    tosteal_tuple = list(zip(tosteal_cards,tosteal_roots))

    #Lines to take from the template datacard
    open_template = open(template_card,"r")
    lines = open_template.readlines()
    systs = lines[15:]#cuts out part manunally written bellow
    gen_systs = [syst for syst in systs if 'Zp2000ND400NS200' not in syst]#removes the signal sample specific systematics
    print("Building template of background systematics for datacard")
    syst_dict = {}
    for syst in gen_systs:
        breaks = syst.split()
        systname = breaks[0]+"$"+breaks[1]
        syst_dict[systname] = {"sig":breaks[-4],"DY":breaks[-3],"TT":breaks[-2],"VV":breaks[-1],"type":breaks[-5]}

    #Histograms to take from the template root file
    temp_root = ROOT.TFile.Open(template_root,'read')
    keys = [key.GetName() for key in list(temp_root.GetListOfKeys())]
    bkg_keys = [key for key in keys if 'Zp2000ND400NS200' not in key]#these are the histograms that need to be added
    bkg_hists = []
    for item in bkg_keys:
        bkg_hists.append(temp_root.Get(item))
        bkg_hists[-1].SetDirectory(0)

    #Beginning stealing
    for tup in tosteal_tuple:
        card = tup[0]
        root = tup[1]
        protoname = card.replace("test_datacards_2024-10-03/Run2_161718_ZllHbbMET_datacard_vvuncbugfix_mumu_","")
        name = protoname.replace("_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.txt","")
        
        print("Gathering info for ",name)
        #Steal parts of root file
        mark_root = ROOT.TFile.Open(root,'read')
        keys = [key.GetName() for key in list(mark_root.GetListOfKeys())]
        sig_keys = [key for key in keys if name in key]#these are the histograms that need to be added
        hsig = mark_root.Get(name)
        sig_yield = str(hsig.Integral())
        sig_hists = []
        for item in sig_keys:
            sig_hists.append(mark_root.Get(item))
            sig_hists[-1].SetDirectory(0)

        #Steal parts of datacard
        open_mark = open(card,"r")
        lines = open_mark.readlines()
        systs = lines[15:]#cuts out part manunally written bellow
        sig_systs = [syst for syst in systs if name in syst]#the signal sample specific systematics
        check_systs = list(set(systs)-set(sig_systs))
        check_systs = [syst for syst in check_systs if "TT_StatsUnc" not in syst]#remove the old ttbar ones
        
        #Update the lines that need to be updated
        print("    Updating signal systemaitcs lines for the datacard")
        for syst in check_systs:
            pieces = syst.split()
            check_name = pieces[0]+"$"+pieces[1]

            if "#" in check_name:
                continue
            if syst_dict[check_name]["type"]== 'lnN' and "-" not in syst_dict[check_name]["sig"]:
                #print("This is the line with signal values we need:")
                #print("    ",syst)
                #print("    value extracted: ",pieces[-4])
                #print("    previous value:  ",syst_dict[check_name]["sig"])
                syst_dict[check_name]["sig"] = pieces[-4]
                #print("    resulting line info: ",syst_dict[check_name]["sig"])

        print("    Writing the new root file")
        hists_to_write = bkg_hists+sig_hists
        hists_to_write.sort()
        rootFileName = go.makeOutFile('Run2_161718_ZllHbbMET',strdescrip+'_'+chan+'_'+name,'.root',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        outf = ROOT.TFile(rootFileName,"recreate")
        for hist in hists_to_write:
            hist.Write()
        outf.Write()
        outf.Close()

        print("    Beginning to write the new datacard")
        #write the the datacard
        prepCardName = go.makeOutFile('Run2_'+yearstr+'_ZllHbbMET','datacard_'+strdescrip+'_'+chan+'_'+name,'.txt',str(zptcut),str(hptcut),str(metcut),str(btagwp))
        datcard = open(prepCardName,"w")

        #Write the card
        datcard.write("imax 1\n")
        datcard.write("jmax 3\n")
        datcard.write("kmax *\n")              
        datcard.write("------------\n")
        datcard.write("shapes * * {0} $PROCESS $PROCESS_$SYSTEMATIC \n".format(rootFileName.split("/")[-1]))
        datcard.write("------------\n")
        datcard.write("bin {0} \n".format(name+"_mumu"))
        datcard.write("observation 26\n")
        datcard.write("------------\n")
        datcard.write("bin {0} {1} {2} {3}\n".format(name+"_mumu",name+"_mumu",name+"_mumu",name+"_mumu"))
        datcard.write("process "+name+" DY TT VV\n")
        datcard.write("process 0 1 2 3\n")
        datcard.write("rate "+sig_yield+" 6.51803520321846 7.19403937458992 0.7364198043942451\n")
        datcard.write("------------\n")
        datcard.write("lumi_13TeV                    lnN 1.016 - - 1.016\n")

        print("    Adding the systematics for the new datacard")
        #find the offset needed to make the datacard readable
        syst_descrips = list(set([x.split("$")[0] for x in syst_dict.keys()]))
        maxlength = max([len(x) for x in syst_descrips])
        for syst in syst_dict.keys():
            nameoffset = maxlength - len(syst.split("$")[0])
            cardstr = "{0} {1} {2} {3} {4} {5}\n".format(syst.split("$")[0]+" "*nameoffset,syst_dict[syst]["type"],syst_dict[syst]["sig"],syst_dict[syst]["DY"],syst_dict[syst]["TT"],syst_dict[syst]["VV"])
            datcard.write(cardstr)

        print("    Adding the signal stats uncertainty lines to datacard")
        for syst in sig_systs:
            datcard.write(syst)

        print("    Finished writing datacard")
        print("    Closing the signal sources.")
        mark_root.Close()
