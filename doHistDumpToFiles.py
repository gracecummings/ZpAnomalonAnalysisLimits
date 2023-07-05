import glob
import ROOT
import gecorg_test as go

files = glob.glob("analysis_output_ZpAnomalon/2022-09-30/Zp2000ND400NS200_expectedsignal*")
sfiles = [f.split("_Zptcut100_Hptcut300_metcut75_btagwp8E-10.root")[0].split("analysis_output_ZpAnomalon/2022-09-30/")[-1] for f in files]
#limdict = {"0.0":"0","0.7492":"lowexp","1.0859":"medexp","1.6097":"highexp"}#2000400200
#limdict = {"0.0":"0","0.8553":"lowexp","1.2695":"medexp","1.9426":"highexp"}#4000800200
#limdict = {"0.0":"0","1.158":"lowexp","1.7422":"medexp","2.6797":"highexp"}#55001800200
#repstr = "tests_default_seed"

#limdict =  {"0.0":"0","0.5488":"lowexp","0.8086":"medexp","1.2244":"highexp"}#2000400200
#limdict =  {"0.0":"0","0.854":"lowexp","1.2422":"medexp","1.8611":"highexp"}#4000800200
#limdict = {"0.0":"0","1.2192":"lowexp","1.7734":"medexp","2.6641":"highexp"}#55001800200

#limdict = {"0.0":"0","0.0782":"lowexp","0.124":"medexp","0.2061":"highexp"}#55001800200
#limdict =  {"0.0":"0","0.5328":"lowexp","0.7852":"medexp","1.1826":"highexp"}#2000400200
#limdict =  {"":"0","":"lowexp","0.161743":"medexp","":"highexp"}#4000800200

######THE 2800 TO 5000 BIN
#limdict = {"0.0":"0","0.2349":"lowexp","0.375":"medexp","0.647":"highexp"}#55001800200
limdict =  {"0.0":"0","0.5328":"lowexp","0.7852":"medexp","1.1951":"highexp"}#2000400200
#limdict =  {"0.0":"0","0.1857":"lowexp","0.2822":"medexp","0.4532":"highexp"}#4000800200


repstr = "2800to5000bin"
for f in files:
    print(f)
    n1= f.split("_Zptcut100_Hptcut300_metcut75_btagwp8E-10.root")[0].split("analysis_output_ZpAnomalon/2022-09-30/")[-1]
    lim = n1.split("_"+repstr)[0].split("expectedsignal")[-1]
    name = n1.replace(lim+"_"+repstr,limdict[lim])+"_2800to5000bin"
    tf = ROOT.TFile(f)
    keys = tf.GetListOfKeys()
    nkeys = [key.GetName() for key in keys]
    keyset = set(nkeys)
    #badkeys = set(["hdrerr","hdrerrge","hrge","hrge","hdrge","hdrerrgehl"])
    badkeys = set(["hdrerr","hdrerrge","hrge"])
    names = {"hdr":"diff","hdrerrhl":"pull","hdrerrgehl":"gepull","hdrge":"gediff","hr":"rdist"}
    pltkeys = keyset-badkeys
    for key in pltkeys:
        hname = key
        h = tf.Get(hname)
        tc = ROOT.TCanvas("tc","tc",1060,800)
        tc.Draw()
        tc.cd()
        h.Draw()
        tc.cd()

        outname = go.makeOutFile(name,names[hname],".png","100","300","75","8E-10")
        #print(outname)
        tc.SaveAs(outname)
        
        
