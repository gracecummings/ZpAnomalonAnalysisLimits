import gecorg_test as go
import glob
import ROOT

#config = configparser.RawConfigParser()
#config.optionxform = str
#fp = open('xsects_2017.ini')#2017 and 2018dy+jets xs same
#config.read_file(fp)
#xspairs = config.items('TT')

#ttbarbkgs17 = glob.glob('analysis_output_ZpAnomalon/2021-10-28/Fall17.TTTo*_upout_metcomp*.root')
#ttbarbkgs18 = glob.glob('analysis_output_ZpAnomalon/2021-10-28/Autumn18.TTTo*_upout_metcomp*.root')

inputs = glob.glob('analysis_output_ZpAnomalon/2021-10-28/*_upout_metcomp*.root')

for f in inputs:
    tf = ROOT.TFile(f)
    hpfmet = tf.Get('h_pfmet')
    hmetclean = tf.Get('h_metclean')
    hgenmet = tf.Get('h_genmet')
    hpfmet.SetLineColor(ROOT.kBlue)
    hmetclean.SetLineColor(ROOT.kRed)
    hgenmet.SetLineColor(ROOT.kBlack)

    leg = ROOT.TLegend(0.50,0.50,0.88,0.88)
    leg.AddEntry(hmetclean,"METClean","l")
    leg.AddEntry(hgenmet,"GenMET","l")
    leg.AddEntry(hpfmet,"PFMET","l")
    

    tc = ROOT.TCanvas("tc","tc",700,600)
    tc.cd()
    hpfmet.Draw("HIST")
    hmetclean.Draw("HISTSAME")
    hgenmet.Draw("HISTSAME")
    leg.Draw()

    protoname = f.split("/")[-1]
    name = protoname.split("_upout")[0]

    pngname = go.makeOutFile(name,"metcomp_btagsf",".png","0.0","250.0","0.0","0.0")
    tc.SaveAs(pngname)
