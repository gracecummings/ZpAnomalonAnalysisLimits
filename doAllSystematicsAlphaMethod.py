import configparser
import subprocess

if __name__=='__main__':


    config = configparser.RawConfigParser()
    config.optionxform = str
    fp = open('systematics.ini')
    config.read_file(fp)
    systs = config.sections()
    #print("this is a test")
    print(systs[1:])

    for syst in systs[1:]:
        if "btag" not in syst:
            continue
        print("Doing DY estimation with "+syst+"up")
        subprocess.run(["python","doAlphaRatioFits_systematics.py","-m", "75.0", "-z", "100.0","-j","300.0","-wp","0.8","-sn",syst,"-sd","up"])
        print("Doing DY estimation with "+syst+"dwn")
        subprocess.run(["python","doAlphaRatioFits_systematics.py","-m", "75.0", "-z", "100.0","-j","300.0","-wp","0.8","-sn",syst,"-sd","dwn"])
