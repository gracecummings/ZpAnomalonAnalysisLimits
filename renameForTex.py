import glob
import subprocess

infiles = glob.glob("mumu_2022-03-31_ProperREOIDSF/h*")
#print(infiles)

for f in infiles:
    protoname = f.split("Zptcut")[0]
    fsplit = f.split("_")
    cuts = [piece.split("cut") for piece in fsplit[-4:]]
    zpt  = cuts[0][-1]
    hpt  = cuts[1][-1]
    met  = cuts[2][-1]
    btagstr = cuts[3][0]
    wp    = btagstr.split(".")[1]
    wpstr = str(wp)+"E-10"
    texname = protoname+"Zptcut"+zpt[:-2]+"_Hptcut"+hpt[:-2]+"_metcut"+met[:-2]+"_btagwp"+wpstr+".png"
    subprocess.run(["mv",f,texname])
