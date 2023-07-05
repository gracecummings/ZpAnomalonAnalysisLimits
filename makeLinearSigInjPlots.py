import ROOT
import matplotlib.pyplot as plt
import numpy as np
import glob

mp = 'Zp2000ND400NS200'
infs = glob.glob("analysis_output_ZpAnomalon/2022-09-30/"+mp+"*_expectedsignal*signalinjectionplots_Zptcut100_Hptcut300_metcut75_btagwp8E-10.root")

pdict = {}
injs = []

for f in infs:
    inj = f.split("_expectedsignal")[-1].split("_")[0]
    pdict[float(inj)] = f
    injs.append(float(inj))

injs = sorted(injs)
rmeas = []
errs  = []

for p in injs:
    print('Injected signal: ',p)
    tf = ROOT.TFile(pdict[p])
    h = tf.Get("hrge")
    f1 = h.GetFunction("f1")
    mu = f1.GetParameter(1)
    sigma = f1.GetParameter(2)
    print("The mean of the guassian: ",mu)
    print('The sigma of the gaussian: ',sigma)
    rmeas.append(mu)
    errs.append(sigma)

x = np.linspace(-.2,round(max(injs)*2,1),5)
    
fig = plt.figure()

plt.plot(x,x,label="y=x")
plt.errorbar(injs,rmeas,yerr=errs,label="r v. r injected, 1 sigma bars",fmt=".k")
plt.legend()
plt.title(mp)
plt.xlabel("signal strength (r), injected")
plt.ylabel("signal strenght (r), measured")
plt.savefig(mp+"_rvrinj_emu_ge_extendedlimit.png")
plt.savefig(mp+"_rvrinj_emu_ge_extendedlimit.pdf")
plt.show()


