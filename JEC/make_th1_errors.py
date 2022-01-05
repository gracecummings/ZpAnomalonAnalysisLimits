import numpy as np
import ROOT
from array import array

infile = "Autumn18_V19_MC_Uncertainty_AK8PFPuppi.txt"

jerar = np.loadtxt(infile,skiprows=1,dtype='double')
etabinlowed = array('d',jerar[:,0])
ptbinlowed  = array('d',jerar[0,3::3])
h18jecunc = ROOT.TH2D("h18jecunc","2018 JEC Uncertainties",len(ptbinlowed)-1,ptbinlowed,len(etabinlowed)-1,etabinlowed)

print(jerar)

#for row in 



#with open(infile) as f:
    #for line in f.readlines():
    
        
