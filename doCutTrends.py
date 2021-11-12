import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("cutgrid_corrweights_1111.txt")
zmap = np.zeros((5,2))

#row structure
#Zptcut metcut zzyieldsb wzyieldsb ttyieldsb dyyieldsb totyieldsb

valcheck = 11
xbinedges = np.zeros(5)
ybinedges = np.zeros(2)

titledict = {
    2:"ZZ sideband yield",
    3:"WZ sideband yield",
    4:"TTbar sideband yield",
    5:"DY sideband yield",
    6:"total background in sideband",
    7:"ZZ signal region yield",
    8:"WZ signal region yield",
    9:"TTbar signal region yield",
    10:"DY signal region yield",
    11:"total background in signal region"
    }

for row in data:
    #print("Zptcut metcut zzyieldsb wzyieldsb ttyieldsb dyyieldsb totyieldsb")
    #print(row)
    print(" zpt: ",row[0])
    
    print(" met: ",row[1])
    print(" int: ",row[valcheck])
    x = int((row[0]-100)/25.0)
    y = int((row[1]-50)/25)
    z = row[valcheck]

    print(" xbin: ",x)
    print(" xbin: ",y)
    zmap[x,y] = z
    xbinedges[x] = row[0]
    ybinedges[y] = row[1]

zerosbotleft = np.flip(zmap.T,axis=0)#zpt along horizontal,met along vertical, (0,0) in bottom left corner

fig, ax = plt.subplots()

print(zerosbotleft)

#for heat maps of cuts
ax.set_title(titledict[valcheck])
c = ax.pcolormesh(xbinedges,ybinedges,zmap.T,shading="auto")#hack for zpt edges (x-axis)
#c = ax.pcolormesh(zmap.T,shading="auto")#working with stupid labels
ax.set_ylabel('met cut')
ax.set_xlabel('z p_{T} cut')
fig.colorbar(c,ax=ax)
plt.show()

#for running plots
