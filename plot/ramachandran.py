#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys


fin=open(sys.argv[1],"r")

phi=[]
psi=[]
frame=[]

for i in fin:
    line=i.split()
    frame.append(float(line[0]))
    phi.append(float(line[1]))
    psi.append(float(line[2]))

nframe=np.array(frame)/frame[-1]

nphi=np.array(phi)
npsi=np.array(psi)

fig, ax = plt.subplots()

cmap=plt.get_cmap("jet")#"Spectral")



plt.xlabel('Phi (deg)')
plt.ylabel('Psi (deg)')

axes = plt.gca()

axes.set_xlim([-180,180])
axes.set_ylim([-180,180])


plt.scatter(nphi,npsi,c=nframe,marker="+",cmap=cmap) #.reversed())



plt.savefig(sys.argv[1].replace(".dat",".png"))

plt.show()


