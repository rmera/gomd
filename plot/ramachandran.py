#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

### How to use:

#ramachandran.py rama.dat [-RGB] [-trajcolor] [alpha]

#If the corresponding goMD run was performed with the rgb option you need to give the -RGB option here.
#In that case you can also give the -trajcolor to ignore those rgb values and color the "points" by 
#frame. Note that the options are cap-sensitive! "-rgb" or "-TrajColor" won't work.

#If the goMD run was _not_ performed with the rgb option, do not give either the -RGB or -trajcolor
#options. The angles will be colored by frame.

#In either case you can give a number between 0 and 1 for alpha or the "opaqueness" of each point in the plot.
#(i.e. lower values make the points more transparent). Use lower values if you have lots of points, so the last 
#ones don't cover the first ones.

#Note that all options must be given _after_ the rama.dat (i.e. the output from goMD).


fin=open(sys.argv[1],"r")

alpha=0.5
rgb=False
trajcolor=False
if "-RGB" in sys.argv: #note the caps sensitivity
    rgb=True
    del(sys.argv[sys.argv.index("-RGB")]) #so it's not a problem when retrieving the alpha

if "-trajcolor" in sys.argv: #note the caps sensitivity
    trajcolor=True
    del(sys.argv[sys.argv.index("-trajcolor")]) #so it's not a problem when retrieving the alpha


if len(sys.argv)>2:
    alpha=float(sys.argv[2])

phi=[]
psi=[]
rgbs=[]
frame=[]

first=True
angles=0
for i in fin:

    line=i.split()
    if first:
        n=0
        if rgb:
            n=5
        else:
            n=2
        angles=(len(line)-1)//n
        skip=n
        for j in range(angles):
            phi.append([])
            psi.append([])
            if rgb:
                rgbs.append([]) 
        first=False
    frame.append(float(line[0]))
    ph=[]
    ps=[]
    for j in range(angles):
        phi[j].append(float(line[1+j*skip]))
        psi[j].append(float(line[2+j*skip]))
        if rgb:
            rgbs[j].append([(float(line[3+j*skip]))/255,(float(line[4+j*skip]))/255,(float(line[5+j*skip]))/255])

nframe=np.array(frame)/frame[-1]

fig, ax = plt.subplots()

cmap=plt.get_cmap("jet")#"Spectral")



plt.xlabel('Phi (deg)')
plt.ylabel('Psi (deg)')

axes = plt.gca()

axes.set_xlim([-180,180])
axes.set_ylim([-180,180])

glyphs=["+","2","x","|","-","1"]
nglyphs=len(glyphs)
if len(phi)>nglyphs:
    glyphs=[]
    #If there are too many residues we just give up on the different glyphs
for i in range(len(phi)):
        glyphs.append("+")

for i,ph in enumerate(phi):
        if not rgb or trajcolor:
            plt.scatter(ph,psi[i],c=nframe,marker=glyphs[i],alpha=alpha,cmap=cmap) #.reversed())
        else:
            plt.scatter(ph,psi[i],c=rgbs[i],marker=glyphs[i],alpha=alpha) #.reversed())



plt.savefig(sys.argv[1].replace(".dat",".png"))

plt.show()


