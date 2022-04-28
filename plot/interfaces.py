#!/usr/bin/env python

#This is on the verge of what is the acceptable size for a "simple" script.
#More functionality, and I'll have to refactor the code into functions and whatnot.

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

p = argparse.ArgumentParser()

p.add_argument("fname", type=str, help="Input file, gomd format")
p.add_argument("totalres", type=int, help="Total number of residues in the protein")
p.add_argument("--co", type=int, help="Chain offset, the number of residues in the first chain, if the structure is a dimer and the residue numbers start from 1 in each chain. Note that trimers and above are not supported unless the residue numers are not restarted on each chain.", default="0")
p.add_argument("--forceyrange", type=str,help="Give two numbers separated by a space. Forces the boundary of the y axis to be twose two numbers",default="")
a = p.parse_args()


x=[]
y=[]
for i in range(1,a.totalres+1):
    x.append(i)
    y.append(0)




yrange=[0,1]
force=a.forceyrange
if force!="":
    f=force.split()
    yrange.append(float(f[0]))
    yrange.append(float(f[1]))


#The following is not very pretty, sorry about that.
prop="Frequencia in Interface"
#if you don't want to tag the different numbers, just give a "" in this option.


extension=a.fname.split(".")[-1]
fin=open(a.fname,"r")


frames=0
for line in fin:
        if line.startswith("@") or line.startswith("&") or line.startswith("#"):
            continue
        frames+=1
        fields=line.split()
        reslist=[]
        for i,v in enumerate(fields[1:]):
            w=int(float(v))
            #assign negative residue IDs to the next monomer, if any.
            if w<0:
                w=(-1*w)+a.co
         #   print(len(y), w-1) #############33
            y[w-1]+=1 #If the residue in the position v-1 (i.e. the residue with ID v) is in the list, we add 1 to its position.

for i,v in enumerate(y):
    y[i]=y[i]/frames  #average fraction of its participation on the interface


print("# Interface file in 'standard' goChem format")
for i,v in enumerate(x):
    print("%3d  %3.1f"%(v,y[i]))
         

glyphs=["b-","r-","g-","k-","c-","m-","k--","b^-","ro-","g.-","c:"]        



# changing the rc parameters and plotting a line plot
plt.rcParams['figure.figsize'] = [10, 5]
# Create plots with pre-defined labels.
fig, ax = plt.subplots()



if yrange:
    axes = plt.gca()
    axes.set_ylim(yrange)

plt.xlabel("Residue number")
aafreq=20
plt.xticks(np.arange(1,len(x),aafreq)) # We want more ticks in this case, as we need to pinpoint each residue
   # plt.rcParams['figure.figsize'] = [8, 4] # we need more space also
plt.ylabel("Fraction of Interface participation")


ax.plot(x,y,glyphs[0])
#ax.legend(loc='upper right')
plt.savefig(a.fname.replace("."+extension,".png"),dpi=600)

plt.show()


