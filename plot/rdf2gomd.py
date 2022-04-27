#!/usr/bin/env python
import sys

nm=sys.argv[1]

fin=open(nm,"r")

mat=[[],[]]

fins=[]
for i in sys.argv[2:]:
    fins.append(open(i,"r"))
    mat.append([])

read=0

#what matters is that, by when we finish reading the file, we'll have mat[0] contain the fields from the second to last line (the
#x-axis),while mat[1] contains those of the last line (y-axis). For the other files, we'll have only the last line, as all 
#x-axis should be the same.
for i in fin:
    if read>0:
        mat[0]=mat[1]
    mat[1]=i.split()
    read+=1
    if read%2==0:
        mat[1]=mat[1][1:] #the columns for the Y axis begin with the frame number, which we remove here.
    for j,v in enumerate(fins):
        mat[2+j]=v.readline().split()[1:]



for i,v in enumerate(mat[0]):
    t=""
    for j,w in enumerate(mat):
        t="%s %s"%(t,w[i])
    print(t)

