#!/usr/bin/env python
import sys

#Usage: Copy the last 2 lines from the goMD rdf output and paste them to a file. Call this script on that file. It will write a plots.py formatted file to stdout, plus one line jump that should not be considered. You can also add the 2 last lines of several goMD rdf outputs (each pair of lines in different files) and give them to this script, which will put them as separate columns in the same plots.py file.

nm=sys.argv[1]

fin=open(nm,"r")

mat=[[],[]]

fins=[]
for i in sys.argv[2:]:
    fins.append(open(i,"r"))
    mat.append([])



#what matters is that, by when we finish reading the file, we'll have mat[0] contain the fields from the second to last line (the
#x-axis),while mat[1] contains those of the last line (y-axis). For the other files, we'll have only the last line, as all 
#x-axis should be the same.
mat[0]=fin.readline().split()
mat[1]=fin.readline().split()[1:]  #the columns for the Y axis begins with the frame number, which we remove here.
fin.close()

for j,v in enumerate(fins):
    v.readline() #This is the x axis, which is info we already have.
    mat[2+j]=v.readline().split()[1:]
    v.close()

print("mat",mat,"end mat") ########################

for i,v in enumerate(mat[0]):
    t=""
    for j,w in enumerate(mat):
        t="%s %s"%(t,w[i])
    print(t)

