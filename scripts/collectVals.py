import os
import subprocess as sp

uRet = wRet = 0

files = os.listdir( "../graph" )
for f in files:
    if f == "cap6000_c":
        continue
    if f == "simple":
        continue
    #print "Collecting values for graph %s..." % f,

    saucyFile = open( "../res/"+f+".saucy", "r" )
    saucy2_5File = open( "../res/"+f+".saucy2-5", "r" )

    for line in saucyFile:
        temp = line.split()
        if temp[0] == "vertices":
            sVerts = temp[2]
        elif temp[0] == "edges":
            sEdges = temp[2]
    
    for line in saucy2_5File:
        temp = line.split()
        if temp[0] == "vertices":
            s2Verts = temp[2]
        elif temp[0] == "edges":
            s2Edges = temp[2]
    
    saucyFile.close()
    saucy2_5File.close()

    print "\\verb " + f + " & " + s2Verts + " & " + s2Edges + " & " + sVerts + " & " + sEdges + "\\\\"

    #print "Done"
