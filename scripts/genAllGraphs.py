import os
import subprocess as sp

uRet = wRet = 0

files = os.listdir( "../mps" )
for f in files:
    print "Generating graph %s..." % f,

    #uCmd = ["../bin/MPS2Graph", "../mps/" + f, "../ungraph/" + f.split(".")[0], "0"]
    wCmd = ["../bin/MPS2Graph", "../mps/" + f, "../graph/" + f.split(".")[0]]

    #uRet = sp.call(uCmd)
    wRet = sp.call(wCmd)

    #if uRet != 0:
    #    print "\nError generating unweighted graph for %s" % f
    if wRet != 0:
        print "\nError generating weighted graph for %s" % f

    print "Done"
