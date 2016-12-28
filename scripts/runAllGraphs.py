import os
import subprocess as sp

files = os.listdir( "../graph" )
for f in files:
    if f == "cap6000_c":
        print "Skipping cap6000_c"
        continue
    if f == "simple":
        print "Skipping simple"
        continue
    #if f == "qp2":
    #    print "Skipping qp2"
    #    continue

    print "Running graph %s..." % f,

    #saucyCmd = ["../bin/saucy", "-t 600", "-r 50", "-s", "../ungraph/" + f]
    saucyCmd = ["../bin/saucy", "-r 1000", "-s", "../ungraph/" + f]
    #saucy0Cmd = ["../bin/saucy0", "-t 600", "-r 50", "-s", "../graph/" + f]
    #saucy2_5Cmd = ["../bin/saucy2-5", "-t 600", "-r 50", "-s", "../graph/" + f]
    saucy2_5Cmd = ["../bin/saucy2-5", "-r 1000", "-s", "../graph/" + f]

    saucyPipe = sp.Popen(saucyCmd, stdout=sp.PIPE, stderr=sp.PIPE)
    #saucy0Pipe = sp.Popen(saucy0Cmd, stdout=sp.PIPE, stderr=sp.PIPE)
    saucy2_5Pipe = sp.Popen(saucy2_5Cmd, stdout=sp.PIPE, stderr=sp.PIPE)

    temp1 = saucyPipe.communicate()
    #temp2 = saucy0Pipe.communicate()
    temp3 = saucy2_5Pipe.communicate()

    if temp1[1]:
        print "\nError in saucy for %s:" % f
        print temp1[1]
    else:
        saucyOut = open( "../res/"+f+".saucy", "w" )
        saucyOut.write( temp1[0] )
        saucyOut.close()

    #if temp2[1]:
    #    print "\nError in saucy0 for %s:" % f
    #    print temp1[1]
    #else:
    #    saucy0Out = open( "../res/"+f+".saucy0", "w" )
    #    saucy0Out.write( temp2[0] )
    #    saucy0Out.close()

    if temp3[1]:
        print "\nError in saucy2-5 for %s:" % f
        print temp1[1]
    else:
        saucy2_5Out = open( "../res/"+f+".saucy2-5", "w" )
        saucy2_5Out.write( temp3[0] )
        saucy2_5Out.close()

    print "Done"
