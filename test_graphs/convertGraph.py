import sys

if __name__ == "__main__":
    f = open( sys.argv[1], "r" )
    edges = []
    for line in f:
        data = line.split()
        if data[0] == "c":
            continue
        elif data[0] == "p":
            V = int(data[2])
            E = int(data[3])
        elif data[0] == "e":
            temp1 = (int(data[1])-1, int(data[2])-1)
            temp2 = (int(data[2])-1, int(data[1])-1)
            if temp1 not in edges and temp2 not in edges:
                edges.append( (int(data[1])-1, int(data[2])-1) )
    f.close()

    f = open( sys.argv[1]+".amorph", "w" )
    f.write( str(V)+" "+str(len(edges))+" 1\n" )
    for i in xrange( V ):
        f.write( str(i)+" 0\n" )
    for e in edges:
        f.write( str(e[0])+" "+str(e[1])+" 1\n")
    f.close()
