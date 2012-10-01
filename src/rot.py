import sys

def usage():
    print "usage: rot.py matrix-file in-pdb out-pdb"

def rotatepdb(rot_fn,infile,outfile):
    rot_fh = open(rot_fn)
    trans  = map(float,rot_fh.readline().split())
    rot    = []
    rot.append(map(float,rot_fh.readline().split()))
    rot.append(map(float,rot_fh.readline().split()))
    rot.append(map(float,rot_fh.readline().split()))
        
    in_fh  = open(infile)
    out_fh = open(outfile,"w")
    for line in in_fh:
        if (line[0:5] == "ATOM "):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            xnew = trans[0] + rot[0][0]*x + rot[0][1]*y + rot[0][2]*z
            ynew = trans[1] + rot[1][0]*x + rot[1][1]*y + rot[1][2]*z
            znew = trans[2] + rot[2][0]*x + rot[2][1]*y + rot[2][2]*z
            out_fh.write("%s% 8.3f% 8.3f% 8.3f%s" % (line[0:30], xnew, ynew, znew, line[54:]))
        else:
            out_fh.write(line)
    in_fh.close()
    out_fh.close()

if len(sys.argv) != 4:
    usage()
    sys.exit(-1)
else:
    rotatepdb(sys.argv[1],sys.argv[2],sys.argv[3])
