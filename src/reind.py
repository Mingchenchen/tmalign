import re

class aa:

    def __init__(self,name, ab3, ab1):
        self.name = str(name).upper()
        self.ab3  = str(ab3).upper()
        self.ab1  = str(ab1).upper()

aas = ( aa("Alanine", "ALA", "A"), \
aa("Cysteine",	"CYS", "C"), \
aa("AsparticAcid",	"ASP", "D"), \
aa("GlutamicAcid",	"GLU", "E"), \
aa("Phynylalanine",	"PHE", "F"), \
aa("Glycine",	"GLY", "G"), \
aa("Histidine",	"HIS", "H"), \
aa("Isoleucine",	"ILE", "I"), \
aa("Lysine",	"LYS", "K"), \
aa("Leucine",	"LEU", "L"), \
aa("Methionine",	"MET", "M"), \
aa("Asparagine",	"ASN", "N"), \
aa("Proline",	"PRO", "P"), \
aa("Glutamine",	"GLN", "Q"), \
aa("Arginine",	"ARG", "R"), \
aa("Serine",	"SER", "S"), \
aa("Threonine",	"THR", "T"), \
aa("Valine",	"VAL", "V"), \
aa("Tryptophan",	"TRP", "W"), \
aa("Tyrosine",	"TYR", "Y")
)
#aa("NO MATCH",	"XXX", "-"), \

ab1to3 = {}
for acid in aas:
  ab1to3[acid.ab1] = acid.ab3

def newind(base,output="reindex.pdb"):
    in_fh    = open(base + ".pdb",'r')
    seq_fh   = open(base + ".seq",'r')
    out_fh   = open(output,'w')
    debug_fh = open(base + "-debug.pdb",'w')
    atom     = 0
    res      = 0
    reslast  = 0
    pdb_struct = {}
    res_struct = []

    seq_ind  = 0

    in_sync = True

    for line in in_fh.readlines():
        if (not line.startswith('ATOM ')):
            continue
        resnew   = line[22:27]
        restype  = line[17:20]
        atomtype = line[21:22]

        if (resnew != reslast): # We've started a new residue in the PDB file
            if (in_sync): # if we were previously in sync
                # output the last residue and increment the residue counter
                pdb_struct[res] = res_struct
                res_struct = []
                res += 1
                # Find the new expected sequence from the FASTA alignment file
                # output the last residue and increment the residue counter
                seq = seq_fh.readline() 
                seq_ind  += 1
                while not re.match("[A-Za-z]",seq[0]):
                    seq = seq_fh.readline()

            in_sync = False # Assume that we are now NOT in sync

            # check to see if this new residue matches the FASTA sequence residue
            if (restype == ab1to3[seq[:-1].upper()]): # then our current residue matches the FASTA residue
                in_sync = True
            else: # Not in sync with residues yet, so only update last residue index and continue
                reslast = resnew
                continue

        atom += 1
        res_struct.append(line[0:5] + "%6i" % (atom) + line[11:22] + "%4i" % (res) + line[27:-1] + " " + line[22:27])
        out_fh.write(line[0:5] + "%6i" % (atom) + line[11:22] + "%4i" % (res) + line[27:-1])
        #assert(res == seq_ind)
        #assert(restype == ab1to3[seq[:-1].upper()])
        debug_fh.write(line[0:5] + "%6i" % (atom) + line[11:22] + "%4i" % (res) + line[27:-1] + " " + line[22:27] + " " + str(seq_ind) + " " +  ab1to3[seq[:-1].upper()] + "\n")
        reslast = resnew

    out_fh.write("TER  %6i      %s %s%4i" % (atom, restype, atomtype, res) + "\n")
    out_fh.write("END\n")

    out_fh.close()
    debug_fh.close()


    return pdb_struct

def getres(input, code, output="residues.pdb"):
    in_fh   = open(input,'r')
    pdb_fh  = open(code.upper() + ".pdb",'r')
    out_fh  = open(output,'w')

    pdb_struct = newind(code.upper() + ".pdb", code.upper() + "-reind.pdb")

    in_tgt  = False
    seq     = ''
    for line in in_fh.readlines():
        if (not in_tgt):
            if (line.startswith('>t|PDB:' + code.lower())):
                in_tgt = True
            continue
        if (len(line) == 1):
            break
        seq += line[:-1]

    residx = 0
    for c in list(seq):
        print residx, c
        if (c == "-"):
            continue

        residx += 1
        if re.match("[A-Z]",c):
            for a in pdb_struct[residx]:
                print a

    print len(seq), 
    print seq
