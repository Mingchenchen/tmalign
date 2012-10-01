#!/bin/bash
# usage: tmalign-set.sh /path/to/pdblist.dat /path/to/refpdb.pdb

if [[ (! -e $1) || (! -e $2) ]]; then
    echo "usage: $0 /path/to/pdblist.dat /abs/path/to/file.pdb"
    exit
fi

mkdir frag tm output 
pushd frag
baseset=`basename $1`
date > ../output/$baseset.log;
time for f in `cat $1`; do
    echo $f >> ../output/$baseset.log;
    ~ijstokes/projects/tmalign/bin/TMalign-mac64 $GRID_SE/vo/sbgrid/$f $2 > ../tm/`basename $f`.tm 2>> ../output/$baseset.log;
done >> ../output/$baseset.log 2>&1;
date >> ../output/$baseset.log;
popd
