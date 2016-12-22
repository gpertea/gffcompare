#!/bin/sh
ver=$(fgrep '#define VERSION ' gffcompare.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=gffcompare-$ver
linpack=$pack.OSX_x86_64
echo "preparing $linpack.tar.gz"
echo "-------------------"
/bin/rm -rf $linpack
/bin/rm -f $linpack.tar.gz
mkdir $linpack
make clean
make release
cp LICENSE gffcompare $linpack/
tar cvfz $linpack.tar.gz $linpack
ls -l $linpack.tar.gz
echo "scp $linpack.tar.gz  salz:~/html/software/stringtie/dl/"
