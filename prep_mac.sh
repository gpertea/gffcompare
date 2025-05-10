#!/bin/bash
ver=$(fgrep '#define VERSION ' gffcompare.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=gffcompare-$ver
macpack=$pack.OSX_x86_64
echo "preparing $macpack.tar.gz"
echo "-------------------"
/bin/rm -rf $macpack
/bin/rm -f $macpack.tar.gz
mkdir $macpack
make clean
make release
cp LICENSE README.md gffcompare trmap $macpack/
tar cvfz $macpack.tar.gz $macpack
ls -l $macpack.tar.gz
echo "scp $macpack.tar.gz  salz:~/html/software/stringtie/dl/"
