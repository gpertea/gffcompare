#!/usr/bin/env bash
set -e
ver=$(fgrep '#define VERSION ' gffcompare.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
srcpack=gffcompare-$ver
source prep_source.sh
linpack=$pack.Linux_x86_64
echo "preparing $linpack.tar.gz"
echo "-------------------"
/bin/rm -rf $linpack
/bin/rm -f $linpack.tar.gz
mkdir $linpack
cd $srcpack
make clean
make static
cp LICENSE README.md gffcompare trmap ../$linpack/
cd ..
tar cvfz $linpack.tar.gz $linpack
ls -l $srcpack.tar.gz $linpack.tar.gz
echo "scp $linpack.tar.gz $srcpack.tar.gz  salz:~/html/software/stringtie/dl/"
echo "(then on the server:)"
echo "perl -i -pe 's/gffcompare\-\d+\.\d+\.\d+\w?\./gffcompare-$ver./g' ~/html/software/stringtie/gff*.shtml"

