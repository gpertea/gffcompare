#!/bin/sh
ver=$(fgrep '#define VERSION ' gffcompare.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=gffcompare-$ver
echo "preparing $pack.tar.gz"
echo "----------------------"
/bin/rm -rf $pack
/bin/rm -f $pack.tar.gz
mkdir $pack
mkdir $pack/gclib
libdir=$pack/gclib/

cp -p Makefile gffcompare.cpp gtf_tracking.{h,cpp} $pack/
cp -p ../gclib/{GVec,GList,GHash}.hh $libdir
cp -p ../gclib/{GArgs,GBase,gdna,GStr,gff,codons,GFaSeqGet,GFastaIndex}.{h,cpp} $libdir
tar cvfz $pack.tar.gz $pack
ls -l $pack.tar.gz

