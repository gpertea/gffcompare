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

cp LICENSE Makefile gffcompare.cpp gtf_tracking.{h,cpp} $pack/
cp ../gclib/{GVec,GList,GHash}.hh $libdir
cp ../gclib/{GArgs,GBase,gdna,GStr,gff,codons,GFaSeqGet,GFastaIndex}.{h,cpp} $libdir
tar cvfz $pack.tar.gz $pack
ls -l $pack.tar.gz

