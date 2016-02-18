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
<<<<<<< HEAD

cp -p gffcompare.cpp gtf_tracking.{h,cpp} $pack/
sed 's|\.\./gclib|./gclib|' Makefile > $pack/Makefile
cp -p ../gclib/{GVec,GList,GHash}.hh $libdir
cp -p ../gclib/{GArgs,GBase,gdna,GStr,gff,codons,GFaSeqGet,GFastaIndex}.{h,cpp} $libdir
=======
cp LICENSE gffcompare.cpp gtf_tracking.{h,cpp} $pack/
sed 's|\.\./gclib|./gclib|' Makefile > $pack/Makefile
cp ../gclib/{GVec,GList,GHash}.hh $libdir
cp ../gclib/{GArgs,GBase,gdna,GStr,gff,codons,GFaSeqGet,GFastaIndex}.{h,cpp} $libdir
>>>>>>> 9c4c98b8b220d400843ab3083e92e7e7c9c5babd
tar cvfz $pack.tar.gz $pack
ls -l $pack.tar.gz

