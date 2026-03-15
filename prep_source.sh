#!/usr/bin/env bash
ver=$(fgrep '#define VERSION ' gffcompare.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=gffcompare-$ver
echo " preparing source $pack.tar.gz"
echo "------------------------------------"
/bin/rm -rf $pack $pack.tar.gz
mkdir -p $pack/gclib
cp Makefile LICENSE README.md gffcompare.cpp gtf_tracking.{h,cpp} trmap.cpp $pack/
GCL=./gclib
if [ ! -f $GCL/GBase.h ]; then
  if [ -d .git ]; then
    git submodule sync -- gclib
    git submodule update --init --checkout $GCL
  else
    echo "Error: $GCL/GBase.h not found"
    exit 1
  fi
fi
cp -p $GCL/{GVec,GList,GIntervalTree,GHashMap,khashl}.hh $GCL/{xxhash,wyhash,GBitVec}.h $pack/gclib/
cp -p $GCL/{GArgs,GBase,gdna,GStr,gff,codons,GFaSeqGet,GFastaIndex}.{h,cpp} $pack/gclib/
tar cvfz $pack.tar.gz $pack
ls -l $pack.tar.gz
