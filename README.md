## GffCompare
* compare and evaluate the accuracy of RNA-Seq transcript assemblers (Cufflinks, Stringtie). 
* collapse (merge) duplicate transcripts from multiple GTF/GFF3 files (e.g. resulted from assembly of different samples)
* classify transcripts from one or multiple GTF/GFF3 files as they relate to reference transcripts provided in a
annotation file (also in GTF/GFF3 format)

The original form of this program is also distributed as part of the Cufflinks suite, under the name "CuffCompare" 
(see manual: http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/). All the options and parameters of CuffCompare
are supported by GffCompare (for now), while new features will likely be added to GffCompare in the future.

This package only depends on my other code library, [GCLib](../../../gclib). In order to build GffCompare from this source package 
the following steps can be taken:
```
  cd /some/build/dir
  git clone https://github.com/gpertea/gclib
  git clone https://github.com/gpertea/gffcompare
  cd gffcompare
  make
```
This should build the **gffcompare** binary in the current directory.

