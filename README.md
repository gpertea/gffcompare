## GffCompare
* compare and evaluate the accuracy of RNA-Seq transcript assemblers (Cufflinks, Stringtie). 
* collapse (merge) duplicate transcripts from multiple GTF/GFF3 files (e.g. resulted from assembly of different samples)
* classify transcripts from one or multiple GTF/GFF3 files as they relate to reference transcripts provided in a
annotation file (also in GTF/GFF3 format)

The original form of this program is also distributed as part of the Cufflinks suite, under the name "CuffCompare" 
(see manual: http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/). Most of the options and parameters of CuffCompare
are supported by GffCompare, while new features will likely be added to GffCompare in the future.

A notable difference from GffCompare is that when a single query GTF/GFF file is given as input, along with a reference annotation (-r option),
gffcompare switches into "annotation mode" and it generates a .annotated.gtf file instead of the .combined.gtf produced by CuffCompare with the 
same parameters. This file has the same general format as CuffCompare's .combined.gtf file (with "class codes" assigned to transcripts as per 
their relationship with the matching/overlapping reference transcript),  but the original transcript IDs are preserved, so gffcompare can thus 
be used as a simple way of annotating a set of transcripts.

Another important difference is that the input transcripts are no longer discarded when they are found to be "intron redundant", i.e. 
contained within other, longer isoforms. CuffCompare had the -G option to prevent collapsing of such intron redundant isoforms into 
their longer "containers", but GffCompare has made this the default mode of operation (hence the -G option is no longer needed 
and is simply ignored when given).


Steps for building this package from source (the only dependency is my other code library, [GCLib](../../../gclib)):
```
  cd /some/build/dir
  git clone https://github.com/gpertea/gclib
  git clone https://github.com/gpertea/gffcompare
  cd gffcompare
  make release
```
This should build the **gffcompare** binary in the current directory.



