## GffCompare
* compare and evaluate the accuracy of RNA-Seq transcript assemblers (Cufflinks, Stringtie). 
* collapse (merge) duplicate transcripts from multiple GTF/GFF3 files (e.g. resulted from assembly of different samples)
* classify transcripts from one or multiple GTF/GFF3 files as they relate to reference transcripts provided in a
annotation file (also in GTF/GFF3 format)

The official documentation and download packages for this utility can be found online here:
http://ccb.jhu.edu/software/stringtie/gffcompare.shtml

For more information about the GFF3/GTF file formats expected by GffCompare please refer to the available online documentation -- a quick review can be found at: http://ccb.jhu.edu/software/stringtie/gff.shtml

The original version of this program was distributed as part of the Cufflinks suite, under the name "CuffCompare".
The overall functionality and most of the options of CuffCompare are still supported by GffCompare, while new functionality is being added to GffCompare only, as it is the program which is actively maintained.

An example of a new feature of GffCompare (compare to its predecessor CuffCompare) is this: when a single query GTF/GFF file is given as input for analysis, along with a reference annotation (-r option), GffCompare switches into *annotation mode* and it generates a *.annotated.gtf* file instead of the *.combined.gtf* produced by CuffCompare with the same parameters. This file has the same general format as CuffCompare's *.combined.gtf* file (with "class codes" assigned to transcripts as per their relationship with the matching/overlapping reference transcript), but the original transcript IDs are preserved, so GffCompare can thus be used as a simple way of annotating a set of transcripts.

Another important difference is that the input transcripts are by default no longer discarded when they are found to be "intron redundant", i.e. contained within other, longer isoforms. CuffCompare had the -G option to prevent collapsing of such intron redundant isoforms into their longer "containers", but GffCompare has made this the default mode of operation (hence the -G option is no longer needed and is simply ignored when given). However please note that transcripts with fully redundant intron chains (i.e. with the same exact intron coordinates, hence the same intron-exon structure except the terminal exon ends) are *still* discarded when GTF/GFF files are loaded

Steps for building this package from source (the only dependency is my other code library, [GCLib](../../../gclib)):
```
  cd /some/build/dir
  git clone https://github.com/gpertea/gclib
  git clone https://github.com/gpertea/gffcompare
  cd gffcompare
  make release
```
This should build the **gffcompare** binary in the current directory.



